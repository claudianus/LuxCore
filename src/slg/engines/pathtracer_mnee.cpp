/***************************************************************************
 * Copyright 1998-2026 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxCoreRender.                                   *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software      *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and      *
 * limitations under the License.                                           *
 ***************************************************************************/

// MNEE (Manifold Next Event Estimation): direct light sampling through a
// single delta specular chain x0 -> x1 -> y.
//
// Reference: Hanika, Droske, Fascione 2015 (MNEE) and Zeltner, Hanika, Gross
// 2020 (Specular Manifold Sampling, SS solver). The implementation follows the
// official SMS reference (dev-tools/reference/manifold_ss.cpp) with two
// adaptations:
//   - the constraint Jacobian dH/dX is computed numerically (finite
//     differences with a re-projection trace) instead of analytically, so
//     curvature (dndu/dndv from the mesh differentials) is included without
//     hand-porting the full derivative algebra;
//   - the specular factor is evaluated with LuxCore's own material code
//     (MirrorMaterial::Kr, GlassMaterial::EvalSpecularTransmission static
//     method) so the estimator is exactly consistent with the renderer's
//     BSDF sampling.
// Validated against a virtual-light truth in dev-tools/mnee_design.md
// (numpy port of the analytic pipeline: E[C_mnee]/E[C_virtual] = 1.00000).
//
// Disjointness (unbiasedness) with the plain direct light estimator: the
// plain estimator contributes 0 whenever a delta specular surface blocks the
// shadow ray (Scene::Intersect stops there), and forward BSDF sampling hits a
// positional delta light with probability 0. MNEE fills exactly the paths
// x0 -> x1 (delta specular) -> y (point/spot/mappoint), so no MIS is needed.
// Mesh/area lights are excluded from MNEE (forward BSDF sampling covers them
// through delta speculars) until a MIS weight is added.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "luxrays/usings.h"
#include "slg/lights/light.h"
#include "slg/usings.h"
#include "slg/engines/pathtracer.h"
#include "slg/materials/mirror.h"
#include "slg/materials/glass.h"
#include "slg/scene/scene.h"

using namespace std;
using namespace luxrays;
using namespace slg;

//------------------------------------------------------------------------------
// Small linear algebra helpers (2x2) and the specular chain vertex
//------------------------------------------------------------------------------

namespace {

struct MneeVec2 {
	float x, y;
};

struct MneeMat2 {
	// [a11 a12]
	// [a21 a22]
	float a11, a12, a21, a22;
};

inline float MneeDet(const MneeMat2 &m) {
	return m.a11 * m.a22 - m.a12 * m.a21;
}

inline MneeMat2 MneeInverse(const MneeMat2 &m, const float det) {
	return MneeMat2{ m.a22 / det, -m.a12 / det, -m.a21 / det, m.a11 / det };
}

inline MneeMat2 MneeMul(const MneeMat2 &A, const MneeMat2 &B) {
	return MneeMat2{
		A.a11 * B.a11 + A.a12 * B.a21, A.a11 * B.a12 + A.a12 * B.a22,
		A.a21 * B.a11 + A.a22 * B.a21, A.a21 * B.a12 + A.a22 * B.a22
	};
}

// A specular chain vertex on the caustic caster surface (Zeltner
// ManifoldVertex restricted to what the SS solver needs).
struct MneeVertex {
	Point p;
	Vector dpdu, dpdv;
	Normal n, gn, dndu, dndv;
	// Orthonormal tangents (Zeltner make_orthonormal result)
	Vector s, t;
	// Relative IOR (interior/exterior). 1.f == conductor (mirror).
	float eta;
};

inline void MneeCoordinateSystem(const Vector &n, Vector &s, Vector &t);

// Orthonormalize the surface parameterization (Zeltner
// ManifoldVertex::make_orthonormal)
inline void MneeOrthonormalize(MneeVertex &v) {
	const float invNorm1 = 1.f / sqrtf(v.dpdu.LengthSquared());
	v.dpdu *= invNorm1;
	v.dndu *= invNorm1;

	const float dp = Dot(v.dpdu, v.dpdv);
	const Vector dpdvTmp = v.dpdv - dp * v.dpdu;
	const Normal dndvTmp = v.dndv - dp * v.dndu;
	const float invNorm2 = 1.f / sqrtf(dpdvTmp.LengthSquared());
	v.dpdv = dpdvTmp * invNorm2;
	v.dndv = dndvTmp * invNorm2;

	v.s = v.dpdu;
	v.t = v.dpdv;
}

inline void MneeInitVertex(MneeVertex &v, const BSDF &bsdf, const float eta) {
	const HitPoint &hitPoint = bsdf.hitPoint;
	v.p = hitPoint.p;
	v.dpdu = hitPoint.dpdu;
	v.dpdv = hitPoint.dpdv;
	v.n = hitPoint.shadeN;
	v.gn = hitPoint.geometryN;
	v.dndu = hitPoint.dndu;
	v.dndv = hitPoint.dndv;
	v.eta = eta;

	// Meshes without UVs have degenerate differentials: fall back to an
	// orthonormal flat parameterization around the shading normal (curvature
	// terms are lost in this case, i.e. non-UV meshes are treated as flat).
	const float dpduLenSq = v.dpdu.LengthSquared();
	const float dpdvLenSq = v.dpdv.LengthSquared();
	if (dpduLenSq < 1e-24f || dpdvLenSq < 1e-24f ||
			Dot(v.dpdu, v.dpdv) * Dot(v.dpdu, v.dpdv) >
			0.9999f * 0.9999f * dpduLenSq * dpdvLenSq) {
		Vector s, t;
		MneeCoordinateSystem(Vector(v.n.x, v.n.y, v.n.z), s, t);
		v.dpdu = s;
		v.dpdv = t;
		v.dndu = Normal();
		v.dndv = Normal();
	}

	MneeOrthonormalize(v);
}

inline void MneeCoordinateSystem(const Vector &n, Vector &s, Vector &t) {
	if (fabsf(n.x) > fabsf(n.y)) {
		const float invNorm = 1.f / sqrtf(n.x * n.x + n.z * n.z);
		s = Vector(-n.z * invNorm, 0.f, n.x * invNorm);
	} else {
		const float invNorm = 1.f / sqrtf(n.y * n.y + n.z * n.z);
		s = Vector(0.f, n.z * invNorm, -n.y * invNorm);
	}
	t = Cross(n, s);
}

} // anonymous namespace

//------------------------------------------------------------------------------
// Half-vector constraint (Zeltner compute_step_halfvector, n_offset = 0)
//
// Residual C = (h.s, h.t) at the specular vertex; the Jacobian dC/dX is
// computed numerically (finite differences with a re-projection trace so the
// normal rotation / curvature is captured). The light-side Jacobian dC/dX2 is
// pure arithmetic on the fake orthonormal light frame (Zeltner point emitter
// handling).
//------------------------------------------------------------------------------

// etaOverride > 0.f forces the half-vector IOR ratio of the constraint
// (0.f = use the physical vertex eta). The geometric term is always
// evaluated with eta = +1 (the area-measure Jacobian of Zeltner's
// geometric_term, validated by the numpy prototype), while the Newton
// constraint may use the physical eta (e.g. -1 for an opposite-side mirror
// reflection).
static bool MneeResidual(const Point &x0p, const Point &lightPos,
		const MneeVertex &vtx, MneeVec2 &C, const float etaOverride = 0.f) {
	Vector wi = x0p - vtx.p;
	float r01 = wi.Length();
	if (r01 < 1e-3f)
		return false;
	wi /= r01;

	Vector wo = lightPos - vtx.p;
	float r12 = wo.Length();
	if (r12 < 1e-3f)
		return false;
	wo /= r12;

	float eta = (etaOverride > 0.f) ? etaOverride : vtx.eta;
	if (Dot(wi, vtx.gn) < 0.f)
		eta = 1.f / eta;
	Vector h = wi + eta * wo;
	if (eta != 1.f)
		h = -h;
	h *= 1.f / h.Length();

	C.x = Dot(vtx.s, h);
	C.y = Dot(vtx.t, h);
	return true;
}

static bool MneeReproject(
		luxrays::IntersectionDeviceRef device, SceneConstRef scene,
		const float time,
		const Point &pOff, const Normal &gn, const float gap,
		BSDF &reprojectBsdf) {
	// Start slightly above the surface (along the geometric normal) to avoid
	// self-intersection: the perturbed point is only ~O(eps^2) off the
	// (curved) surface, for a flat one it is exactly on it.
	const Point origin = pOff + gap * Vector(gn.x, gn.y, gn.z);
	Ray reprojectRay(origin, Vector(-gn.x, -gn.y, -gn.z), 0.f, .1f, time);
	RayHit reprojectHit;
	Spectrum reprojectThru;
	PathVolumeInfo reprojectVol;
	if (!scene.Intersect(IntersectionDevicePtr(&device),
			INDIRECT_RAY, &reprojectVol, .5f, &reprojectRay,
			&reprojectHit, &reprojectBsdf, &reprojectThru, nullptr,
			nullptr, false))
		return false;
	return true;
}

static bool MneeConstraintWithJacobian(
		luxrays::IntersectionDeviceRef device, SceneConstRef scene,
		const float time,
		const Point &x0p, const Point &lightPos,
		const MneeVertex &vtx,
		MneeVec2 &C, MneeMat2 &jac, const float etaOverride = 0.f) {
	if (!MneeResidual(x0p, lightPos, vtx, C, etaOverride))
		return false;

	const float eps = Max(1e-5f, 1e-4f * Distance(x0p, vtx.p));
	const float C0[2] = { C.x, C.y };
	for (int k = 0; k < 2; ++k) {
		// Perturb along the (orthonormal) tangent, re-project onto the same
		// surface so the normal rotation (curvature) is included
		const Point pPert = vtx.p + ((k == 0) ? eps * vtx.dpdu : eps * vtx.dpdv);
		BSDF pertBsdf;
		if (!MneeReproject(device, scene, time, pPert, vtx.gn, .5f * eps, pertBsdf))
			return false;

		MneeVertex vPert;
		MneeInitVertex(vPert, pertBsdf, vtx.eta);

		MneeVec2 CP;
		if (!MneeResidual(x0p, lightPos, vPert, CP, etaOverride))
			return false;

		if (k == 0) {
			jac.a11 = (CP.x - C0[0]) / eps;
			jac.a21 = (CP.y - C0[1]) / eps;
		} else {
			jac.a12 = (CP.x - C0[0]) / eps;
			jac.a22 = (CP.y - C0[1]) / eps;
		}
	}

	return true;
}

static MneeMat2 MneeLightJacobian(const Point &x0p, const Point &lightPos,
		const MneeVertex &vtx, const float eps, const float etaOverride = 0.f) {
	// Fake light vertex frame: orthonormal s/t built on the light->x1
	// direction (Zeltner emitter_interaction_to_vertex, point emitter branch)
	Vector dLight = vtx.p - lightPos;
	const float r = dLight.Length();
	if (r < 1e-3f)
		return MneeMat2{ 0.f, 0.f, 0.f, 0.f };

	Vector s2, t2;
	MneeCoordinateSystem(dLight, s2, t2);

	const Vector wi = Normalize(x0p - vtx.p);
	float eta = (etaOverride > 0.f) ? etaOverride : vtx.eta;
	if (Dot(wi, vtx.gn) < 0.f)
		eta = 1.f / eta;

	Vector h = wi + eta * Normalize(lightPos - vtx.p);
	if (eta != 1.f)
		h = -h;
	h *= 1.f / h.Length();

	const float C0[2] = { Dot(vtx.s, h), Dot(vtx.t, h) };
	float CP[2][2];
	for (int k = 0; k < 2; ++k) {
		const Point lightPosP = lightPos + ((k == 0) ? eps * s2 : eps * t2);
		const Vector woP = Normalize(lightPosP - vtx.p);
		Vector hP = wi + eta * woP;
		if (eta != 1.f)
			hP = -hP;
		hP *= 1.f / hP.Length();
		CP[k][0] = Dot(vtx.s, hP);
		CP[k][1] = Dot(vtx.t, hP);
	}

	return MneeMat2{
		(CP[0][0] - C0[0]) / eps, (CP[1][0] - C0[0]) / eps,
		(CP[0][1] - C0[1]) / eps, (CP[1][1] - C0[1]) / eps
	};
}

// Analytic constraint Jacobians and geometric term (Zeltner
// geometric_term structure, evaluated with the physical vertex eta; flat and
// curved surfaces via the dndu/dndv curvature terms of the mesh differentials).
// Returns dw0/dx1 * |det(inv(J1) * J2)|. If jacVertex (the constraint
// Jacobian dC/dX used by the Newton step) is not null it receives J1.
static float MneeGeometricTermWithJacobians(const Point &x0p,
		const Point &lightPos, const MneeVertex &vtx, MneeMat2 *jacVertex,
		float *det1Out = nullptr, float *det2Out = nullptr) {
	Vector wi = x0p - vtx.p;
	const float r01 = wi.Length();
	if (r01 < 1e-3f)
		return 0.f;
	wi /= r01;

	Vector wo = lightPos - vtx.p;
	const float r12 = wo.Length();
	if (r12 < 1e-3f)
		return 0.f;
	wo /= r12;

	float eta = vtx.eta;
	if (Dot(wi, vtx.gn) < 0.f)
		eta = 1.f / eta;
	Vector h = wi + eta * wo;
	if (eta != 1.f)
		h = -h;
	const float ilh = 1.f / h.Length();
	h *= ilh;
	const float ilo = (1.f / r12) * eta * ilh;
	const float ili = (1.f / r01) * ilh;

	const Vector &s = vtx.s;
	const Vector &t = vtx.t;

	// dC/dx1 (vertex side, with curvature terms)
	Vector dhDdu = -vtx.dpdu * (ili + ilo) +
			wi * (Dot(wi, vtx.dpdu) * ili) +
			wo * (Dot(wo, vtx.dpdu) * ilo);
	Vector dhDdv = -vtx.dpdv * (ili + ilo) +
			wi * (Dot(wi, vtx.dpdv) * ili) +
			wo * (Dot(wo, vtx.dpdv) * ilo);
	dhDdu -= h * Dot(dhDdu, h);
	dhDdv -= h * Dot(dhDdv, h);
	if (eta != 1.f) {
		dhDdu = -dhDdu;
		dhDdv = -dhDdv;
	}
	const float dotHN = Dot(h, vtx.n);
	const float dotHDndu = Dot(h, vtx.dndu);
	const float dotHDndv = Dot(h, vtx.dndv);
	const float dotDpduN = Dot(vtx.dpdu, vtx.n);
	const float dotDpdvN = Dot(vtx.dpdv, vtx.n);
	const MneeMat2 j1{
		Dot(dhDdu, s) - Dot(vtx.dpdu, vtx.dndu) * dotHN - dotDpduN * dotHDndu,
		Dot(dhDdv, s) - Dot(vtx.dpdu, vtx.dndv) * dotHN - dotDpduN * dotHDndv,
		Dot(dhDdu, t) - Dot(vtx.dpdv, vtx.dndu) * dotHN - dotDpdvN * dotHDndu,
		Dot(dhDdv, t) - Dot(vtx.dpdv, vtx.dndv) * dotHN - dotDpdvN * dotHDndv
	};

	// dC/dx2 (fake light frame, from the light toward the vertex)
	Vector dLight = vtx.p - lightPos;
	const float rl = dLight.Length();
	if (rl < 1e-3f)
		return 0.f;
	dLight /= rl;
	Vector s2, t2;
	MneeCoordinateSystem(dLight, s2, t2);
	Vector dhDdu2 = ilo * (s2 - wo * Dot(wo, s2));
	Vector dhDdv2 = ilo * (t2 - wo * Dot(wo, t2));
	dhDdu2 -= h * Dot(dhDdu2, h);
	dhDdv2 -= h * Dot(dhDdv2, h);
	if (eta != 1.f) {
		dhDdu2 = -dhDdu2;
		dhDdv2 = -dhDdv2;
	}
	const MneeMat2 j2{
		Dot(dhDdu2, s), Dot(dhDdv2, s),
		Dot(dhDdu2, t), Dot(dhDdv2, t)
	};

	const float det1 = MneeDet(j1);
	const float det2 = MneeDet(j2);

	if (jacVertex)
		*jacVertex = j1;
	if (det1Out)
		*det1Out = det1;
	if (det2Out)
		*det2Out = det2;
	if (fabs(det1) < 1e-9f || fabs(det2) < 1e-15f)
		return 0.f;

	const float dx1Dx2 = fabsf(det2 / det1);
	const Vector d01 = x0p - vtx.p;
	const float r01sq = d01.LengthSquared();
	const float dw0Dx1 = fabsf(Dot(d01, vtx.gn)) / (sqrtf(r01sq) * r01sq);
	return dw0Dx1 * dx1Dx2;
}

static float MneeGeometricTerm(const Point &x0p, const Point &lightPos,
		const MneeVertex &vtx, float *det1Out = nullptr,
		float *det2Out = nullptr) {
	return MneeGeometricTermWithJacobians(x0p, lightPos, vtx, nullptr,
			det1Out, det2Out);
}

// Term-by-term diagnostic of the assembly. Enabled with LUX_MNEE_DEBUG=1: for
// every accepted contribution it prints the receiver, the solved specular
// vertex and the light position together with each factor of the estimate, so
// an external analytic model (dev-tools/mnee_glass_term_model.py) can be
// diffed against the live code term by term.
static bool MneeDebugEnabled() {
	static const bool enabled = (getenv("LUX_MNEE_DEBUG") != nullptr);
	return enabled;
}
// Instrument (revert): per-REJECTED-attempt dump (x0, light) for the
// Newton-failure energy measurement (does the -5% glass shortfall sit in
// rejected attempts?). Pair with LUX_MNEE_DEBUG accepted lines.
static bool MneeRejXEnabled() {
	static const bool enabled = (getenv("LUX_MNEE_REJX") != nullptr);
	return enabled;
}
#define MNEE_REJX(why) do { \
		if (MneeRejXEnabled()) { \
			printf("MNEE_REJX %s x0=%.9g %.9g %.9g y=%.9g %.9g %.9g\n", why, \
				x0p.x, x0p.y, x0p.z, lightPos.x, lightPos.y, lightPos.z); \
			fflush(stdout); \
		} \
	} while (0)

//------------------------------------------------------------------------------
// PathTracer::MNEEDirectSampling
//------------------------------------------------------------------------------


bool PathTracer::MNEEDirectSampling(
		luxrays::IntersectionDeviceRef device,
		SceneConstRef scene,
		const float time,
		const EyePathInfo &pathInfo,
		const luxrays::Spectrum &pathThroughput,
		const BSDF &bsdf,
		LightSourceConstRef light, const float lightPickPdf, const float risScale,
		const luxrays::Ray &shadowRay, const float directPdfW0,
		const luxrays::RayHit &shadowRayHit,
		const BSDF &shadowBsdf, PathVolumeInfo &volInfo,
		const float u1, const float u2, const float u3, const float u4,
		SampleResult *sampleResult) const {
	// The occluder material defines the specular event
	const MaterialType seedMatType = shadowBsdf.GetMaterialType();
	const Point &x0p = bsdf.hitPoint.p;

	// Light position: the shadow ray maxt has been rewritten by Scene::Intersect
	// to the occluder distance, so recover the light distance from directPdfW
	// (= squared distance to the light, preserved by Illuminate).
	const Point lightPos = shadowRay.o + shadowRay.d * sqrtf(directPdfW0);
	observer_ptr<const MirrorMaterial> mirrorMat = nullptr;
	observer_ptr<const GlassMaterial> glassMat = nullptr;
	float etaVertex = 1.f;
	if (seedMatType == MIRROR) {
		mirrorMat = dynamic_observer_cast<const MirrorMaterial>(shadowBsdf.GetMaterial());
		if (!mirrorMat)
			{ return false; }

		// Generalized half-vector IOR ratio for a conductor: +1 when the two
		// endpoints are on the same side of the surface (h = wi + wo, the
		// same-side reflection case), -1 when they are on opposite sides
		// (h = wi - wo).
		//
		// Only the same-side relation is a reflection: the reflected ray
		// leaves the surface with the normal component of its direction
		// flipped, so it reaches the receiver on the same side of the tangent
		// plane the light is on. Endpoints on opposite sides would need the
		// surface to transmit, which a mirror does not do; solving that case
		// injects light where light tracing, BIDIR and deep path tracing all
		// measure exactly zero (dev-tools/sota_p1_mnee_mirror_physics_test.py,
		// "plane" case).
		const Normal &gn1s = shadowBsdf.hitPoint.geometryN;
		const Vector toX0 = x0p - shadowBsdf.hitPoint.p;
		const Vector toY = lightPos - shadowBsdf.hitPoint.p;
		etaVertex = (Dot(toX0, gn1s) * Dot(toY, gn1s) > 0.f) ? 1.f : -1.f;
		if (etaVertex != 1.f)
			return false;
	} else if (seedMatType == GLASS) {
		glassMat = dynamic_observer_cast<const GlassMaterial>(shadowBsdf.GetMaterial());
		if (!glassMat)
			{ return false; }
		// Dispersive glass: the manifold uses a single IOR ratio, skip
		if (glassMat->GetCauchyB() &&
				glassMat->GetCauchyB()->GetFloatValue(shadowBsdf.hitPoint) > 0.f)
			return false;

		const float nc = ExtractExteriorIors(shadowBsdf.hitPoint, glassMat->GetExteriorIOR());
		const float nt = ExtractInteriorIors(shadowBsdf.hitPoint, glassMat->GetInteriorIOR());
		if (nt <= 0.f || nc <= 0.f)
			return false;
		etaVertex = nt / nc;
	} else
		{ return false; }

	// Count as an attempt (past the material gate)
	static int attemptCheck = 0;
	++attemptCheck;

	//------------------------------------------------------------------------------
	// Newton solve (Zeltner newton_solver, n_offset = 0, step_scale = 1)
	//------------------------------------------------------------------------------

	MneeVertex vtx;
	MneeInitVertex(vtx, shadowBsdf, etaVertex);

	// The shadow-ray seed lies exactly on the x0->y line, where the
	// generalized half-vector degenerates (h = wi + eta*wo = (1 - eta)*wi).
	// For a mirror (eta == 1) the line seed carries no information at all and
	// the Newton iteration diverges from there. Seed the reflection chain
	// with the first hit of the ray from x0 toward the light mirrored across
	// the tangent plane at the shadow hit: for a flat mirror this is the
	// exact solution, for curved reflectors the best local approximation
	// (Zeltner's "Modified MNEE" idea, with the mirrored light instead of the
	// shape bbox center). For eta != 1 the line seed is informative
	// ((1 - eta) * wi != 0) and is kept.
	if (etaVertex == 1.f) {
		const Normal &gn1 = shadowBsdf.hitPoint.geometryN;
		const Point x1Line = shadowBsdf.hitPoint.p;
		const Vector x1ToLight = lightPos - x1Line;
		const float proj = 2.f * Dot(x1ToLight, gn1);
		const Vector mirroredLight = Vector(lightPos.x, lightPos.y, lightPos.z) -
				proj * Vector(gn1.x, gn1.y, gn1.z);
		const Vector dSeed = Normalize(mirroredLight -
				Vector(x0p.x, x0p.y, x0p.z));
		Ray seedRay(bsdf.GetRayOrigin(dSeed), dSeed, 0.f,
				numeric_limits<float>::infinity(), time);
		RayHit seedHit;
		BSDF seedBsdf;
		Spectrum seedThru;
		PathVolumeInfo seedVol = volInfo;
		if (scene.Intersect(IntersectionDevicePtr(&device),
				INDIRECT_RAY, &seedVol, .5f, &seedRay,
				&seedHit, &seedBsdf, &seedThru, nullptr, nullptr, false) &&
				seedHit.meshIndex == shadowRayHit.meshIndex) {
			MneeInitVertex(vtx, seedBsdf, etaVertex);
		}
	}

	// Small deterministic tangent offset so the initial half-vector never
	// degenerates exactly (e.g. an off-axis light still on the x0->center
	// line). The offset only selects the Newton basin; the solve itself is
	// pulled to the exact constraint solution.
	const float seedShift = Max(1e-4f, 1e-3f * Distance(x0p, vtx.p));
	vtx.p += (vtx.dpdu + vtx.dpdv) * (seedShift / sqrtf(2.f));

	bool solved = false;
	float beta = 1.f;
	BSDF finalBsdf = shadowBsdf;
	MneeVec2 residual{ 0.f, 0.f };
	MneeMat2 jac{ 0.f, 0.f, 0.f, 0.f };

	u_int iteration = 0;
	while (iteration < mneeMaxIterations) {
		// Constraint residual and analytic Jacobian of the current vertex
		if (!MneeResidual(x0p, lightPos, vtx, residual)) {
			break;
		}
		const float g = MneeGeometricTermWithJacobians(x0p, lightPos, vtx, &jac);
	
		if (sqrtf(residual.x * residual.x + residual.y * residual.y) < 3e-4f) {
			solved = true;
			break;
		}

		const float det = MneeDet(jac);
		if (fabs(det) < 1e-9f) {
			break;
		}
		const MneeMat2 invJac = MneeInverse(jac, det);
		MneeVec2 dX{ invJac.a11 * residual.x + invJac.a12 * residual.y,
				invJac.a21 * residual.x + invJac.a22 * residual.y };
		// Clamp the step near singular configurations (det(J1) -> 0 along the
		// mirror axis): the Newton direction is still the descent direction
		// but its magnitude explodes, so cap it to a fraction of the vertex
		// distance and let the line search find the residual decrease.
		const float dXNorm = sqrtf(dX.x * dX.x + dX.y * dX.y);
		const float dXMax = .25f * Distance(x0p, vtx.p);
		if (dXNorm > dXMax) {
			dX.x *= dXMax / dXNorm;
			dX.y *= dXMax / dXNorm;
		}

		// Line search: accept the step only when the residual decreases and
		// the re-projection stays on the same mesh (Zeltner's beta
		// backtracking, extended with a residual decrease check; the
		// half-vector residual is strongly nonlinear for coarse seeds, so
		// this keeps the iteration inside the Newton basin).
		const float resNorm = sqrtf(residual.x * residual.x + residual.y * residual.y);
		bool stepAccepted = false;
		while (beta > 1e-2f) {
			const Point pProp = vtx.p - beta * (vtx.dpdu * dX.x + vtx.dpdv * dX.y);
			const Vector dProp = Normalize(pProp - x0p);

			Ray propRay(bsdf.GetRayOrigin(dProp), dProp, 0.f,
					numeric_limits<float>::infinity(), time);
			RayHit propHit;
			BSDF propBsdf;
			Spectrum propThru;
			PathVolumeInfo propVol = volInfo;
			if (!scene.Intersect(IntersectionDevicePtr(&device),
					INDIRECT_RAY, &propVol, .5f, &propRay,
					&propHit, &propBsdf, &propThru, nullptr, nullptr, false))
				break;

			if (propHit.meshIndex != shadowRayHit.meshIndex) {
				beta *= .5f;
				++iteration;
				continue;
			}

			MneeVertex vProp;
			MneeInitVertex(vProp, propBsdf, etaVertex);
			MneeVec2 resProp;
			if (!MneeResidual(x0p, lightPos, vProp, resProp)) {
				beta *= .5f;
				++iteration;
				continue;
			}
			const float resPropNorm = sqrtf(resProp.x * resProp.x + resProp.y * resProp.y);

			if (resPropNorm < resNorm) {
				beta = Min(1.f, 2.f * beta);
				vtx = vProp;
				finalBsdf = propBsdf;
				stepAccepted = true;
				break;
			}

			beta *= .5f;
			++iteration;
		}
		if (!stepAccepted)
			break;
		++iteration;
	}

	if (!solved) {
		MNEE_REJX("newton");
		return false;
	}

	//------------------------------------------------------------------------------
	// Post-solve validity check (Zeltner newton_solver tail): the half-vector
	// formulation can converge to a solution of the wrong specular mode
	//------------------------------------------------------------------------------

	const Vector wi = Normalize(x0p - vtx.p);
	const Vector wo = Normalize(lightPos - vtx.p);
	const float cosX = Dot(vtx.gn, wi);
	const float cosY = Dot(vtx.gn, wo);
	const bool refraction = (cosX * cosY < 0.f);
	if (mirrorMat) {
		// Mirror: only a same-side reflection is physical (see the etaVertex
		// gate above). Reject solutions with the opposite side relation - the
		// half-vector formulation can converge to such a mode, and accepting
		// it produced light where none exists.
		if (refraction) {
			MNEE_REJX("mode-mirror");
			return false;
		}
	} else if (!refraction) {
		// Glass: only refraction solutions are supported (Zeltner SS handles
		// dielectric transmission; external dielectric reflection is out of
		// scope)
		MNEE_REJX("mode-glass");
		return false;
	}

	//------------------------------------------------------------------------------
	// Specular factor at the solved vertex, with LuxCore's own material code
	// (mirror: Kr; glass: (1 - F) * eta^2 through the renderer's static
	// evaluation, in the eye-path convention: hitPoint.fromLight == false and
	// localFixedDir = wi)
	//------------------------------------------------------------------------------

	Spectrum specFactor;
	BSDFEvent specEvent;
	if (mirrorMat) {
		specFactor = mirrorMat->GetKr()->GetSpectrumValue(finalBsdf.hitPoint).Clamp(0.f, 1.f);
		specEvent = SPECULAR | REFLECT;
	} else {
		const Spectrum kt = glassMat->GetKt()->GetSpectrumValue(finalBsdf.hitPoint).Clamp(0.f, 1.f);
		const float nc = ExtractExteriorIors(finalBsdf.hitPoint, glassMat->GetExteriorIOR());
		const float nt = ExtractInteriorIors(finalBsdf.hitPoint, glassMat->GetInteriorIOR());

		const Vector localFixedDir = finalBsdf.GetFrame().ToLocal(wi);
		Vector localSampledDir;
		specFactor = GlassMaterial::EvalSpecularTransmission(finalBsdf.hitPoint,
				localFixedDir, 0.f, kt, nc, nt, 0.f, &localSampledDir);
		specEvent = SPECULAR | TRANSMIT;
	}
	if (specFactor.Black()) {
		MNEE_REJX("spec");
		return false;
	}

	//------------------------------------------------------------------------------
	// Geometric term (analytic): dw0/dx1 * |det(inv(dc1/dx1) * dc1/dx2)| with
	// Zeltner's curvature structure evaluated with the physical constraint
	// eta. Validated numerically:
	//  - same-side reflection (eta = +1): E[G]/virtual-light truth = 1.00000
	//    (numpy prototype, dev-tools/mnee_design.md);
	//  - opposite-side mirror reflection (eta = -1): G = 1/r'^2 with the
	//    virtual-light truth (derive/curved checks in the session notes).
	// No clamping: the true Jacobian exceeds 1 near glancing configurations
	// and clamping would bias the estimate.
	//------------------------------------------------------------------------------
	float mneeDet1 = 0.f, mneeDet2 = 0.f;
	const float geometricTerm = MneeGeometricTerm(x0p, lightPos, vtx,
			&mneeDet1, &mneeDet2);
	if (geometricTerm <= 0.f || isnan(geometricTerm) || isinf(geometricTerm)) {
		MNEE_REJX("geoterm");
		return false;
	}

	//------------------------------------------------------------------------------
	// Second segment x1 -> y: Illuminate at the specular vertex and check
	// visibility (the chain is valid only if y is directly visible from x1)
	//------------------------------------------------------------------------------

	// Volume state after the specular event at x1
	PathVolumeInfo volSeg2 = volInfo;
	volSeg2.Update(specEvent, finalBsdf);

	Ray shadowRay2;
	float directPdfW2;
	const Spectrum lightRadiance2 = light.Illuminate(scene, finalBsdf, time,
			u1, u2, u3, shadowRay2, directPdfW2);
	if (lightRadiance2.Black()) {
		MNEE_REJX("seg2black");
		return false;
	}
	verify(!isnan(directPdfW2) && !isinf(directPdfW2));

	RayHit occlHit;
	BSDF occlBsdf;
	Spectrum seg2Throughput;
	PathVolumeInfo volSeg2Trace = volSeg2;
	if (scene.Intersect(IntersectionDevicePtr(&device),
			SHADOW_RAY, &volSeg2Trace, u4, &shadowRay2,
			&occlHit, &occlBsdf, &seg2Throughput, nullptr, nullptr, true)) {
		return false;
	}

	//------------------------------------------------------------------------------
	// Contribution assembly
	//------------------------------------------------------------------------------

	// Receiver BSDF toward x1 (includes the cosine, like the plain estimator)
	BSDFEvent receiverEvent;
	float bsdfPdfW;
	const Spectrum bsdfEval0 = bsdf.Evaluate(Normalize(vtx.p - x0p),
			&receiverEvent, &bsdfPdfW);
	if (bsdfEval0.Black()) {
		return false;
	}
	verify(!isnan(bsdfPdfW) && !isinf(bsdfPdfW));

	// Path depth: one vertex for the receiver event, one for the specular
	// event at x1 (bookkeeping only; the visibility flags are used by the
	// MIS weight in the plain estimator, MNEE has no MIS partner)
	PathDepthInfo mneeDepthInfo = pathInfo.depth;
	mneeDepthInfo.IncDepths(receiverEvent);
	mneeDepthInfo.IncDepths(specEvent);

	// MNEE light weight. The second-segment measure conversion depends on the
	// constraint that was solved:
	//
	//  - plain half-vector (vtx eta == 1, i.e. the same-side reflection seed
	//    h = wi + wo): the analytic geometric term converts to the light's
	//    area measure and the r12^2 factor rides here. Validated against an
	//    independent virtual-light scene at 0.013% (dev-tools/
	//    sota_p1_mnee_test.py, flat mirror replacing the mirror with the
	//    mirrored point light).
	//  - generalized half-vector (eta != 1: dielectric transmission with
	//    eta = nt/nc, and the opposite-side mirror law with eta = -1): the
	//    analytic geometric term already carries the full conversion and
	//    directPdfW2 must NOT be multiplied. Validated per pixel against
	//    converged light tracing plus the analytic per-point radiance
	//    (dev-tools/mnee_glass_pixel_compare.py: including r12^2 gives a 6.2x
	//    shape error over one image, excluding it agrees to 1.3x, the same
	//    spread with which light tracing itself matches the analytic
	//    radiance).
	//
	// risScale keeps ReSTIR RIS unbiased.
	const bool plainHalfVector = (etaVertex == 1.f);
	const Spectrum lightWeight = lightRadiance2 *
			((plainHalfVector ? directPdfW2 : 1.f) * risScale / lightPickPdf);
	const Spectrum incomingRadiance = bsdfEval0 * specFactor * geometricTerm *
			lightWeight * seg2Throughput;
	verify(!incomingRadiance.IsNaN() && !incomingRadiance.IsInf());

	if (MneeDebugEnabled()) {
		const Vector dbgWi = Normalize(x0p - vtx.p);
		const Vector dbgWo = Normalize(lightPos - vtx.p);
		const float dbgR01 = Distance(x0p, vtx.p);
		const float dbgR12 = Distance(lightPos, vtx.p);
		// Analytic per-point radiance leaving x0 toward the camera for this
		// path (the quantity a light-transport estimator measures):
		//   (kd/pi)*|cos| * (1-F)*eta_t^2 * I / r12^2
		// built from the live factors, so the two can be compared per pixel
		// without knowing the camera projection.
		const float dbgTruth = bsdfEval0.Filter() * specFactor.Filter() *
				lightRadiance2.Filter() / (dbgR12 * dbgR12);
		printf("MNEE_DBG film=%.4g %.4g x0=%.9g %.9g %.9g | x1=%.9g %.9g %.9g | "
				"y=%.9g %.9g %.9g | r01=%.9g r12=%.9g eta=%.9g dotWiGn=%.9g "
				"dotWiShadeN=%.9g cosWiGn=%.9g cosWoGn=%.9g | det1=%.9g "
				"det2=%.9g G=%.9g dw0dx1=%.9g | spec=%.9g bsdf0=%.9g lr=%.9g "
				"dpdf=%.9g pick=%.9g seg2=%.9g in=%.9g truth=%.9g\n",
				sampleResult->filmX, sampleResult->filmY,
				x0p.x, x0p.y, x0p.z, vtx.p.x, vtx.p.y, vtx.p.z,
				lightPos.x, lightPos.y, lightPos.z,
				dbgR01, dbgR12, etaVertex,
				Dot(dbgWi, vtx.gn), Dot(dbgWi, vtx.n),
				fabsf(Dot(dbgWi, vtx.gn)), fabsf(Dot(dbgWo, vtx.gn)),
				mneeDet1, mneeDet2, geometricTerm,
				fabsf(Dot(x0p - vtx.p, vtx.gn)) / (dbgR01 * dbgR01),
				specFactor.Filter(), bsdfEval0.Filter(),
				lightRadiance2.Filter(), directPdfW2, lightPickPdf,
				seg2Throughput.Filter(), incomingRadiance.Filter(), dbgTruth);
		fflush(stdout);
	}

	sampleResult->AddDirectLight(light.GetID(), specEvent, pathThroughput,
			incomingRadiance, 1.f);

	return true;
}

//------------------------------------------------------------------------------
// MNEE, multi-specular chain (N >= 2 delta specular vertices)
//
// The single vertex solver above covers eye -> x0 -> x1 -> y. Closed glass
// slabs and glass balls need more: the light reaches the receiver through two
// or more refractions, so the shadow ray from x0 is stopped by a specular
// surface that cannot see the light either (the next face is in the way).
//
// This solver generalizes the chain to N vertices
//
//   eye -> x0 -> x1 -> ... -> xN -> y
//
// with one generalized half-vector constraint per vertex. The numerics are
// ported from the numpy prototype that was verified before the port
// (dev-tools/mnee_ms_proto.py, notes in dev-tools/mnee_design.md section 4b):
//
//   - the constraint at vertex i involves only x_{i-1}, x_i and x_{i+1}, so
//     the finite difference Jacobian is exactly block tridiagonal (verified
//     there: no leakage into the far blocks);
//   - the Newton step solves that system with the block Thomas recursion
//     (verified against the dense solution), and the matrix right hand side
//     variant of the same recursion against the dense determinant (2.6e-14);
//   - the chain geometric term is |det(dx_1/dy)|, i.e. the (1, N) block of the
//     inverse constraint Jacobian times the light's own Jacobian; for N = 1 it
//     reduces to the single vertex term |det(J2) / det(J1)|;
//   - Newton converges with Snell's law satisfied at every vertex (1.3e-09).
//
// The topology is discovered by tracing the straight ray from x0 toward the
// light and collecting the delta specular surfaces it pierces (the reference's
// mnee_init seeding); the solve then bends the vertices onto the refraction
// manifold. The delta mode of a vertex is fixed by its material (a mirror
// reflects, glass transmits) and the constraint IOR ratio follows the side of
// the previous point, exactly as in the single vertex solver.
//
// Disjointness with the plain estimator follows the single vertex argument:
// the shadow ray from x0 is blocked, and the chain cannot be sampled by
// forward BSDF sampling of a positional delta light. Light tracing and BIDIR
// do measure these paths, so they are the unbiased reference (see the slab
// case of dev-tools/sota_p1_mnee_glass_test.py).
//------------------------------------------------------------------------------

// One vertex of the MNEE chain
struct MneeChainVertex {
	MneeVertex v;
	// The material's relative IOR (interior/exterior), 1 for a conductor. The
	// constraint flips it when the previous point lies on the -geometryN side,
	// exactly like MneeResidual does for a single vertex.
	float etaVertex;
	Spectrum specFactor;
	BSDFEvent specEvent;
	BSDF bsdf;
	observer_ptr<const MirrorMaterial> mirrorMat;
	observer_ptr<const GlassMaterial> glassMat;
};

static const u_int MNEE_MS_MAX_VERTICES = 4;

// Block tridiagonal constraint Jacobian: the 2x2 blocks of x_{i-1}, x_i, x_{i+1}
struct MneeJacobianBlock {
	MneeMat2 prev, cur, next;
};

// 2x2 helpers on top of the anonymous namespace ones
inline MneeMat2 MneeMatSub(const MneeMat2 &A, const MneeMat2 &B) {
	return MneeMat2{ A.a11 - B.a11, A.a12 - B.a12, A.a21 - B.a21, A.a22 - B.a22 };
}

inline MneeVec2 MneeMatVec(const MneeMat2 &A, const MneeVec2 &v) {
	return MneeVec2{ A.a11 * v.x + A.a12 * v.y, A.a21 * v.x + A.a22 * v.y };
}

// The generalized half-vector constraint at one chain vertex: h = wi + eta * wo
// must be parallel to the surface normal, with eta the material's relative IOR
// flipped when the previous point lies on the -geometryN side (LuxCore's
// exterior IOR sits on the +geometryN side; the same convention as
// MneeResidual). wi points at the previous chain point, wo at the next one.
// Returns false for degenerate configurations.
static bool MneeChainResidual(const Point &pPrev, const Point &pNext,
		const MneeVertex &v, const float etaVertex, MneeVec2 &C) {
	Vector wi = pPrev - v.p;
	const float r0 = wi.Length();
	if (r0 < 1e-4f)
		return false;
	wi *= 1.f / r0;

	Vector wo = pNext - v.p;
	const float r1 = wo.Length();
	if (r1 < 1e-4f)
		return false;
	wo *= 1.f / r1;

	float eta = etaVertex;
	if (Dot(wi, v.gn) < 0.f)
		eta = 1.f / eta;

	Vector h = wi + eta * wo;
	if (eta != 1.f)
		h = -h;
	const float l = h.Length();
	// Collinear segments make the half-vector degenerate (for eta == 1,
	// wi == -wo gives h == 0): the constraint carries no information there.
	if (l < 1e-5f)
		return false;
	h *= 1.f / l;

	C.x = Dot(v.s, h);
	C.y = Dot(v.t, h);
	return true;
}

// Residuals of the whole chain (no scene access)
static bool MneeChainResiduals(const Point &x0p, const Point &lightPos,
		const MneeChainVertex *chain, const u_int n,
		MneeVec2 *residual, float &maxResidual) {
	Point pts[MNEE_MS_MAX_VERTICES + 2];
	pts[0] = x0p;
	for (u_int i = 0; i < n; ++i)
		pts[i + 1] = chain[i].v.p;
	pts[n + 1] = lightPos;

	maxResidual = 0.f;
	for (u_int i = 0; i < n; ++i) {
		if (!MneeChainResidual(pts[i], pts[i + 2], chain[i].v, chain[i].etaVertex,
				residual[i]))
			return false;
		maxResidual = Max(maxResidual, sqrtf(residual[i].x * residual[i].x +
				residual[i].y * residual[i].y));
	}

	return true;
}

// Finite difference constraint Jacobian of the whole chain. Perturbing vertex i
// changes the constraints at i-1, i and i+1 only, so the result is exactly
// block tridiagonal. The perturbed vertex is re-projected onto its surface, so
// the normal rotation (curvature) is included, exactly like the single vertex
// solver's MneeConstraintWithJacobian.
static bool MneeChainJacobian(
		luxrays::IntersectionDeviceRef device, SceneConstRef scene,
		const float time, const Point &x0p, const Point &lightPos,
		const MneeChainVertex *chain, const u_int n,
		MneeJacobianBlock *blocks, MneeVec2 *residual, float &maxResidual) {
	if (!MneeChainResiduals(x0p, lightPos, chain, n, residual, maxResidual))
		return false;

	const MneeMat2 zero{ 0.f, 0.f, 0.f, 0.f };
	for (u_int i = 0; i < n; ++i) {
		blocks[i].prev = zero;
		blocks[i].cur = zero;
		blocks[i].next = zero;
	}

	// Note: the finite difference Jacobian is computed even when the residuals
	// are already zero (i.e. the straight line seed is the exact solution, as
	// it is for a flat interface at normal incidence). The geometric term needs
	// the Jacobian regardless, so returning early here would drop exactly those
	// contributions.
	Point pts[MNEE_MS_MAX_VERTICES + 2];
	pts[0] = x0p;
	for (u_int i = 0; i < n; ++i)
		pts[i + 1] = chain[i].v.p;
	pts[n + 1] = lightPos;

	for (u_int i = 0; i < n; ++i) {
		const float eps = Max(1e-5f, 1e-4f * Distance(pts[i], pts[i + 1]));
		for (u_int k = 0; k < 2; ++k) {
			const Point pPert = chain[i].v.p +
					((k == 0) ? eps * chain[i].v.dpdu : eps * chain[i].v.dpdv);
			BSDF pertBsdf;
			if (!MneeReproject(device, scene, time, pPert, chain[i].v.gn, .5f * eps, pertBsdf))
				return false;

			MneeVertex pertV;
			MneeInitVertex(pertV, pertBsdf, chain[i].etaVertex);

			Point ptsP[MNEE_MS_MAX_VERTICES + 2];
			for (u_int t = 0; t < n + 2; ++t)
				ptsP[t] = pts[t];
			ptsP[i + 1] = pertV.p;

			for (int j = (int)i - 1; j <= (int)i + 1; ++j) {
				if ((j < 0) || (j >= (int)n))
					continue;

				MneeVec2 CP;
				const MneeVertex &vj = (j == (int)i) ? pertV : chain[j].v;
				if (!MneeChainResidual(ptsP[j], ptsP[j + 2], vj, chain[j].etaVertex, CP))
					return false;

				MneeMat2 *block;
				if (j == (int)i - 1)
					block = &blocks[j].next;
				else if (j == (int)i)
					block = &blocks[i].cur;
				else
					block = &blocks[j].prev;

				const float dCx = (CP.x - residual[j].x) / eps;
				const float dCy = (CP.y - residual[j].y) / eps;
				if (k == 0) {
					block->a11 = dCx;
					block->a21 = dCy;
				} else {
					block->a12 = dCx;
					block->a22 = dCy;
				}
			}
		}
	}

	return true;
}

// Block Thomas decomposition (Kaplanyan 2014 supplement fig. 2, the reference's
// invert_tridiagonal_step). Returns false on a singular diagonal block.
static bool MneeTridiagonalInvert(const MneeJacobianBlock *blocks, const u_int n,
		MneeMat2 *tmp, MneeMat2 *invLambda) {
	if (n == 0)
		return false;

	float det = MneeDet(blocks[0].cur);
	if (fabsf(det) < 1e-12f)
		return false;
	invLambda[0] = MneeInverse(blocks[0].cur, det);
	tmp[0] = blocks[0].prev;

	for (u_int i = 1; i < n; ++i) {
		tmp[i] = MneeMul(blocks[i].prev, invLambda[i - 1]);
		const MneeMat2 m = MneeMatSub(blocks[i].cur,
				MneeMul(tmp[i], blocks[i - 1].next));
		det = MneeDet(m);
		if (fabsf(det) < 1e-12f)
			return false;
		invLambda[i] = MneeInverse(m, det);
	}

	return true;
}

// Newton step of the chain (vector right hand side)
static bool MneeTridiagonalSolve(const MneeJacobianBlock *blocks, const u_int n,
		const MneeVec2 *rhs, MneeVec2 *dx) {
	MneeMat2 tmp[MNEE_MS_MAX_VERTICES], invLambda[MNEE_MS_MAX_VERTICES];
	if (!MneeTridiagonalInvert(blocks, n, tmp, invLambda))
		return false;

	dx[0] = rhs[0];
	for (u_int i = 1; i < n; ++i) {
		const MneeVec2 back = MneeMatVec(tmp[i], dx[i - 1]);
		dx[i] = MneeVec2{ rhs[i].x - back.x, rhs[i].y - back.y };
	}
	dx[n - 1] = MneeMatVec(invLambda[n - 1], dx[n - 1]);
	for (int i = (int)n - 2; i >= 0; --i) {
		const MneeVec2 back = MneeMatVec(blocks[i].next, dx[i + 1]);
		dx[i] = MneeMatVec(invLambda[i], MneeVec2{ dx[i].x - back.x,
				dx[i].y - back.y });
	}

	return true;
}

// The same verified recursion with a matrix right hand side: a unit block at
// the last vertex yields the (1, N) block of the inverse constraint Jacobian,
// which is the chain's response to a perturbation of the light.
static bool MneeTridiagonalSolveMatrixRhs(const MneeJacobianBlock *blocks,
		const u_int n, MneeMat2 &dxFirst) {
	MneeMat2 tmp[MNEE_MS_MAX_VERTICES], invLambda[MNEE_MS_MAX_VERTICES];
	if (!MneeTridiagonalInvert(blocks, n, tmp, invLambda))
		return false;

	const MneeMat2 zero{ 0.f, 0.f, 0.f, 0.f };
	MneeMat2 d[MNEE_MS_MAX_VERTICES];
	for (u_int i = 0; i < n; ++i)
		d[i] = zero;
	d[n - 1] = MneeMat2{ 1.f, 0.f, 0.f, 1.f };

	for (u_int i = 1; i < n; ++i)
		d[i] = MneeMatSub(d[i], MneeMul(tmp[i], d[i - 1]));
	d[n - 1] = MneeMul(invLambda[n - 1], d[n - 1]);
	for (int i = (int)n - 2; i >= 0; --i)
		d[i] = MneeMul(invLambda[i], MneeMatSub(d[i],
				MneeMul(blocks[i].next, d[i + 1])));

	dxFirst = d[0];
	return true;
}

// dC_last/dy: how the last vertex's constraint reacts to moving the light along
// its (fake) tangent frame, the point emitter branch of Zeltner's
// emitter_interaction_to_vertex (the same construction as MneeLightJacobian).
static MneeMat2 MneeChainLightJacobian(const Point &pPrev, const Point &lightPos,
		const MneeVertex &v, const float etaVertex, const float eps) {
	const MneeMat2 zero{ 0.f, 0.f, 0.f, 0.f };

	Vector dLight = v.p - lightPos;
	const float rl = dLight.Length();
	if (rl < 1e-3f)
		return zero;
	dLight *= 1.f / rl;
	Vector s2, t2;
	MneeCoordinateSystem(dLight, s2, t2);

	const Vector wi = Normalize(pPrev - v.p);
	float eta = etaVertex;
	if (Dot(wi, v.gn) < 0.f)
		eta = 1.f / eta;

	Vector h = wi + eta * Normalize(lightPos - v.p);
	if (eta != 1.f)
		h = -h;
	const float l = h.Length();
	if (l < 1e-6f)
		return zero;
	h *= 1.f / l;
	const float C0x = Dot(v.s, h), C0y = Dot(v.t, h);

	float CP[2][2];
	for (u_int k = 0; k < 2; ++k) {
		const Point lightPosP = lightPos + ((k == 0) ? eps * s2 : eps * t2);
		Vector hP = wi + eta * Normalize(lightPosP - v.p);
		if (eta != 1.f)
			hP = -hP;
		hP *= 1.f / hP.Length();
		CP[k][0] = Dot(v.s, hP);
		CP[k][1] = Dot(v.t, hP);
	}

	return MneeMat2{ (CP[0][0] - C0x) / eps, (CP[1][0] - C0x) / eps,
			(CP[0][1] - C0y) / eps, (CP[1][1] - C0y) / eps };
}

// Initialize a chain vertex from a delta specular hit. Returns false for
// materials the chain solver does not model.
static bool MneeChainVertexInit(MneeChainVertex &cv, const BSDF &bsdf) {
	const MaterialType type = bsdf.GetMaterialType();

	cv.bsdf = bsdf;
	cv.mirrorMat = nullptr;
	cv.glassMat = nullptr;
	cv.specEvent = SPECULAR;
	cv.specFactor = Spectrum(1.f);

	if (type == MIRROR) {
		cv.mirrorMat = dynamic_observer_cast<const MirrorMaterial>(bsdf.GetMaterial());
		if (!cv.mirrorMat)
			return false;
		cv.etaVertex = 1.f;
		cv.specFactor = cv.mirrorMat->GetKr()->GetSpectrumValue(bsdf.hitPoint).Clamp(0.f, 1.f);
		cv.specEvent |= REFLECT;
	} else if (type == GLASS) {
		cv.glassMat = dynamic_observer_cast<const GlassMaterial>(bsdf.GetMaterial());
		if (!cv.glassMat)
			return false;
		// Dispersive glass: the chain uses a single IOR ratio, skip
		if (cv.glassMat->GetCauchyB() &&
				cv.glassMat->GetCauchyB()->GetFloatValue(bsdf.hitPoint) > 0.f)
			return false;

		const float nc = ExtractExteriorIors(bsdf.hitPoint, cv.glassMat->GetExteriorIOR());
		const float nt = ExtractInteriorIors(bsdf.hitPoint, cv.glassMat->GetInteriorIOR());
		if (nt <= 0.f || nc <= 0.f)
			return false;
		cv.etaVertex = nt / nc;
		cv.specEvent |= TRANSMIT;
	} else
		return false;

	MneeInitVertex(cv.v, bsdf, cv.etaVertex);
	return true;
}

// Chain topology: trace the straight ray from x0 toward the light and collect
// the delta specular surfaces it pierces (the reference's mnee_init seeding).
// The first vertex is the shadow ray's blocker, which the caller already has.
// Returns the number of collected vertices; the solver only handles 2 or more
// (a single vertex is the SS solver's job).
static u_int MneeChainDiscover(
		luxrays::IntersectionDeviceRef device, SceneConstRef scene,
		const float time, const Point &x0p, const Point &lightPos,
		const float u4, PathVolumeInfo volInfo,
		const BSDF &firstBsdf, MneeChainVertex *chain, const u_int maxVertices) {
	if (maxVertices < 2)
		return 0;
	if (!MneeChainVertexInit(chain[0], firstBsdf))
		return 0;
	u_int n = 1;

	const Vector dir = Normalize(lightPos - x0p);
	for (u_int guard = 0; (guard < MNEE_MS_MAX_VERTICES) && (n < maxVertices); ++guard) {
		Ray ray(chain[n - 1].bsdf.GetRayOrigin(dir), dir, 0.f,
				numeric_limits<float>::infinity(), time);
		RayHit hit;
		BSDF hitBsdf;
		Spectrum through;
		PathVolumeInfo vol = volInfo;
		if (!scene.Intersect(IntersectionDevicePtr(&device), INDIRECT_RAY, &vol, u4,
				&ray, &hit, &hitBsdf, &through, nullptr, nullptr, false))
			break;      // the ray escaped: the last vertex may see the light

		if (!hitBsdf.IsDelta() || !(hitBsdf.GetEventTypes() & SPECULAR))
			break;      // a non specular surface ends the chain

		MneeChainVertex cv;
		if (!MneeChainVertexInit(cv, hitBsdf))
			break;
		chain[n++] = cv;
	}

	return n;
}

bool PathTracer::MNEEMultiDirectSampling(
		luxrays::IntersectionDeviceRef device,
		SceneConstRef scene,
		const float time,
		const EyePathInfo &pathInfo,
		const luxrays::Spectrum &pathThroughput,
		const BSDF &bsdf,
		LightSourceConstRef light, const float lightPickPdf, const float risScale,
		const luxrays::Ray &shadowRay, const float directPdfW0,
		const luxrays::RayHit &shadowRayHit,
		const BSDF &shadowBsdf, PathVolumeInfo &volInfo,
		const float u1, const float u2, const float u3, const float u4,
		SampleResult *sampleResult) const {
	const Point &x0p = bsdf.hitPoint.p;
	// Light position: the shadow ray maxt has been rewritten by Scene::Intersect
	// to the occluder distance, so recover the light distance from directPdfW
	// (= squared distance to the light, preserved by Illuminate).
	const Point lightPos = shadowRay.o + shadowRay.d * sqrtf(directPdfW0);

	// Every exit of this function is a sample the estimator does not cover, and
	// coverage (not the per-chain value, which an independent analytic model
	// reproduces to 1-3%) is what limits a curved caster. LUX_MNEE_REJ prints one
	// line per rejected attempt so the reasons can be counted: pair it with
	// LUX_MNEE_DEBUG and attempts = accepted + rejected.
	static const bool rejDebug = (getenv("LUX_MNEE_REJ") != nullptr);
	auto rej = [&](const char *why) {
		if (rejDebug)
			printf("MNEE_MS_REJ %s film=%.4g %.4g\n", why, sampleResult->filmX,
					sampleResult->filmY);
		return false;
	};

	//--------------------------------------------------------------------------
	// Chain topology (straight line seed)
	//--------------------------------------------------------------------------
	MneeChainVertex chain[MNEE_MS_MAX_VERTICES];
	const u_int maxVertices = Min(mneeMaxSpecular, MNEE_MS_MAX_VERTICES);
	const u_int n = MneeChainDiscover(device, scene, time, x0p, lightPos, u4,
			volInfo, shadowBsdf, chain, maxVertices);
	if (n < 2)
		return rej("chain<2");

	//--------------------------------------------------------------------------
	// Newton solve on the whole chain (block tridiagonal step, line search with
	// re-projection onto the shapes).
	//--------------------------------------------------------------------------
	MneeJacobianBlock blocks[MNEE_MS_MAX_VERTICES];
	MneeVec2 residual[MNEE_MS_MAX_VERTICES];
	float maxResidual = 0.f;
	bool solved = false;
	float beta = 1.f;
	u_int iteration = 0;
	const char *failWhy = "iterations";

	// Per-iteration trace (LUX_MNEE_ITER): shows whether the Newton stalls at a
	// fixed residual, creeps down, or oscillates. The outer caustic of a curved
	// caster is where the chain solve fails, and the failure mode decides what to
	// change (step size policy vs the walk parameterization itself).
	static const bool iterDebug = (getenv("LUX_MNEE_ITER") != nullptr);

	while (iteration < mneeMaxIterations) {
		if (!MneeChainJacobian(device, scene, time, x0p, lightPos, chain, n,
				blocks, residual, maxResidual)) {
			failWhy = "jacobian";
			break;
		}
		if (iterDebug) {
			printf("MNEE_MS_IT film=%.4g %.4g it=%u res=%.6g beta=%.4g "
					"x0=%.9g %.9g %.9g",
					sampleResult->filmX, sampleResult->filmY, iteration,
					maxResidual, beta, x0p.x, x0p.y, x0p.z);
			for (u_int k = 0; k < n; ++k)
				printf(" | x%u=%.9g %.9g %.9g C=(%.4g %.4g)", k + 1,
						chain[k].v.p.x, chain[k].v.p.y, chain[k].v.p.z,
						residual[k].x, residual[k].y);
			printf("\n");
		}
		if (maxResidual < 1e-5f) {
			solved = true;
			break;
		}

		MneeVec2 dx[MNEE_MS_MAX_VERTICES];
		if (!MneeTridiagonalSolve(blocks, n, residual, dx)) {
			failWhy = "tridiagonal";
			break;
		}

		bool stepAccepted = false;
		while (beta > 1e-2f) {
			MneeChainVertex trial[MNEE_MS_MAX_VERTICES];
			bool projected = true;
			for (u_int i = 0; i < n; ++i) {
				const float eps = Max(1e-5f, 1e-4f * Distance(x0p, chain[i].v.p));
				const Point pProp = chain[i].v.p - beta *
						(chain[i].v.dpdu * dx[i].x + chain[i].v.dpdv * dx[i].y);
				BSDF propBsdf;
				if (!MneeReproject(device, scene, time, pProp, chain[i].v.gn,
						.5f * eps, propBsdf) ||
						(propBsdf.GetMaterialType() != chain[i].bsdf.GetMaterialType())) {
					projected = false;
					break;
				}
				// Keep the material data and the constraint mode of the seed:
				// only the surface position and differentials move.
				trial[i] = chain[i];
				MneeInitVertex(trial[i].v, propBsdf, chain[i].etaVertex);
				trial[i].bsdf = propBsdf;
			}

			if (projected) {
				MneeVec2 trialRes[MNEE_MS_MAX_VERTICES];
				float trialMax = 0.f;
				if (MneeChainResiduals(x0p, lightPos, trial, n, trialRes, trialMax) &&
						(trialMax < maxResidual)) {
					for (u_int i = 0; i < n; ++i)
						chain[i] = trial[i];
					beta = Min(1.f, 2.f * beta);
					stepAccepted = true;
					break;
				}
			}

			beta *= .5f;
			++iteration;
		}
		if (!stepAccepted) {
			failWhy = "no-step";
			break;
		}
		++iteration;
	}

	if (!solved)
		return rej(failWhy);

	//--------------------------------------------------------------------------
	// Post-solve validity: each vertex must have solved the mode its material
	// implies (a mirror reflects, so the two segments stay on one side of the
	// surface; glass transmits, so they are on opposite sides), and every
	// specular factor must be valid at the solved configuration.
	//--------------------------------------------------------------------------
	Point pts[MNEE_MS_MAX_VERTICES + 2];
	pts[0] = x0p;
	for (u_int i = 0; i < n; ++i)
		pts[i + 1] = chain[i].v.p;
	pts[n + 1] = lightPos;

	Spectrum specProduct(1.f);
	bool plainHalfVector = true;
	for (u_int i = 0; i < n; ++i) {
		const Vector wi = Normalize(pts[i] - pts[i + 1]);
		const Vector wo = Normalize(pts[i + 2] - pts[i + 1]);
		const float cosI = Dot(chain[i].v.gn, wi);
		const float cosO = Dot(chain[i].v.gn, wo);

		if (chain[i].etaVertex == 1.f) {
			if (cosI * cosO < 0.f)
				return rej("mirror-side");
			specProduct *= chain[i].specFactor;
		} else {
			if (cosI * cosO > 0.f)
				return rej("dielectric-same-side");

			const Spectrum kt = chain[i].glassMat->GetKt()->
					GetSpectrumValue(chain[i].bsdf.hitPoint).Clamp(0.f, 1.f);
			const float nc = ExtractExteriorIors(chain[i].bsdf.hitPoint,
					chain[i].glassMat->GetExteriorIOR());
			const float nt = ExtractInteriorIors(chain[i].bsdf.hitPoint,
					chain[i].glassMat->GetInteriorIOR());
			const Vector localFixedDir = chain[i].bsdf.GetFrame().ToLocal(wi);
			Vector localSampledDir;
			const Spectrum trans = GlassMaterial::EvalSpecularTransmission(
					chain[i].bsdf.hitPoint, localFixedDir, 0.f, kt, nc, nt, 0.f,
					&localSampledDir);
			if (trans.Black())
				return rej("tir");
			specProduct *= trans;
			plainHalfVector = false;
		}
	}
	if (specProduct.Black())
		return rej("spec-black");

	//--------------------------------------------------------------------------
	// Last segment xN -> y: Illuminate at the last vertex, then check visibility
	//--------------------------------------------------------------------------
	PathVolumeInfo volLast = volInfo;
	for (u_int i = 0; i < n; ++i)
		volLast.Update(chain[i].specEvent, chain[i].bsdf);

	Ray shadowRay2;
	float directPdfW2;
	const Spectrum lightRadiance2 = light.Illuminate(scene, chain[n - 1].bsdf, time,
			u1, u2, u3, shadowRay2, directPdfW2);
	if (lightRadiance2.Black())
		return rej("light-black");
	verify(!isnan(directPdfW2) && !isinf(directPdfW2));

	RayHit occlHit;
	BSDF occlBsdf;
	Spectrum segThroughput;
	PathVolumeInfo volTrace = volLast;
	if (scene.Intersect(IntersectionDevicePtr(&device), SHADOW_RAY, &volTrace, u4,
			&shadowRay2, &occlHit, &occlBsdf, &segThroughput, nullptr, nullptr, true))
		return rej("occluded");

	//--------------------------------------------------------------------------
	// Chain geometric term: dw0_dx1 * |det(dx_1 / dy)| with
	//   dx_1 / dy = (A^-1)_{1,N} * dC_N/dy
	// (A = the block tridiagonal constraint Jacobian; only the last constraint
	// involves the light). For N = 1 this is the single vertex term.
	//--------------------------------------------------------------------------
	MneeJacobianBlock geoBlocks[MNEE_MS_MAX_VERTICES];
	MneeVec2 geoResidual[MNEE_MS_MAX_VERTICES];
	float geoMax = 0.f;
	if (!MneeChainJacobian(device, scene, time, x0p, lightPos, chain, n,
			geoBlocks, geoResidual, geoMax))
		return rej("geo-jacobian");

	MneeMat2 dxFirst;
	if (!MneeTridiagonalSolveMatrixRhs(geoBlocks, n, dxFirst))
		return rej("geo-tridiagonal");

	const float epsLight = Max(1e-5f, 1e-4f * Distance(x0p, chain[n - 1].v.p));
	const MneeMat2 lightJac = MneeChainLightJacobian(pts[n - 1], lightPos,
			chain[n - 1].v, chain[n - 1].etaVertex, epsLight);
	const MneeMat2 dxDy = MneeMul(dxFirst, lightJac);

	const Vector d01 = x0p - chain[0].v.p;
	const float r01sq = d01.LengthSquared();
	if (r01sq < 1e-6f)
		return false;
	const float dw0Dx1 = fabsf(Dot(d01, chain[0].v.gn)) / (sqrtf(r01sq) * r01sq);
	const float geometricTerm = dw0Dx1 * fabsf(MneeDet(dxDy));
	if (geometricTerm <= 0.f || isnan(geometricTerm) || isinf(geometricTerm))
		return false;

	//--------------------------------------------------------------------------
	// Contribution. The r12^2 (Illuminate directPdfW) measure factor belongs to
	// the light weight only when every vertex solved the plain half-vector
	// (eta == 1, a pure reflection chain); with any dielectric vertex the
	// analytic geometric term carries the full conversion (see the single
	// vertex solver for the validation of this rule).
	//--------------------------------------------------------------------------
	BSDFEvent receiverEvent;
	float bsdfPdfW;
	const Spectrum bsdfEval0 = bsdf.Evaluate(Normalize(chain[0].v.p - x0p),
			&receiverEvent, &bsdfPdfW);
	if (bsdfEval0.Black())
		return false;

	PathDepthInfo mneeDepthInfo = pathInfo.depth;
	mneeDepthInfo.IncDepths(receiverEvent);
	for (u_int i = 0; i < n; ++i)
		mneeDepthInfo.IncDepths(chain[i].specEvent);

	const Spectrum lightWeight = lightRadiance2 *
			((plainHalfVector ? directPdfW2 : 1.f) * risScale / lightPickPdf);
	const Spectrum incomingRadiance = bsdfEval0 * specProduct * geometricTerm *
			lightWeight * segThroughput;
	verify(!incomingRadiance.IsNaN() && !incomingRadiance.IsInf());

	if (MneeDebugEnabled()) {
		const float dbgR12 = Distance(lightPos, chain[n - 1].v.p);
		// The analytic per-point radiance a light transport estimator measures
		// for this chain: every specular factor times the light's radiance at
		// the last vertex (see the single vertex debug print).
		const float dbgTruth = bsdfEval0.Filter() * specProduct.Filter() *
				lightRadiance2.Filter() / (dbgR12 * dbgR12);
		printf("MNEE_MS_DBG film=%.4g %.4g n=%u x0=%.9g %.9g %.9g",
				sampleResult->filmX, sampleResult->filmY, n,
				x0p.x, x0p.y, x0p.z);
		for (u_int i = 0; i < n; ++i) {
			// The eta side rule keys off Dot(wi, gn), so the debug output has to
			// show the geometric normal itself: a mesh whose geometric normals
			// point against the surface's outward direction silently inverts
			// every vertex's constraint IOR (see dev-tools/mnee_design.md 4e).
			const Vector dbgWi = Normalize(pts[i] - pts[i + 1]);
			const Vector dbgWo = Normalize(pts[i + 2] - pts[i + 1]);
			printf(" | x%u=%.9g %.9g %.9g gn=%.4g %.4g %.4g eta=%.9g "
					"dWiGn=%.4g cWiGn=%.4g cWoGn=%.4g", i + 1,
					chain[i].v.p.x, chain[i].v.p.y, chain[i].v.p.z,
					chain[i].v.gn.x, chain[i].v.gn.y, chain[i].v.gn.z,
					chain[i].etaVertex, Dot(dbgWi, chain[i].v.gn),
					fabsf(Dot(dbgWi, chain[i].v.gn)),
					fabsf(Dot(dbgWo, chain[i].v.gn)));
		}
		printf(" | y=%.9g %.9g %.9g residual=%.3e spec=%.9g G=%.9g in=%.9g "
				"truth=%.9g\n", lightPos.x, lightPos.y, lightPos.z,
				maxResidual, specProduct.Filter(), geometricTerm,
				incomingRadiance.Filter(), dbgTruth);
		fflush(stdout);
	}

	sampleResult->AddDirectLight(light.GetID(), chain[n - 1].specEvent,
			pathThroughput, incomingRadiance, 1.f);

	return true;
}
// vim: autoindent noexpandtab tabstop=4 shiftwidth=4
