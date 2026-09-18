/***************************************************************************
 * Copyright 1998-2020 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxCoreRender.                                   *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

// Hair BSDF after Chiang et al. 2019 "A Practical and Controllable Hair and
// Fur Model". Ported from pbrt-v3's HairBSDF (BSD-3 licensed, (c) Matt Pharr,
// Greg Humphreys, Wenzel Jakob) and adapted to the LuxCore material API:
//   - the fiber tangent comes from the strands tessellation via vertex AOV
//     layers (object space, transformed to the shading frame per hit)
//   - the azimuthal hit offset h is computed geometrically from the shading
//     normal and the outgoing direction (works for true cylinder geometry)
//   - LuxCore's "f * cos" Evaluate convention means we return pbrt's fsum
//     without the final 1/cos division

#include <array>
#include <numeric>

#include "slg/materials/hairmat.h"
#include "slg/textures/fresnel/fresneltexture.h"
#include "luxrays/utils/utils.h"

using namespace luxrays;
using namespace slg;

namespace {

constexpr int pMax = 3;
constexpr float SqrtPiOver8 = 0.626657069f; // sqrt(Pi / 8)

inline float SafeSqrt(const float x) { return sqrtf(Max(0.f, x)); }
inline float SafeASin(const float x) { return asinf(Clamp(x, -1.f, 1.f)); }
inline float Fract(const float x) { return x - floorf(x); }

inline float I0(float x) {
	float val = 0.f;
	float x2i = 1.f;
	int64_t ifact = 1;
	int i4 = 1;
	// I0(x) ~= Sum_i x^(2i) / (4^i (i!)^2)
	for (int i = 0; i < 10; ++i) {
		if (i > 1) ifact *= i;
		val += x2i / (i4 * Sqr(ifact));
		x2i *= x * x;
		i4 *= 4;
	}
	return val;
}

inline float LogI0(float x) {
	if (x > 12.f)
		return x + 0.5f * (-logf(2.f * M_PI) + logf(1.f / x) + 1.f / (8.f * x));
	else
		return logf(I0(x));
}

inline float Mp(float cosThetaI, float cosThetaO, float sinThetaI,
		float sinThetaO, float v) {
	const float a = cosThetaI * cosThetaO / v;
	const float b = sinThetaI * sinThetaO / v;
	const float mp = (v <= .1f)
		? expf(LogI0(a) - b - 1.f / v + 0.6931f + logf(1.f / (2.f * v)))
		: (expf(-b) * I0(a)) / (sinhf(1.f / v) * 2.f * v);
	return mp;
}

inline float Phi(const int p, const float gammaO, const float gammaT) {
	return 2.f * p * gammaT - 2.f * gammaO + p * M_PI;
}

inline float Logistic(float x, const float s) {
	x = fabsf(x);
	return expf(-x / s) / (s * Sqr(1.f + expf(-x / s)));
}

inline float LogisticCDF(const float x, const float s) {
	return 1.f / (1.f + expf(-x / s));
}

inline float TrimmedLogistic(const float x, const float s, const float a, const float b) {
	return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
}

inline float Np(float phi, const int p, const float s,
		const float gammaO, const float gammaT) {
	float dphi = phi - Phi(p, gammaO, gammaT);
	// Remap dphi to [-pi, pi]
	while (dphi > M_PI) dphi -= 2.f * M_PI;
	while (dphi < -M_PI) dphi += 2.f * M_PI;
	return TrimmedLogistic(dphi, s, -M_PI, M_PI);
}

inline float SampleTrimmedLogistic(const float u, const float s,
		const float a, const float b) {
	const float k = LogisticCDF(b, s) - LogisticCDF(a, s);
	const float x = -s * logf(1.f / (u * k + LogisticCDF(a, s)) - 1.f);
	return Clamp(x, a, b);
}

// Attenuation A_p for each scattering order p = 0..pMax
inline void Ap(const float cosThetaO, const float eta, const float h,
		const Spectrum &T, Spectrum *ap) {
	// p = 0 attenuation at initial cylinder intersection
	const float cosGammaO = SafeSqrt(1.f - h * h);
	const float cosTheta = cosThetaO * cosGammaO;
	const float f = FresnelTexture::CauchyEvaluate(eta, cosTheta);
	ap[0] = f;

	// p = 1 attenuation
	ap[1] = Sqr(1.f - f) * T;

	// Attenuation up to pMax
	for (int p = 2; p < pMax; ++p)
		ap[p] = ap[p - 1] * T * f;

	// Remaining orders of scattering
	ap[pMax] = ap[pMax - 1] * f * T / (Spectrum(1.f) - T * f);
}

// Rotate (sinThetaO, cosThetaO) by the lobe-dependent scale tilt
inline void TiltThetaO(const int p, const float *sin2k, const float *cos2k,
		const float sinThetaO, const float cosThetaO,
		float *sinThetaOp, float *cosThetaOp) {
	if (p == 0) {
		*sinThetaOp = sinThetaO * cos2k[1] - cosThetaO * sin2k[1];
		*cosThetaOp = cosThetaO * cos2k[1] + sinThetaO * sin2k[1];
	} else if (p == 1) {
		*sinThetaOp = sinThetaO * cos2k[0] + cosThetaO * sin2k[0];
		*cosThetaOp = cosThetaO * cos2k[0] - sinThetaO * sin2k[0];
	} else if (p == 2) {
		*sinThetaOp = sinThetaO * cos2k[2] + cosThetaO * sin2k[2];
		*cosThetaOp = cosThetaO * cos2k[2] - sinThetaO * sin2k[2];
	} else {
		*sinThetaOp = sinThetaO;
		*cosThetaOp = cosThetaO;
	}
	*cosThetaOp = fabsf(*cosThetaOp);
}

// Split one stratified float into two (pbrt-style demux) so the BSDF gets
// its four sample dimensions from the three values LuxCore supplies.
inline void DemuxFloat(const float u, float *a, float *b) {
	*a = Fract(u * 65536.f);
	*b = Fract(u * 4294967296.f);
}

// Shared per-hit setup for the hair frame basis (tangent + radial dir).
// Returns e1 (radial) and e2 = cross(tangent, e1).
inline void HairFrame(const Vector &T, Vector *e1, Vector *e2) {
	Vector nPerp = Vector(0.f, 0.f, 1.f) - T * T.z;
	if (nPerp.LengthSquared() < 1e-12f)
		nPerp = Vector(1.f, 0.f, 0.f) - T * T.x;
	*e1 = Normalize(nPerp);
	*e2 = Cross(T, *e1);
}

} // anonymous namespace

//------------------------------------------------------------------------------
// HairMaterial
//------------------------------------------------------------------------------

HairMaterial::HairMaterial(TextureConstPtr frontTransp, TextureConstPtr backTransp,
		TextureConstPtr emitted, TextureConstPtr bump,
		TextureConstPtr sigmaATex, TextureConstPtr colorTex,
		TextureConstPtr eumelaninTex, TextureConstPtr pheomelaninTex,
		TextureConstPtr etaTex, TextureConstPtr betaMTex, TextureConstPtr betaNTex,
		TextureConstPtr alphaTex) :
		Material(frontTransp, backTransp, emitted, bump),
		sigmaA(sigmaATex), color(colorTex),
		eumelanin(eumelaninTex), pheomelanin(pheomelaninTex),
		eta(etaTex), betaM(betaMTex), betaN(betaNTex), alpha(alphaTex) {
}

Spectrum HairMaterial::EvalSigmaA(const HitPoint &hitPoint, const float betaN) const {
	if (sigmaA)
		return sigmaA->GetSpectrumValue(hitPoint).Clamp();

	if (color) {
		const Spectrum c = color->GetSpectrumValue(hitPoint).Clamp(0.f, 1.f);
		// SigmaAFromReflectance (Chiang 2019)
		const float denom = 5.969f - 0.215f * betaN + 2.532f * Sqr(betaN) -
				10.73f * (betaN * betaN * betaN) + 5.574f * powf(betaN, 4.f) +
				0.245f * powf(betaN, 5.f);
		Spectrum sa;
		for (u_int i = 0; i < COLOR_SAMPLES; ++i)
			sa.c[i] = Sqr(logf(Max(c.c[i], 1e-5f)) / denom);
		return sa;
	}

	// Melanin concentration model (default: brown hair 1.3 / 0)
	const float ce = eumelanin ? Max(0.f, eumelanin->GetFloatValue(hitPoint)) : 1.3f;
	const float cp = pheomelanin ? Max(0.f, pheomelanin->GetFloatValue(hitPoint)) : 0.f;
	const float eu[3] = {0.419f, 0.697f, 1.37f};
	const float ph[3] = {0.187f, 0.4f, 1.05f};
	float sa[3];
	for (int i = 0; i < 3; ++i)
		sa[i] = ce * eu[i] + cp * ph[i];
	return Spectrum(sa);
}

bool HairMaterial::SetupContext(const HitPoint &hitPoint,
		const Vector &localEyeDir, HairContext *ctx) const {
	// Strand tangent: vertex AOV (object space) -> world -> local frame
	const float tx = hitPoint.GetVertexAOV(HAIR_TANGENT_X_DATA_INDEX);
	const float ty = hitPoint.GetVertexAOV(HAIR_TANGENT_Y_DATA_INDEX);
	const float tz = hitPoint.GetVertexAOV(HAIR_TANGENT_Z_DATA_INDEX);
	Vector tObj(tx, ty, tz);
	if (tObj.LengthSquared() < 1e-12f) {
		// No strand data (e.g. material used on a regular mesh): fall back
		// to dpdv which is the fiber direction for strand-like UV layouts.
		tObj = hitPoint.dpdv;
	}
	const Vector tWorld = Normalize(hitPoint.localToWorld * tObj);
	ctx->tangent = Normalize(hitPoint.GetFrame().ToLocal(tWorld));

	const float bm = Clamp(betaM->GetFloatValue(hitPoint), 1e-2f, 1.f);
	const float bn = Clamp(betaN->GetFloatValue(hitPoint), 1e-2f, 1.f);
	const float a = alpha->GetFloatValue(hitPoint);
	ctx->eta = eta->GetFloatValue(hitPoint);

	// Longitudinal variance from beta_m
	ctx->v[0] = Sqr(0.726f * bm + 0.812f * Sqr(bm) + 3.7f * powf(bm, 20.f));
	ctx->v[1] = .25f * ctx->v[0];
	ctx->v[2] = 4.f * ctx->v[0];
	ctx->v[3] = ctx->v[2];

	// Azimuthal logistic scale from beta_n
	ctx->s = SqrtPiOver8 *
			(0.265f * bn + 1.194f * Sqr(bn) + 5.372f * powf(bn, 22.f));

	// Scale tilt (alpha) terms
	ctx->sin2kAlpha[0] = sinf(Radians(a));
	ctx->cos2kAlpha[0] = SafeSqrt(1.f - Sqr(ctx->sin2kAlpha[0]));
	for (int i = 1; i < 3; ++i) {
		ctx->sin2kAlpha[i] = 2.f * ctx->cos2kAlpha[i - 1] * ctx->sin2kAlpha[i - 1];
		ctx->cos2kAlpha[i] = Sqr(ctx->cos2kAlpha[i - 1]) - Sqr(ctx->sin2kAlpha[i - 1]);
	}

	// Azimuthal offset h: sine of the angle between the hit radial
	// direction and the ray-azimuth plane. e1_wo is the azimuthal direction
	// of the outgoing ray projected perpendicular to the fiber.
	Vector e1, e2;
	HairFrame(ctx->tangent, &e1, &e2);
	Vector e1wo = localEyeDir - ctx->tangent * Dot(localEyeDir, ctx->tangent);
	if (e1wo.LengthSquared() < 1e-12f)
		e1wo = e1; // outgoing ray parallel to the fiber: h = 0
	else
		e1wo = Normalize(e1wo);
	const Vector e2wo = Cross(ctx->tangent, e1wo);
	const float h = Clamp(Dot(e1, e2wo), -1.f, 1.f);
	ctx->gammaO = SafeASin(h);

	ctx->sigma_a = EvalSigmaA(hitPoint, bn);
	return true;
}

Spectrum HairMaterial::f(const HitPoint &hitPoint, const HairContext &ctx,
		const Vector &localLightDir, const Vector &localEyeDir) const {
	const Vector &T = ctx.tangent;
	Vector e1, e2;
	HairFrame(T, &e1, &e2);

	// Hair coordinate terms for wo (outgoing/eye)
	const float sinThetaO = Dot(localEyeDir, T);
	const float cosThetaO = SafeSqrt(1.f - Sqr(sinThetaO));
	const float phiO = atan2f(Dot(localEyeDir, e2), Dot(localEyeDir, e1));

	// Hair coordinate terms for wi (light)
	const float sinThetaI = Dot(localLightDir, T);
	const float cosThetaI = SafeSqrt(1.f - Sqr(sinThetaI));
	const float phiI = atan2f(Dot(localLightDir, e2), Dot(localLightDir, e1));

	// cos(thetaT) for refracted ray
	const float sinThetaT = sinThetaO / ctx.eta;
	const float cosThetaT = SafeSqrt(1.f - Sqr(sinThetaT));

	// gammaT for refracted ray
	const float h = sinf(ctx.gammaO);
	const float etap = sqrtf(ctx.eta * ctx.eta - Sqr(sinThetaO)) / cosThetaO;
	const float sinGammaT = h / etap;
	const float cosGammaT = SafeSqrt(1.f - Sqr(sinGammaT));
	const float gammaT = SafeASin(sinGammaT);

	// Transmittance T of a single path through the cylinder
	const Spectrum Tt = Exp(-ctx.sigma_a * (2.f * cosGammaT / cosThetaT));

	// Evaluate hair BSDF
	const float phi = phiI - phiO;
	Spectrum ap[pMax + 1];
	Ap(cosThetaO, ctx.eta, h, Tt, ap);

	Spectrum fsum(0.f);
	for (int p = 0; p < pMax; ++p) {
		float sinThetaOp, cosThetaOp;
		TiltThetaO(p, ctx.sin2kAlpha, ctx.cos2kAlpha,
				sinThetaO, cosThetaO, &sinThetaOp, &cosThetaOp);
		fsum += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, ctx.v[p]) *
				ap[p] * Np(phi, p, ctx.s, ctx.gammaO, gammaT);
	}

	// Remaining terms after pMax
	fsum += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, ctx.v[pMax]) *
			ap[pMax] / (2.f * M_PI);
	// LuxCore wants f * cos: pbrt's BSDF value is fsum/|cosThetaI| so the
	// cosine-multiplied value is exactly fsum — no division needed.
	return fsum;
}

void HairMaterial::ComputeApPdf(const HairContext &ctx, const float cosThetaO,
		float *apPdf) const {
	const float sinThetaO = SafeSqrt(1.f - cosThetaO * cosThetaO);

	const float sinThetaT = sinThetaO / ctx.eta;
	const float cosThetaT = SafeSqrt(1.f - Sqr(sinThetaT));

	const float h = sinf(ctx.gammaO);
	const float etap = sqrtf(ctx.eta * ctx.eta - Sqr(sinThetaO)) / cosThetaO;
	const float sinGammaT = h / etap;
	const float cosGammaT = SafeSqrt(1.f - Sqr(sinGammaT));

	const Spectrum T = Exp(-ctx.sigma_a * (2.f * cosGammaT / cosThetaT));
	Spectrum ap[pMax + 1];
	Ap(cosThetaO, ctx.eta, h, T, ap);

	float sumY = 0.f;
	for (int i = 0; i <= pMax; ++i)
		sumY += ap[i].Y();
	for (int i = 0; i <= pMax; ++i)
		apPdf[i] = (sumY > 0.f) ? (ap[i].Y() / sumY) : 0.f;
}

Spectrum HairMaterial::Albedo(const HitPoint &hitPoint) const {
	HairContext ctx;
	if (!SetupContext(hitPoint, Vector(0.f, 0.f, 1.f), &ctx))
		return Spectrum();
	float apPdf[pMax + 1];
	ComputeApPdf(ctx, 1.f, apPdf);
	float sum = 0.f;
	for (int i = 0; i <= pMax; ++i)
		sum += apPdf[i];
	return Spectrum(Clamp(sum * 0.5f, 0.f, 1.f));
}

Spectrum HairMaterial::Evaluate(const HitPoint &hitPoint,
		const Vector &localLightDir, const Vector &localEyeDir,
		BSDFEvent *event, float *directPdfW, float *reversePdfW) const {
	HairContext ctx;
	if (!SetupContext(hitPoint, localEyeDir, &ctx)) {
		*event = NONE;
		return Spectrum();
	}

	const Spectrum result = f(hitPoint, ctx, localLightDir, localEyeDir);
	if (result.Black()) {
		*event = NONE;
		return Spectrum();
	}

	*event = GLOSSY | REFLECT | TRANSMIT;

	if (directPdfW) {
		float apPdf[pMax + 1];
		const float sinThetaO = Dot(localEyeDir, ctx.tangent);
		const float cosThetaO = SafeSqrt(1.f - Sqr(sinThetaO));
		const float sinThetaI = Dot(localLightDir, ctx.tangent);
		const float cosThetaI = SafeSqrt(1.f - Sqr(sinThetaI));
		ComputeApPdf(ctx, cosThetaO, apPdf);

		Vector e1, e2;
		HairFrame(ctx.tangent, &e1, &e2);
		const float phiO = atan2f(Dot(localEyeDir, e2), Dot(localEyeDir, e1));
		const float phiI = atan2f(Dot(localLightDir, e2), Dot(localLightDir, e1));
		const float phi = phiI - phiO;

		const float h = sinf(ctx.gammaO);
		const float etap = sqrtf(ctx.eta * ctx.eta - Sqr(sinThetaO)) / cosThetaO;
		const float sinGammaT = h / etap;
		const float gammaT = SafeASin(sinGammaT);

		float pdf = 0.f;
		for (int p = 0; p < pMax; ++p) {
			float sinThetaOp, cosThetaOp;
			TiltThetaO(p, ctx.sin2kAlpha, ctx.cos2kAlpha,
					sinThetaO, cosThetaO, &sinThetaOp, &cosThetaOp);
			pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, ctx.v[p]) *
					apPdf[p] * Np(phi, p, ctx.s, ctx.gammaO, gammaT);
		}
		pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, ctx.v[pMax]) *
				apPdf[pMax] * (1.f / (2.f * M_PI));
		*directPdfW = pdf;
	}
	if (reversePdfW)
		*reversePdfW = *directPdfW; // the hair BSDF is reciprocal here

	return result;
}

Spectrum HairMaterial::Sample(const HitPoint &hitPoint,
		const Vector &localFixedDir, Vector *localSampledDir,
		const float u0, const float u1, const float passThroughEvent,
		float *pdfW, BSDFEvent *event) const {
	HairContext ctx;
	if (!SetupContext(hitPoint, localFixedDir, &ctx)) {
		*event = NONE;
		return Spectrum();
	}

	const Vector &T = ctx.tangent;
	Vector e1, e2;
	HairFrame(T, &e1, &e2);

	const float sinThetaO = Dot(localFixedDir, T);
	const float cosThetaO = SafeSqrt(1.f - Sqr(sinThetaO));
	const float phiO = atan2f(Dot(localFixedDir, e2), Dot(localFixedDir, e1));

	// Four sample dimensions from the three supplied values
	float u[2][2];
	DemuxFloat(u0, &u[0][0], &u[0][1]);
	DemuxFloat(u1, &u[1][0], &u[1][1]);
	(void)passThroughEvent;

	// Choose the lobe p from the Ap pdf
	float apPdf[pMax + 1];
	ComputeApPdf(ctx, cosThetaO, apPdf);
	int p;
	float sel = u[0][0];
	for (p = 0; p < pMax; ++p) {
		if (sel < apPdf[p]) break;
		sel -= apPdf[p];
	}

	float sinThetaOp, cosThetaOp;
	TiltThetaO(p, ctx.sin2kAlpha, ctx.cos2kAlpha,
			sinThetaO, cosThetaO, &sinThetaOp, &cosThetaOp);

	// Sample M_p -> thetaI
	const float u10 = Max(u[1][0], 1e-5f);
	const float cosTheta =
			1.f + ctx.v[p] * logf(u10 + (1.f - u10) * expf(-2.f / ctx.v[p]));
	const float sinTheta = SafeSqrt(1.f - Sqr(cosTheta));
	const float cosPhi = cosf(2.f * M_PI * u[1][1]);
	const float sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhi * cosThetaOp;
	const float cosThetaI = SafeSqrt(1.f - Sqr(sinThetaI));

	// Sample N_p -> dphi
	const float h = sinf(ctx.gammaO);
	const float etap = sqrtf(ctx.eta * ctx.eta - Sqr(sinThetaO)) / cosThetaO;
	const float sinGammaT = h / etap;
	const float gammaT = SafeASin(sinGammaT);
	float dphi;
	if (p < pMax)
		dphi = Phi(p, ctx.gammaO, gammaT) +
				SampleTrimmedLogistic(u[0][1], ctx.s, -M_PI, M_PI);
	else
		dphi = 2.f * M_PI * u[0][1];

	const float phiI = phiO + dphi;
	*localSampledDir = sinThetaI * T +
			cosThetaI * cosf(phiI) * e1 + cosThetaI * sinf(phiI) * e2;

	// Pdf for the sampled direction
	float pdf = 0.f;
	for (int pp = 0; pp < pMax; ++pp) {
		float stp, ctp;
		TiltThetaO(pp, ctx.sin2kAlpha, ctx.cos2kAlpha,
				sinThetaO, cosThetaO, &stp, &ctp);
		pdf += Mp(cosThetaI, ctp, sinThetaI, stp, ctx.v[pp]) *
				apPdf[pp] * Np(dphi, pp, ctx.s, ctx.gammaO, gammaT);
	}
	pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, ctx.v[pMax]) *
			apPdf[pMax] * (1.f / (2.f * M_PI));
	*pdfW = pdf;
	if (pdf <= 0.f)
		return Spectrum();

	*event = (p == 0) ? (GLOSSY | REFLECT) : (GLOSSY | TRANSMIT);

	return f(hitPoint, ctx, *localSampledDir, localFixedDir) / pdf;
}

void HairMaterial::Pdf(const HitPoint &hitPoint,
		const Vector &localLightDir, const Vector &localEyeDir,
		float *directPdfW, float *reversePdfW) const {
	HairContext ctx;
	if (!SetupContext(hitPoint, localEyeDir, &ctx)) {
		if (directPdfW) *directPdfW = 0.f;
		if (reversePdfW) *reversePdfW = 0.f;
		return;
	}

	const Vector &T = ctx.tangent;
	Vector e1, e2;
	HairFrame(T, &e1, &e2);

	const float sinThetaO = Dot(localEyeDir, T);
	const float cosThetaO = SafeSqrt(1.f - Sqr(sinThetaO));
	const float phiO = atan2f(Dot(localEyeDir, e2), Dot(localEyeDir, e1));
	const float sinThetaI = Dot(localLightDir, T);
	const float cosThetaI = SafeSqrt(1.f - Sqr(sinThetaI));
	const float phiI = atan2f(Dot(localLightDir, e2), Dot(localLightDir, e1));

	const float h = sinf(ctx.gammaO);
	const float etap = sqrtf(ctx.eta * ctx.eta - Sqr(sinThetaO)) / cosThetaO;
	const float sinGammaT = h / etap;
	const float gammaT = SafeASin(sinGammaT);

	float apPdf[pMax + 1];
	ComputeApPdf(ctx, cosThetaO, apPdf);

	const float phi = phiI - phiO;
	float pdf = 0.f;
	for (int p = 0; p < pMax; ++p) {
		float sinThetaOp, cosThetaOp;
		TiltThetaO(p, ctx.sin2kAlpha, ctx.cos2kAlpha,
				sinThetaO, cosThetaO, &sinThetaOp, &cosThetaOp);
		pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, ctx.v[p]) *
				apPdf[p] * Np(phi, p, ctx.s, ctx.gammaO, gammaT);
	}
	pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, ctx.v[pMax]) *
			apPdf[pMax] * (1.f / (2.f * M_PI));

	if (directPdfW) *directPdfW = pdf;
	if (reversePdfW) *reversePdfW = pdf;
}

void HairMaterial::AddReferencedTextures(std::unordered_set<const Texture *> &referencedTexs) const {
	Material::AddReferencedTextures(referencedTexs);
	if (sigmaA) sigmaA->AddReferencedTextures(referencedTexs);
	if (color) color->AddReferencedTextures(referencedTexs);
	if (eumelanin) eumelanin->AddReferencedTextures(referencedTexs);
	if (pheomelanin) pheomelanin->AddReferencedTextures(referencedTexs);
	eta->AddReferencedTextures(referencedTexs);
	betaM->AddReferencedTextures(referencedTexs);
	betaN->AddReferencedTextures(referencedTexs);
	alpha->AddReferencedTextures(referencedTexs);
}

void HairMaterial::UpdateTextureReferences(TextureConstRef oldTex, TextureRef newTex) {
	Material::UpdateTextureReferences(oldTex, newTex);
	if (sigmaA == &oldTex) sigmaA = &newTex;
	if (color == &oldTex) color = &newTex;
	if (eumelanin == &oldTex) eumelanin = &newTex;
	if (pheomelanin == &oldTex) pheomelanin = &newTex;
	if (eta == &oldTex) eta = &newTex;
	if (betaM == &oldTex) betaM = &newTex;
	if (betaN == &oldTex) betaN = &newTex;
	if (alpha == &oldTex) alpha = &newTex;
}

PropertiesUPtr HairMaterial::ToProperties(const ImageMapCache &imgMapCache,
		const bool useRealFileName) const {
	PropertiesUPtr props = std::make_unique<Properties>();

	props->Set(Material::ToProperties(imgMapCache, useRealFileName));
	props->Set(Property("type")("hairmat"));
	if (sigmaA) props->Set(Property("sigma_a")(sigmaA->GetSDLValue()));
	if (color) props->Set(Property("color")(color->GetSDLValue()));
	if (eumelanin) props->Set(Property("eumelanin")(eumelanin->GetSDLValue()));
	if (pheomelanin) props->Set(Property("pheomelanin")(pheomelanin->GetSDLValue()));
	props->Set(Property("eta")(eta->GetSDLValue()));
	props->Set(Property("beta_m")(betaM->GetSDLValue()));
	props->Set(Property("beta_n")(betaN->GetSDLValue()));
	props->Set(Property("alpha")(alpha->GetSDLValue()));

	return props;
}
