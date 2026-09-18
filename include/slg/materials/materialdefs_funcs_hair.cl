#line 2 "materialdefs_funcs_hair.cl"

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

//------------------------------------------------------------------------------
// Hair material (Chiang et al. 2019, ported from pbrt-v3 HairBSDF)
// Mirror of src/slg/materials/hairmat.cpp for the device path.
//------------------------------------------------------------------------------

#define HAIR_PMAX 3
#define HAIR_SQRTPIOVER8 0.626657069f
// Vertex AOV layers carrying the object-space strand tangent (must match
// HAIR_TANGENT_*_DATA_INDEX in slg/shapes/strands.h)
#define HAIR_TANGENT_X_DATA_INDEX 4
#define HAIR_TANGENT_Y_DATA_INDEX 5
#define HAIR_TANGENT_Z_DATA_INDEX 6

OPENCL_FORCE_INLINE float Hair_SafeSqrt(const float x) { return sqrt(fmax(0.f, x)); }
OPENCL_FORCE_INLINE float Hair_SafeASin(const float x) { return asin(clamp(x, -1.f, 1.f)); }
OPENCL_FORCE_INLINE float Hair_Fract(const float x) { return x - floor(x); }

OPENCL_FORCE_INLINE float Hair_I0(const float x) {
	float val = 0.f;
	float x2i = 1.f;
	int ifact = 1;
	int i4 = 1;
	for (int i = 0; i < 10; ++i) {
		if (i > 1) ifact *= i;
		val += x2i / (i4 * ifact * ifact);
		x2i *= x * x;
		i4 *= 4;
	}
	return val;
}

OPENCL_FORCE_INLINE float Hair_LogI0(const float x) {
	if (x > 12.f)
		return x + 0.5f * (-log(2.f * M_PI_F) + log(1.f / x) + 1.f / (8.f * x));
	else
		return log(Hair_I0(x));
}

OPENCL_FORCE_INLINE float Hair_Mp(const float cosThetaI, const float cosThetaO,
		const float sinThetaI, const float sinThetaO, const float v) {
	const float a = cosThetaI * cosThetaO / v;
	const float b = sinThetaI * sinThetaO / v;
	return (v <= .1f)
		? exp(Hair_LogI0(a) - b - 1.f / v + 0.6931f + log(1.f / (2.f * v)))
		: (exp(-b) * Hair_I0(a)) / (sinh(1.f / v) * 2.f * v);
}

OPENCL_FORCE_INLINE float Hair_Phi(const int p, const float gammaO, const float gammaT) {
	return 2.f * p * gammaT - 2.f * gammaO + p * M_PI_F;
}

OPENCL_FORCE_INLINE float Hair_Logistic(float x, const float s) {
	x = fabs(x);
	return exp(-x / s) / (s * (1.f + exp(-x / s)) * (1.f + exp(-x / s)));
}

OPENCL_FORCE_INLINE float Hair_LogisticCDF(const float x, const float s) {
	return 1.f / (1.f + exp(-x / s));
}

OPENCL_FORCE_INLINE float Hair_TrimmedLogistic(const float x, const float s,
		const float a, const float b) {
	return Hair_Logistic(x, s) / (Hair_LogisticCDF(b, s) - Hair_LogisticCDF(a, s));
}

OPENCL_FORCE_INLINE float Hair_Np(float phi, const int p, const float s,
		const float gammaO, const float gammaT) {
	float dphi = phi - Hair_Phi(p, gammaO, gammaT);
	// Remap dphi to [-pi, pi]
	while (dphi > M_PI_F) dphi -= 2.f * M_PI_F;
	while (dphi < -M_PI_F) dphi += 2.f * M_PI_F;
	return Hair_TrimmedLogistic(dphi, s, -M_PI_F, M_PI_F);
}

OPENCL_FORCE_INLINE float Hair_SampleTrimmedLogistic(const float u, const float s,
		const float a, const float b) {
	const float k = Hair_LogisticCDF(b, s) - Hair_LogisticCDF(a, s);
	const float x = -s * log(1.f / (u * k + Hair_LogisticCDF(a, s)) - 1.f);
	return clamp(x, a, b);
}

OPENCL_FORCE_INLINE void Hair_Ap(const float cosThetaO, const float eta, const float h,
		const float3 T, __private float3 *ap) {
	const float cosGammaO = Hair_SafeSqrt(1.f - h * h);
	const float cosTheta = cosThetaO * cosGammaO;
	const float f = FresnelCauchy_Evaluate(eta, cosTheta);
	ap[0] = f;

	ap[1] = (1.f - f) * (1.f - f) * T;

	for (int p = 2; p < HAIR_PMAX; ++p)
		ap[p] = ap[p - 1] * T * f;

	ap[HAIR_PMAX] = ap[HAIR_PMAX - 1] * f * T /
			(MAKE_FLOAT3(1.f, 1.f, 1.f) - T * f);
}

OPENCL_FORCE_INLINE void Hair_TiltThetaO(const int p,
		__private const float *sin2k, __private const float *cos2k,
		const float sinThetaO, const float cosThetaO,
		__private float *sinThetaOp, __private float *cosThetaOp) {
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
	*cosThetaOp = fabs(*cosThetaOp);
}

OPENCL_FORCE_INLINE void Hair_DemuxFloat(const float u, __private float *a, __private float *b) {
	*a = Hair_Fract(u * 65536.f);
	*b = Hair_Fract(u * 4294967296.f);
}

// Per-hit hair state
typedef struct {
	float3 tangent;   // strand tangent in the local shading frame
	float gammaO;
	float eta;
	float3 sigma_a;
	float v[4];
	float s;
	float sin2kAlpha[3];
	float cos2kAlpha[3];
} HairContext;

// Hair frame basis: e1 = radial (shading normal projected perp. to the
// tangent), e2 = cross(tangent, e1). Local shading normal is +Z.
OPENCL_FORCE_INLINE void Hair_Frame(const float3 T, __private float3 *e1, __private float3 *e2) {
	float3 nPerp = MAKE_FLOAT3(0.f, 0.f, 1.f) - T * T.z;
	if (dot(nPerp, nPerp) < 1e-12f)
		nPerp = MAKE_FLOAT3(1.f, 0.f, 0.f) - T * T.x;
	*e1 = normalize(nPerp);
	*e2 = cross(T, *e1);
}

OPENCL_FORCE_INLINE float3 Hair_EvalSigmaA(__global const Material* restrict material,
		__global const HitPoint *hitPoint, const float betaN
		MATERIALS_PARAM_DECL) {
	if (material->hair.sigmaATexIndex != NULL_INDEX)
		return fmax(Texture_GetSpectrumValue(material->hair.sigmaATexIndex,
				hitPoint TEXTURES_PARAM), BLACK);

	if (material->hair.colorTexIndex != NULL_INDEX) {
		const float3 c = Spectrum_Clamp(Texture_GetSpectrumValue(material->hair.colorTexIndex,
				hitPoint TEXTURES_PARAM));
		// SigmaAFromReflectance (Chiang 2019)
		const float denom = 5.969f - 0.215f * betaN + 2.532f * betaN * betaN -
				10.73f * betaN * betaN * betaN + 5.574f * pow(betaN, 4.f) +
				0.245f * pow(betaN, 5.f);
		float3 sa;
		sa.x = log(fmax(c.x, 1e-5f)) / denom; sa.x *= sa.x;
		sa.y = log(fmax(c.y, 1e-5f)) / denom; sa.y *= sa.y;
		sa.z = log(fmax(c.z, 1e-5f)) / denom; sa.z *= sa.z;
		return sa;
	}

	// Melanin concentration model (default: brown hair 1.3 / 0)
	const float ce = (material->hair.eumelaninTexIndex != NULL_INDEX) ?
			fmax(0.f, Texture_GetFloatValue(material->hair.eumelaninTexIndex,
			hitPoint TEXTURES_PARAM)) : 1.3f;
	const float cp = (material->hair.pheomelaninTexIndex != NULL_INDEX) ?
			fmax(0.f, Texture_GetFloatValue(material->hair.pheomelaninTexIndex,
			hitPoint TEXTURES_PARAM)) : 0.f;
	return MAKE_FLOAT3(
		ce * 0.419f + cp * 0.187f,
		ce * 0.697f + cp * 0.4f,
		ce * 1.37f + cp * 1.05f);
}

OPENCL_FORCE_INLINE void Hair_SetupContext(__global const Material* restrict material,
		__global const HitPoint *hitPoint, const float3 localEyeDir,
		__private HairContext *ctx
		MATERIALS_PARAM_DECL) {
	// Strand tangent: vertex AOV (object space) -> world -> local frame
	float3 tObj = MAKE_FLOAT3(
			HitPoint_GetVertexAOV(hitPoint, HAIR_TANGENT_X_DATA_INDEX EXTMESH_PARAM),
			HitPoint_GetVertexAOV(hitPoint, HAIR_TANGENT_Y_DATA_INDEX EXTMESH_PARAM),
			HitPoint_GetVertexAOV(hitPoint, HAIR_TANGENT_Z_DATA_INDEX EXTMESH_PARAM));
	if (dot(tObj, tObj) < 1e-12f) {
		// No strand data: fall back to dpdv (fiber direction for
		// strand-like UV layouts)
		tObj = VLOAD3F(&hitPoint->dpdv.x);
	}
	const float3 tWorld = normalize(Transform_ApplyVector(&hitPoint->localToWorld, tObj));
	Frame frame;
	Frame_Set_Private(&frame, VLOAD3F(&hitPoint->dpdu.x),
			VLOAD3F(&hitPoint->dpdv.x), VLOAD3F(&hitPoint->shadeN.x));
	ctx->tangent = normalize(Frame_ToLocal_Private(&frame, tWorld));

	const float bm = clamp(Texture_GetFloatValue(material->hair.betaMTexIndex,
			hitPoint TEXTURES_PARAM), 1e-2f, 1.f);
	const float bn = clamp(Texture_GetFloatValue(material->hair.betaNTexIndex,
			hitPoint TEXTURES_PARAM), 1e-2f, 1.f);
	const float a = Texture_GetFloatValue(material->hair.alphaTexIndex,
			hitPoint TEXTURES_PARAM);
	ctx->eta = Texture_GetFloatValue(material->hair.etaTexIndex,
			hitPoint TEXTURES_PARAM);

	// Longitudinal variance from beta_m
	ctx->v[0] = (0.726f * bm + 0.812f * bm * bm + 3.7f * pow(bm, 20.f));
	ctx->v[0] *= ctx->v[0];
	ctx->v[1] = .25f * ctx->v[0];
	ctx->v[2] = 4.f * ctx->v[0];
	ctx->v[3] = ctx->v[2];

	// Azimuthal logistic scale from beta_n
	ctx->s = HAIR_SQRTPIOVER8 *
			(0.265f * bn + 1.194f * bn * bn + 5.372f * pow(bn, 22.f));

	// Scale tilt (alpha) terms
	ctx->sin2kAlpha[0] = sin(a * M_PI_F / 180.f);
	ctx->cos2kAlpha[0] = Hair_SafeSqrt(1.f - ctx->sin2kAlpha[0] * ctx->sin2kAlpha[0]);
	for (int i = 1; i < 3; ++i) {
		ctx->sin2kAlpha[i] = 2.f * ctx->cos2kAlpha[i - 1] * ctx->sin2kAlpha[i - 1];
		ctx->cos2kAlpha[i] = ctx->cos2kAlpha[i - 1] * ctx->cos2kAlpha[i - 1] -
				ctx->sin2kAlpha[i - 1] * ctx->sin2kAlpha[i - 1];
	}

	// Azimuthal offset h from the hit radial direction vs. the outgoing
	// ray azimuth plane
	float3 e1, e2;
	Hair_Frame(ctx->tangent, &e1, &e2);
	float3 e1wo = localEyeDir - ctx->tangent * dot(localEyeDir, ctx->tangent);
	if (dot(e1wo, e1wo) < 1e-12f)
		e1wo = e1; // outgoing ray parallel to the fiber: h = 0
	else
		e1wo = normalize(e1wo);
	const float3 e2wo = cross(ctx->tangent, e1wo);
	const float h = clamp(dot(e1, e2wo), -1.f, 1.f);
	ctx->gammaO = Hair_SafeASin(h);

	ctx->sigma_a = Hair_EvalSigmaA(material, hitPoint, bn MATERIALS_PARAM);
}

OPENCL_FORCE_INLINE float3 Hair_f(__global const HitPoint *hitPoint,
		__private const HairContext *ctx, const float3 localLightDir, const float3 localEyeDir) {
	const float3 T = ctx->tangent;
	float3 e1, e2;
	Hair_Frame(T, &e1, &e2);

	const float sinThetaO = dot(localEyeDir, T);
	const float cosThetaO = Hair_SafeSqrt(1.f - sinThetaO * sinThetaO);
	const float phiO = atan2(dot(localEyeDir, e2), dot(localEyeDir, e1));

	const float sinThetaI = dot(localLightDir, T);
	const float cosThetaI = Hair_SafeSqrt(1.f - sinThetaI * sinThetaI);
	const float phiI = atan2(dot(localLightDir, e2), dot(localLightDir, e1));

	const float sinThetaT = sinThetaO / ctx->eta;
	const float cosThetaT = Hair_SafeSqrt(1.f - sinThetaT * sinThetaT);

	const float h = sin(ctx->gammaO);
	const float etap = sqrt(ctx->eta * ctx->eta - sinThetaO * sinThetaO) / cosThetaO;
	const float sinGammaT = h / etap;
	const float cosGammaT = Hair_SafeSqrt(1.f - sinGammaT * sinGammaT);
	const float gammaT = Hair_SafeASin(sinGammaT);

	const float3 Tt = Spectrum_Exp(ctx->sigma_a * (-2.f * cosGammaT / cosThetaT));

	const float phi = phiI - phiO;
	float3 ap[HAIR_PMAX + 1];
	Hair_Ap(cosThetaO, ctx->eta, h, Tt, ap);

	float3 fsum = BLACK;
	for (int p = 0; p < HAIR_PMAX; ++p) {
		float sinThetaOp, cosThetaOp;
		Hair_TiltThetaO(p, ctx->sin2kAlpha, ctx->cos2kAlpha,
				sinThetaO, cosThetaO, &sinThetaOp, &cosThetaOp);
		fsum += Hair_Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, ctx->v[p]) *
				ap[p] * Hair_Np(phi, p, ctx->s, ctx->gammaO, gammaT);
	}
	fsum += Hair_Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, ctx->v[HAIR_PMAX]) *
			ap[HAIR_PMAX] * (1.f / (2.f * M_PI_F));
	// LuxCore wants f * cos: pbrt's BSDF value is fsum/|cosThetaI| so the
	// cosine-multiplied value is exactly fsum.
	return fsum;
}

OPENCL_FORCE_INLINE void Hair_ComputeApPdf(__private const HairContext *ctx,
		const float cosThetaO, __private float *apPdf) {
	const float sinThetaO = Hair_SafeSqrt(1.f - cosThetaO * cosThetaO);
	const float sinThetaT = sinThetaO / ctx->eta;
	const float cosThetaT = Hair_SafeSqrt(1.f - sinThetaT * sinThetaT);

	const float h = sin(ctx->gammaO);
	const float etap = sqrt(ctx->eta * ctx->eta - sinThetaO * sinThetaO) / cosThetaO;
	const float sinGammaT = h / etap;
	const float cosGammaT = Hair_SafeSqrt(1.f - sinGammaT * sinGammaT);

	const float3 T = Spectrum_Exp(ctx->sigma_a * (-2.f * cosGammaT / cosThetaT));
	float3 ap[HAIR_PMAX + 1];
	Hair_Ap(cosThetaO, ctx->eta, h, T, ap);

	float sumY = 0.f;
	for (int i = 0; i <= HAIR_PMAX; ++i)
		sumY += Spectrum_Y(ap[i]);
	for (int i = 0; i <= HAIR_PMAX; ++i)
		apPdf[i] = (sumY > 0.f) ? (Spectrum_Y(ap[i]) / sumY) : 0.f;
}

OPENCL_FORCE_INLINE float Hair_EvalPdf(__global const HitPoint *hitPoint,
		__private const HairContext *ctx, const float3 localLightDir, const float3 localEyeDir) {
	const float3 T = ctx->tangent;
	float3 e1, e2;
	Hair_Frame(T, &e1, &e2);

	const float sinThetaO = dot(localEyeDir, T);
	const float cosThetaO = Hair_SafeSqrt(1.f - sinThetaO * sinThetaO);
	const float sinThetaI = dot(localLightDir, T);
	const float cosThetaI = Hair_SafeSqrt(1.f - sinThetaI * sinThetaI);
	const float phiO = atan2(dot(localEyeDir, e2), dot(localEyeDir, e1));
	const float phiI = atan2(dot(localLightDir, e2), dot(localLightDir, e1));
	const float phi = phiI - phiO;

	const float h = sin(ctx->gammaO);
	const float etap = sqrt(ctx->eta * ctx->eta - sinThetaO * sinThetaO) / cosThetaO;
	const float sinGammaT = h / etap;
	const float gammaT = Hair_SafeASin(sinGammaT);

	float apPdf[HAIR_PMAX + 1];
	Hair_ComputeApPdf(ctx, cosThetaO, apPdf);

	float pdf = 0.f;
	for (int p = 0; p < HAIR_PMAX; ++p) {
		float sinThetaOp, cosThetaOp;
		Hair_TiltThetaO(p, ctx->sin2kAlpha, ctx->cos2kAlpha,
				sinThetaO, cosThetaO, &sinThetaOp, &cosThetaOp);
		pdf += Hair_Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, ctx->v[p]) *
				apPdf[p] * Hair_Np(phi, p, ctx->s, ctx->gammaO, gammaT);
	}
	pdf += Hair_Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, ctx->v[HAIR_PMAX]) *
			apPdf[HAIR_PMAX] * (1.f / (2.f * M_PI_F));
	return pdf;
}

OPENCL_FORCE_INLINE void HairMaterial_Albedo(__global const Material* restrict material,
		__global const HitPoint *hitPoint,
		__global float *evalStack, uint *evalStackOffset
		MATERIALS_PARAM_DECL) {
	HairContext ctx;
	Hair_SetupContext(material, hitPoint, MAKE_FLOAT3(0.f, 0.f, 1.f), &ctx MATERIALS_PARAM);
	float apPdf[HAIR_PMAX + 1];
	Hair_ComputeApPdf(&ctx, 1.f, apPdf);
	float sum = 0.f;
	for (int i = 0; i <= HAIR_PMAX; ++i)
		sum += apPdf[i];
	EvalStack_PushFloat3(MAKE_FLOAT3(clamp(sum * 0.5f, 0.f, 1.f),
			clamp(sum * 0.5f, 0.f, 1.f), clamp(sum * 0.5f, 0.f, 1.f)));
}

OPENCL_FORCE_INLINE void HairMaterial_GetInteriorVolume(__global const Material* restrict material,
		__global const HitPoint *hitPoint,
		__global float *evalStack, uint *evalStackOffset
		MATERIALS_PARAM_DECL) {
	DefaultMaterial_GetInteriorVolume(material, hitPoint, evalStack, evalStackOffset MATERIALS_PARAM);
}

OPENCL_FORCE_INLINE void HairMaterial_GetExteriorVolume(__global const Material* restrict material,
		__global const HitPoint *hitPoint,
		__global float *evalStack, uint *evalStackOffset
		MATERIALS_PARAM_DECL) {
	DefaultMaterial_GetExteriorVolume(material, hitPoint, evalStack, evalStackOffset MATERIALS_PARAM);
}

OPENCL_FORCE_INLINE void HairMaterial_GetPassThroughTransparency(__global const Material* restrict material,
		__global const HitPoint *hitPoint,
		__global float *evalStack, uint *evalStackOffset
		MATERIALS_PARAM_DECL) {
	DefaultMaterial_GetPassThroughTransparency(material, hitPoint, evalStack, evalStackOffset MATERIALS_PARAM);
}

OPENCL_FORCE_INLINE void HairMaterial_GetEmittedRadiance(__global const Material* restrict material,
		__global const HitPoint *hitPoint,
		__global float *evalStack, uint *evalStackOffset
		MATERIALS_PARAM_DECL) {
	DefaultMaterial_GetEmittedRadiance(material, hitPoint, evalStack, evalStackOffset MATERIALS_PARAM);
}

OPENCL_FORCE_INLINE void HairMaterial_Evaluate(__global const Material* restrict material,
		__global const HitPoint *hitPoint,
		__global float *evalStack, uint *evalStackOffset
		MATERIALS_PARAM_DECL) {
	float3 lightDir, eyeDir;
	EvalStack_PopFloat3(eyeDir);
	EvalStack_PopFloat3(lightDir);

	HairContext ctx;
	Hair_SetupContext(material, hitPoint, eyeDir, &ctx MATERIALS_PARAM);

	const float3 result = Hair_f(hitPoint, &ctx, lightDir, eyeDir);
	if (Spectrum_IsBlack(result)) {
		MATERIAL_EVALUATE_RETURN_BLACK;
	}

	const BSDFEvent event = GLOSSY | REFLECT | TRANSMIT;
	const float directPdfW = Hair_EvalPdf(hitPoint, &ctx, lightDir, eyeDir);

	EvalStack_PushFloat3(result);
	EvalStack_PushBSDFEvent(event);
	EvalStack_PushFloat(directPdfW);
}

OPENCL_FORCE_INLINE void HairMaterial_Sample(__global const Material* restrict material,
		__global const HitPoint *hitPoint,
		__global float *evalStack, uint *evalStackOffset
		MATERIALS_PARAM_DECL) {
	float u0, u1, passThroughEvent;
	EvalStack_PopFloat(passThroughEvent);
	EvalStack_PopFloat(u1);
	EvalStack_PopFloat(u0);
	float3 fixedDir;
	EvalStack_PopFloat3(fixedDir);

	HairContext ctx;
	Hair_SetupContext(material, hitPoint, fixedDir, &ctx MATERIALS_PARAM);

	const float3 T = ctx.tangent;
	float3 e1, e2;
	Hair_Frame(T, &e1, &e2);

	const float sinThetaO = dot(fixedDir, T);
	const float cosThetaO = Hair_SafeSqrt(1.f - sinThetaO * sinThetaO);
	const float phiO = atan2(dot(fixedDir, e2), dot(fixedDir, e1));

	// Four sample dimensions from the three supplied values
	float u[2][2];
	Hair_DemuxFloat(u0, &u[0][0], &u[0][1]);
	Hair_DemuxFloat(u1, &u[1][0], &u[1][1]);

	// Choose the lobe p from the Ap pdf
	float apPdf[HAIR_PMAX + 1];
	Hair_ComputeApPdf(&ctx, cosThetaO, apPdf);
	int p;
	float sel = u[0][0];
	for (p = 0; p < HAIR_PMAX; ++p) {
		if (sel < apPdf[p]) break;
		sel -= apPdf[p];
	}

	float sinThetaOp, cosThetaOp;
	Hair_TiltThetaO(p, ctx.sin2kAlpha, ctx.cos2kAlpha,
			sinThetaO, cosThetaO, &sinThetaOp, &cosThetaOp);

	// Sample M_p -> thetaI
	const float u10 = fmax(u[1][0], 1e-5f);
	const float cosTheta =
			1.f + ctx.v[p] * log(u10 + (1.f - u10) * exp(-2.f / ctx.v[p]));
	const float sinTheta = Hair_SafeSqrt(1.f - cosTheta * cosTheta);
	const float cosPhi = cos(2.f * M_PI_F * u[1][1]);
	const float sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhi * cosThetaOp;
	const float cosThetaI = Hair_SafeSqrt(1.f - sinThetaI * sinThetaI);

	// Sample N_p -> dphi
	const float h = sin(ctx.gammaO);
	const float etap = sqrt(ctx.eta * ctx.eta - sinThetaO * sinThetaO) / cosThetaO;
	const float sinGammaT = h / etap;
	const float gammaT = Hair_SafeASin(sinGammaT);
	float dphi;
	if (p < HAIR_PMAX)
		dphi = Hair_Phi(p, ctx.gammaO, gammaT) +
				Hair_SampleTrimmedLogistic(u[0][1], ctx.s, -M_PI_F, M_PI_F);
	else
		dphi = 2.f * M_PI_F * u[0][1];

	const float phiI = phiO + dphi;
	const float3 sampledDir = sinThetaI * T +
			cosThetaI * cos(phiI) * e1 + cosThetaI * sin(phiI) * e2;

	// Pdf for the sampled direction
	float pdf = 0.f;
	for (int pp = 0; pp < HAIR_PMAX; ++pp) {
		float stp, ctp;
		Hair_TiltThetaO(pp, ctx.sin2kAlpha, ctx.cos2kAlpha,
				sinThetaO, cosThetaO, &stp, &ctp);
		pdf += Hair_Mp(cosThetaI, ctp, sinThetaI, stp, ctx.v[pp]) *
				apPdf[pp] * Hair_Np(dphi, pp, ctx.s, ctx.gammaO, gammaT);
	}
	pdf += Hair_Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, ctx.v[HAIR_PMAX]) *
			apPdf[HAIR_PMAX] * (1.f / (2.f * M_PI_F));

	if (pdf <= 0.f) {
		MATERIAL_SAMPLE_RETURN_BLACK;
	}

	const BSDFEvent event = (p == 0) ? (GLOSSY | REFLECT) : (GLOSSY | TRANSMIT);

	const float3 result = Hair_f(hitPoint, &ctx, sampledDir, fixedDir) / pdf;
	if (Spectrum_IsBlack(result)) {
		MATERIAL_SAMPLE_RETURN_BLACK;
	}

	EvalStack_PushFloat3(result);
	EvalStack_PushFloat3(sampledDir);
	EvalStack_PushFloat(pdf);
	EvalStack_PushBSDFEvent(event);
}

//------------------------------------------------------------------------------
// Material specific EvalOp
//------------------------------------------------------------------------------

OPENCL_FORCE_NOT_INLINE void HairMaterial_EvalOp(
		__global const Material* restrict material,
		const MaterialEvalOpType evalType,
		__global float *evalStack,
		uint *evalStackOffset,
		__global const HitPoint *hitPoint
		MATERIALS_PARAM_DECL) {
	switch (evalType) {
		case EVAL_ALBEDO:
			HairMaterial_Albedo(material, hitPoint, evalStack, evalStackOffset MATERIALS_PARAM);
			break;
		case EVAL_GET_INTERIOR_VOLUME:
			HairMaterial_GetInteriorVolume(material, hitPoint, evalStack, evalStackOffset MATERIALS_PARAM);
			break;
		case EVAL_GET_EXTERIOR_VOLUME:
			HairMaterial_GetExteriorVolume(material, hitPoint, evalStack, evalStackOffset MATERIALS_PARAM);
			break;
		case EVAL_GET_EMITTED_RADIANCE:
			HairMaterial_GetEmittedRadiance(material, hitPoint, evalStack, evalStackOffset MATERIALS_PARAM);
			break;
		case EVAL_GET_PASS_TROUGH_TRANSPARENCY:
			HairMaterial_GetPassThroughTransparency(material, hitPoint, evalStack, evalStackOffset MATERIALS_PARAM);
			break;
		case EVAL_EVALUATE:
			HairMaterial_Evaluate(material, hitPoint, evalStack, evalStackOffset MATERIALS_PARAM);
			break;
		case EVAL_SAMPLE:
			HairMaterial_Sample(material, hitPoint, evalStack, evalStackOffset MATERIALS_PARAM);
			break;
		default:
			// Something wrong here
			break;
	}
}
// vim: autoindent noexpandtab tabstop=4 shiftwidth=4
