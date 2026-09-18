#line 2 "texture_whitenoise_funcs.cl"

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
// White noise texture
//
// Mirrors slg::WhiteNoiseTexture: spatial hash of a 3D seed -> pseudo-random.
//------------------------------------------------------------------------------

// Mirrors WhiteNoiseTexture::SeedFromVector: bit-mix the component bit
// patterns then an integer avalanche. as_uint reinterprets bits, so it is
// deterministic for all float values and bit-identical to the CPU union
// reinterpretation (arithmetic float->uint casts are undefined for negatives).
OPENCL_FORCE_INLINE uint WhiteNoiseTexture_SeedFromVector(const float3 v) {
	uint h = as_uint(v.x) * 0x85ebca6bu ^ as_uint(v.y) * 0xc2b2ae35u ^
			as_uint(v.z) * 0x27d4eb2fu;
	h ^= h >> 16u; h *= 0x85ebca6bu;
	h ^= h >> 13u; h *= 0xc2b2ae35u;
	h ^= h >> 16u;
	return h;
}

OPENCL_FORCE_INLINE float WhiteNoiseTexture_ConstEvaluateFloat(const float3 v,
		const uint seedOffset) {
	Seed seed;
	Rnd_Init(WhiteNoiseTexture_SeedFromVector(v) + seedOffset, &seed);

	return Rnd_FloatValue(&seed);
}

OPENCL_FORCE_INLINE float3 WhiteNoiseTexture_ConstEvaluateSpectrum(const float3 v,
		const uint seedOffset) {
	Seed seed;
	Rnd_Init(WhiteNoiseTexture_SeedFromVector(v) + seedOffset, &seed);

	return MAKE_FLOAT3(Rnd_FloatValue(&seed), Rnd_FloatValue(&seed),
			Rnd_FloatValue(&seed));
}

OPENCL_FORCE_NOT_INLINE void WhiteNoiseTexture_EvalOp(
		__global const Texture* restrict texture,
		const TextureEvalOpType evalType,
		__global float *evalStack,
		uint *evalStackOffset,
		__global const HitPoint *hitPoint,
		const float sampleDistance
		TEXTURES_PARAM_DECL) {
	switch (evalType) {
		case EVAL_FLOAT: {
			float3 seed;
			EvalStack_PopFloat3(seed);

			const float eval = WhiteNoiseTexture_ConstEvaluateFloat(seed,
					texture->whiteNoiseTex.seedOffset);
			EvalStack_PushFloat(eval);
			break;
		}
		case EVAL_SPECTRUM: {
			float3 seed;
			EvalStack_PopFloat3(seed);

			const float3 eval = WhiteNoiseTexture_ConstEvaluateSpectrum(seed,
					texture->whiteNoiseTex.seedOffset);
			EvalStack_PushFloat3(eval);
			break;
		}
		case EVAL_BUMP_GENERIC_OFFSET_U:
			Texture_EvalOpGenericBumpOffsetU(evalStack, evalStackOffset,
					hitPoint, sampleDistance);
			break;
		case EVAL_BUMP_GENERIC_OFFSET_V:
			Texture_EvalOpGenericBumpOffsetV(evalStack, evalStackOffset,
					hitPoint, sampleDistance);
			break;
		case EVAL_BUMP:
			Texture_EvalOpGenericBump(evalStack, evalStackOffset,
					hitPoint, sampleDistance);
			break;
		default:
			// Something wrong here
			break;
	}
}
// vim: autoindent noexpandtab tabstop=4 shiftwidth=4
