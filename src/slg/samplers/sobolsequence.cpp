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

#include <math.h>

#include "slg/samplers/sobolsequence.h"

using namespace std;
using namespace luxrays;
using namespace slg;

//------------------------------------------------------------------------------
// SobolSequence
//------------------------------------------------------------------------------

SobolSequence::SobolSequence() : directions(NULL) {
	rngPass = 0;
	rng0 = 0.f;
	rng1 = 0.f;
	blueNoiseEnable = false;
	blueNoiseSeed = 0;
}

SobolSequence::~SobolSequence() {
	delete[] directions;
}

void SobolSequence::RequestSamples(const u_int size) {
	directions = new u_int[size * SOBOL_BITS];
	GenerateDirectionVectors(directions, size);
}

u_int SobolSequence::SobolDimension(const u_int index, const u_int dimension) const {
	const u_int offset = dimension * SOBOL_BITS;
	u_int result = 0;
	u_int i = index;

	for (u_int j = 0; i; i >>= 1, j++) {
		if (i & 1)
			result ^= directions[offset + j];
	}

	return result;
}

u_int SobolSequence::BlueNoiseHash(u_int x) {
	// murmur3 32-bit finalizer (must match the GPU kernel version)
	x ^= x >> 16;
	x *= 0x7feb352du;
	x ^= x >> 15;
	x *= 0x846ca68bu;
	x ^= x >> 16;
	return x;
}

float SobolSequence::GetSample(const u_int pass, const u_int index) {
	u_int iResult;
	float shift;

	if (blueNoiseEnable) {
		// Blue-noise dithered sampling (Heitz et al. 2019): per-pixel
		// constant, per-dimension hashed digital shift + Cranley-Patterson
		// offset. The Sobol index is used unscrambled so a pixel walks its
		// own stratified prefix of the sequence across passes while
		// neighboring pixels are decorrelated by the per-pixel seed.
		const u_int dimSeed = BlueNoiseHash(blueNoiseSeed ^ (index * 0x9e3779b9u + 0x85ebca6bu));
		iResult = SobolDimension(pass, index) ^ dimSeed;
		shift = BlueNoiseHash(dimSeed ^ 0xc2b2ae35u) * (1.f / 4294967296.f);
	} else {
		// I scramble pass too in order avoid correlations visible with LIGHTCPU and BIDIRCPU
		iResult = SobolDimension(pass + rngPass, index);

		// Cranley-Patterson rotation to reduce visible regular patterns
		shift = (index & 1) ? rng0 : rng1;
	}

	const float fResult = iResult * (1.f / 0xffffffffu);
	const float val = fResult + shift;

	return val - floorf(val);
}
// vim: autoindent noexpandtab tabstop=4 shiftwidth=4
