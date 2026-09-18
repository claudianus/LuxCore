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

#ifndef _SLG_SOBOL_SEQUENCE_H
#define	_SLG_SOBOL_SEQUENCE_H

#include "luxrays/core/randomgen.h"

#include "slg/slg.h"
#include "slg/samplers/sampler.h"

namespace slg {

//------------------------------------------------------------------------------
// SobolSequence
//------------------------------------------------------------------------------

class SobolSequence {
public:
	SobolSequence();
	~SobolSequence();
	
	void RequestSamples(const u_int size);
	float GetSample(const u_int pass, const u_int index);

	u_int rngPass;
	float rng0, rng1;

	// Blue-noise dithered sampling (Heitz et al. 2019): when enabled, each
	// dimension of each pixel is randomized with a hashed digital shift plus
	// a hashed Cranley-Patterson offset. The seed is constant per pixel
	// (across passes) so a pixel keeps its own stratified Sobol prefix while
	// neighboring pixels use decorrelated dithers.
	void SetBlueNoiseSeed(const u_int seed) {
		blueNoiseSeed = seed;
		blueNoiseEnable = true;
	}
	void DisableBlueNoise() { blueNoiseEnable = false; }

	// murmur3 32-bit finalizer (must match the GPU kernel version)
	static u_int BlueNoiseHash(u_int x);

	static void GenerateDirectionVectors(u_int *vectors, const u_int dimensions);
private:
	u_int SobolDimension(const u_int index, const u_int dimension) const;

	u_int *directions;
	bool blueNoiseEnable;
	u_int blueNoiseSeed;
};

}

#endif	/* _SLG_SOBOL_SEQUENCE_H */
// vim: autoindent noexpandtab tabstop=4 shiftwidth=4
