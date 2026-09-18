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

#include "luxrays/core/randomgen.h"

#include "slg/textures/whitenoise.h"

using namespace std;
using namespace luxrays;
using namespace slg;

//------------------------------------------------------------------------------
// White noise texture
//------------------------------------------------------------------------------

float WhiteNoiseTexture::GetFloatValue(const HitPoint &hitPoint) const {
	const Spectrum seed = GetTexture().GetSpectrumValue(hitPoint);

	TauswortheRandomGenerator rnd(
			SeedFromVector(seed.c[0], seed.c[1], seed.c[2]) + seedOffset);

	return rnd.floatValue();
}

Spectrum WhiteNoiseTexture::EvalSpectrumValue(const HitPoint &hitPoint) const {
	const Spectrum seed = GetTexture().GetSpectrumValue(hitPoint);

	// Three decorrelated draws give an uncorrelated RGB triplet.
	TauswortheRandomGenerator rnd(
			SeedFromVector(seed.c[0], seed.c[1], seed.c[2]) + seedOffset);

	return Spectrum(rnd.floatValue(), rnd.floatValue(), rnd.floatValue());
}

PropertiesUPtr WhiteNoiseTexture::ToProperties(const ImageMapCache &imgMapCache, const bool useRealFileName) const {
	auto props = std::make_unique<Properties>();

	const string name = GetName();
	props->Set(Property("scene.textures." + name + ".type")("whitenoise"));
	props->Set(Property("scene.textures." + name + ".texture")(GetTexture().GetSDLValue()));
	props->Set(Property("scene.textures." + name + ".seed")(seedOffset));

	return props;
}
// vim: autoindent noexpandtab tabstop=4 shiftwidth=4
