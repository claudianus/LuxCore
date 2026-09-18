/***************************************************************************
 * Copyright 1998-2020 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxCoreRender.                                   *
 *                                                                         *
 *   Licensed under the Apache License, Version 2.0 (the "License");       *
 *   you may not use this file except in compliance with the License.      *
 *   You may obtain a copy of the License at                               *
 *                                                                         *
 *   Unless required by applicable law or agreed to in writing, software   *
 *   distributed under the License is distributed on an "AS IS" BASIS,     *
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or       *
 *   implied.                                                              *
 *                                                                         *
 *   See the License for the specific language governing permissions and   *
 *   limitations under the License.                                        *
 ***************************************************************************/

#include <cmath>
#include <algorithm>

#include "luxrays/core/color/spectral.h"
#include "luxrays/core/color/spectrumwavelengths.h"
#include "luxrays/core/color/spds/data/rgbE_32.h"
#include "luxrays/core/color/spds/data/rgbD65_32.h"

using namespace luxrays;

namespace {
	bool g_enabled = false;
	thread_local PathWavelengths g_sw;
	thread_local bool g_hasSW = false;

	// Linear-interpolated sample of a raw SPD table, same math as SPD::Sample
	inline float SampleTable(const float *data, const u_int n,
			const float start, const float end, const float lambda) {
		if (n <= 1 || lambda < start || lambda > end)
			return 0.f;
		const float x = (lambda - start) * ((n - 1) / (end - start));
		const u_int b0 = Floor2UInt(std::max(x, 0.f));
		const u_int b1 = std::min(b0 + 1, n - 1);
		return std::lerp(data[b0], data[b1], x - b0);
	}

	struct BasisSet {
		const float *white, *cyan, *magenta, *yellow, *red, *green, *blue;
		u_int n; float start, end;
	};

	const BasisSet reflBasis = {
		refrgb2spect_white, refrgb2spect_cyan, refrgb2spect_magenta,
		refrgb2spect_yellow, refrgb2spect_red, refrgb2spect_green, refrgb2spect_blue,
		refrgb2spect_bins, refrgb2spect_start, refrgb2spect_end
	};
	const BasisSet illumBasis = {
		illumrgb2spect_white, illumrgb2spect_cyan, illumrgb2spect_magenta,
		illumrgb2spect_yellow, illumrgb2spect_red, illumrgb2spect_green,
		illumrgb2spect_blue,
		illumrgb2spect_bins, illumrgb2spect_start, illumrgb2spect_end
	};

	// Smits RGB->SPD decomposition evaluated directly at the path wavelengths:
	// rgb is decomposed into (white, secondary, primary) basis weights, then
	// each basis is sampled at every bin. The result is then renormalized so
	// its CIE luminance over the sampled bins matches rgb.Y() (metameric
	// round-trip: flat/achromatic inputs reproduce themselves exactly under
	// the film projector). Negative bins are clamped so sample validity
	// checks don't reject legitimately-upsampled values.
	Spectrum Upsample(const Spectrum &rgb, const PathWavelengths &sw,
			const BasisSet &basis, const float yTarget) {
		const float r = rgb.c[0], g = rgb.c[1], b = rgb.c[2];

		// Achromatic input: a flat SPD is the only metameric-exact choice —
		// the Smits white basis is NOT perfectly flat and would inject a
		// chromatic bias (observed: cyan cast on neutral walls). Keep the
		// round-trip exact: flat RGB -> flat bins -> flat projected RGB.
		if (r == g && g == b) {
			Spectrum out;
			for (u_int i = 0; i < SPECTRAL_BINS; ++i)
				out.c[i] = (sw.aliveMask & (1U << i)) ? std::max(r, 0.f) : 0.f;
			return out;
		}

		float wW; const float *sec; float wSec; const float *prim; float wPrim;
		if (r <= g && r <= b) {
			wW = r;
			sec = basis.cyan;
			if (g <= b) { wSec = g - r; prim = basis.blue;  wPrim = b - g; }
			else        { wSec = b - r; prim = basis.green; wPrim = g - b; }
		} else if (g <= r && g <= b) {
			wW = g;
			sec = basis.magenta;
			if (r <= b) { wSec = r - g; prim = basis.blue; wPrim = b - r; }
			else        { wSec = b - g; prim = basis.red;  wPrim = r - b; }
		} else {
			wW = b;
			sec = basis.yellow;
			if (r <= g) { wSec = r - b; prim = basis.green; wPrim = g - r; }
			else        { wSec = g - b; prim = basis.red;   wPrim = r - g; }
		}

		Spectrum out;
		float lum = 0.f, nY = 0.f;
		for (u_int i = 0; i < SPECTRAL_BINS; ++i) {
			const float lambda = sw.w[i];
			const float cy = SpectrumWavelengths::spd_ciey.Sample(lambda);
			// The basis is evaluated and the luminance matched over the
			// FULL drawn wavelength set — the material's spectral shape is
			// a property of the material, not of the path state. After a
			// dispersive collapse only the write is masked.
			float v = wW * SampleTable(basis.white, basis.n, basis.start, basis.end, lambda)
				+ wSec * SampleTable(sec, basis.n, basis.start, basis.end, lambda)
				+ wPrim * SampleTable(prim, basis.n, basis.start, basis.end, lambda);
			v = std::max(v, 0.f);
			lum += v * cy;
			nY += cy;
			out.c[i] = (sw.aliveMask & (1U << i)) ? v : 0.f;
		}

		// Normalize so the projected luminance (sum bins*ciey / sum ciey)
		// matches the input luminance. Note the CIE SPDs are pre-scaled by
		// 683*range/samples, so the raw sum must be divided by nY to keep
		// the bins in the same magnitude range as the RGB input.
		if (lum > 0.f && yTarget > 0.f && nY > 0.f) {
			const float s = (yTarget * nY) / lum;
			for (u_int i = 0; i < SPECTRAL_BINS; ++i)
				out.c[i] *= s;
		}

		return out;
	}
}

void Spectral::SetEnabled(const bool enabled) { g_enabled = enabled; }
bool Spectral::IsEnabled() { return g_enabled; }

void Spectral::SetPathWavelengths(const PathWavelengths &sw) {
	g_sw = sw;
	g_hasSW = true;
}
void Spectral::ClearPathWavelengths() { g_hasSW = false; }
const PathWavelengths *Spectral::Current() {
	return (g_enabled && g_hasSW) ? &g_sw : nullptr;
}

float Spectral::CollapseToHero() {
	if (!(g_enabled && g_hasSW))
		return 1.f;
	const u_int heroMask = 1u << g_sw.hero;
	if (g_sw.aliveMask == heroMask)
		return 1.f; // already collapsed
	g_sw.aliveMask = heroMask;
	// Uniform pick of the hero bin out of SPECTRAL_BINS: the surviving
	// estimate carries the whole path's spectral weight.
	return (float)SPECTRAL_BINS;
}

Spectrum Spectral::Reflectance(const Spectrum &rgb) {
	const PathWavelengths *sw = Current();
	return sw ? Reflectance(rgb, *sw) : rgb;
}
Spectrum Spectral::Emission(const Spectrum &rgb) {
	const PathWavelengths *sw = Current();
	return sw ? Emission(rgb, *sw) : rgb;
}
Spectrum Spectral::Reflectance(const Spectrum &rgb, const PathWavelengths &sw) {
	return Upsample(rgb, sw, reflBasis, rgb.Y());
}
Spectrum Spectral::Emission(const Spectrum &rgb, const PathWavelengths &sw) {
	return Upsample(rgb, sw, illumBasis, rgb.Y());
}

Spectrum Spectral::EvaluateSPD(const SPD &spd) {
	const PathWavelengths *sw = Current();
	return sw ? EvaluateSPD(spd, *sw) : Spectrum();
}
Spectrum Spectral::EvaluateSPD(const SPD &spd, const PathWavelengths &sw) {
	Spectrum out(0.f);
	for (u_int i = 0; i < SPECTRAL_BINS; ++i)
		if (sw.aliveMask & (1U << i))
			out.c[i] = spd.Sample(sw.w[i]);
	return out;
}

Spectrum Spectral::WithLuminance(const Spectrum &bins, const PathWavelengths &sw,
		const float yTarget) {
	float lum = 0.f, nY = 0.f;
	for (u_int i = 0; i < SPECTRAL_BINS; ++i) {
		const float cy = SpectrumWavelengths::spd_ciey.Sample(sw.w[i]);
		nY += cy; // full drawn set (dead bins are zeros, not missing samples)
		if (sw.aliveMask & (1U << i))
			lum += bins.c[i] * cy;
	}
	if (lum <= 0.f || yTarget <= 0.f || nY <= 0.f)
		return bins;
	return bins * (yTarget * nY / lum);
}

Spectrum Spectral::ProjectToRGB(const Spectrum &bins, const PathWavelengths &sw) {
	// Accumulate XYZ under the CIE matching functions at the live bins;
	// the white-point normalizer always spans the full drawn wavelength
	// set so a post-collapse single bin keeps its chromaticity.
	float X = 0.f, Y = 0.f, Z = 0.f;
	float nX = 0.f, nY = 0.f, nZ = 0.f;
	for (u_int i = 0; i < SPECTRAL_BINS; ++i) {
		const float lambda = sw.w[i];
		const float cx = SpectrumWavelengths::spd_ciex.Sample(lambda);
		const float cy = SpectrumWavelengths::spd_ciey.Sample(lambda);
		const float cz = SpectrumWavelengths::spd_ciez.Sample(lambda);
		nX += cx;
		nY += cy;
		nZ += cz;
		if (!(sw.aliveMask & (1U << i)))
			continue;
		X += bins.c[i] * cx;
		Y += bins.c[i] * cy;
		Z += bins.c[i] * cz;
	}
	if (nY <= 0.f)
		return Spectrum(0.f);

	// Normalize by the sampled white point so a flat spectrum reproduces its
	// value exactly (achromatic-invariant projection)
	const XYZColor whiteXYZ(nX, nY, nZ);
	const RGBColor whiteRGB = ColorSystem::DefaultColorSystem.ToRGB(whiteXYZ);
	const RGBColor rgb = ColorSystem::DefaultColorSystem.ToRGB(XYZColor(X, Y, Z));

	return Spectrum(
			(whiteRGB.c[0] != 0.f) ? rgb.c[0] / whiteRGB.c[0] : 0.f,
			(whiteRGB.c[1] != 0.f) ? rgb.c[1] / whiteRGB.c[1] : 0.f,
			(whiteRGB.c[2] != 0.f) ? rgb.c[2] / whiteRGB.c[2] : 0.f);
}
