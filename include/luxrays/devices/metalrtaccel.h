/***************************************************************************
 * Copyright 1998-2020 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxCoreRender.                                   *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

#ifndef _LUXRAYS_METALRTACCEL_H
#define	_LUXRAYS_METALRTACCEL_H

// Native Metal hardware ray tracing (MTLAccelerationStructure) intersection
// kernel. This is the HWRT counterpart of the software MBVH kernel: it builds
// a primitive acceleration structure per unique mesh and a single instance
// acceleration structure over the MBVH leaf references, then resolves rays
// with raytracing::intersector instead of traversing BVH nodes in software.
//
// Everything is implemented in metalrtaccel.mm (ObjC++); this header only
// exposes a factory so no ObjC types leak into the C++ headers.

#include "luxrays/core/hardwareintersectiondevice.h"

namespace luxrays {

class MBVHAccel;

// Returns a HardwareIntersectionKernel backed by native Metal ray tracing
// when the device and the scene support it, nullptr otherwise (caller is
// expected to fall back to accel->NewHardwareIntersectionKernel()).
//
// Currently unsupported (=> nullptr, software fallback):
//  - devices without MTLDevice.supportsRaytracing
//  - motion blur (MBVH leaves carrying a motion system)
//  - MBVH leaves built over more than one mesh
//
// The LUXRAYS_METAL_HWRT environment variable can force-disable the path
// with "0" (any other value or unset leaves it enabled).
HardwareIntersectionKernelUPtr NewMetalRTKernelIfPossible(
	HardwareIntersectionDevice &device, const MBVHAccel &accel);

}

#endif	/* _LUXRAYS_METALRTACCEL_H */
// vim: autoindent noexpandtab tabstop=4 shiftwidth=4
