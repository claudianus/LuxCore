/***************************************************************************
 * Copyright 1998-2020 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxCoreRender.                                   *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 *   You may not use this file except in compliance with the License.       *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 *   under the License is distributed on an "AS IS" BASIS,                 *
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 *   See the License for the specific language governing permissions and    *
 *   limitations under the License.                                         *
 ***************************************************************************/

// ObjC imports FIRST (see metaldevice.mm for the reason)
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#if defined(__APPLE__) && !defined(LUXRAYS_DISABLE_METAL)

#include "luxrays/devices/metalintersectiondevice.h"
#include "luxrays/devices/metalrtaccel.h"
#include "luxrays/accelerators/mbvhaccel.h"

#include <cstdlib>

using namespace std;

namespace luxrays {

//------------------------------------------------------------------------------
// MetalIntersectionDevice
//------------------------------------------------------------------------------

MetalIntersectionDevice::MetalIntersectionDevice(
	ContextConstRef context,
	MetalDeviceDescriptionConstRef desc,
	const size_t devIndex
) :
	Device(context, devIndex), MetalDevice(context, desc, devIndex),
	HardwareIntersectionDevice(), kernel(nullptr) {
}

MetalIntersectionDevice::~MetalIntersectionDevice() {
}

void MetalIntersectionDevice::SetDataSet(DataSetSPtr newDataSet) {
	IntersectionDevice::SetDataSet(newDataSet);

	if (dataSet) {
		const AcceleratorType accelType = dataSet->GetAcceleratorType();
		if (accelType != ACCEL_AUTO) {
			accel = dataSet->GetAccelerator(accelType);
		} else {
			// Mirror the OpenCL device selection: the MBVH accelerator +
			// kernel handles both instanced meshes and motion blur (it is
			// built through the generic HardwareDevice interface, so it
			// works on Metal exactly like the plain BVHKernel does).
			//
			// When native HWRT isn't disabled, always pick MBVH even for
			// non-instanced scenes: MetalRTKernel needs the MBVH leaf
			// metadata (uniqueLeafs/bvhLeafs) to build the Metal
			// acceleration structures, and the software MBVH kernel is
			// also a correct fallback.
			const char *hwrtEnv = getenv("LUXRAYS_METAL_HWRT");
			const bool hwrtEnabled = !(hwrtEnv && string(hwrtEnv) == "0");
			if (hwrtEnabled || dataSet->RequiresInstanceSupport() ||
					dataSet->RequiresMotionBlurSupport())
				accel = dataSet->GetAccelerator(ACCEL_MBVH);
			else
				accel = dataSet->GetAccelerator(ACCEL_BVH);
		}
	}
}

void MetalIntersectionDevice::Update() {
	kernel->Update(dataSet);
}

void MetalIntersectionDevice::Start() {
	MetalDevice::Start();

	// Native Metal HWRT first (MTLAccelerationStructure + intersector):
	// needs the MBVH accelerator for its leaf metadata. Falls back to the
	// software kernel (translated MBVH/BVH via cl2msl) when unsupported.
	kernel = nullptr;
	if (accel->GetType() == ACCEL_MBVH)
		kernel = NewMetalRTKernelIfPossible(*this,
				static_cast<const MBVHAccel &>(*accel));

	if (!kernel) {
		// Compile required kernel (BVHKernel/MBVHKernel: the same hardware
		// intersection kernel the OpenCL devices use - vertex/node buffer
		// packing + Accelerator_Intersect_RayBuffer through the cl2msl
		// pipeline)
		kernel = accel->NewHardwareIntersectionKernel(*this);
	}
}

void MetalIntersectionDevice::Stop() {
	kernel.reset();

	MetalDevice::Stop();
}

void MetalIntersectionDevice::EnqueueTraceRayBuffer(HardwareDeviceBuffer *rayBuff,
			HardwareDeviceBuffer *rayHitBuff,
			const unsigned int rayCount) {
	// Enqueue the intersection kernel
	kernel->EnqueueTraceRayBuffer(rayBuff, rayHitBuff, rayCount);
	statsTotalDataParallelRayCount += rayCount;
}

}

#endif
// vim: autoindent noexpandtab tabstop=4 shiftwidth=4
