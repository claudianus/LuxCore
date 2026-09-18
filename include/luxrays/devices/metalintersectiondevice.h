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

#ifndef _LUXRAYS_METALINTERSECTIONDEVICE_H
#define _LUXRAYS_METALINTERSECTIONDEVICE_H

#include "luxrays/devices/metaldevice.h"
#include "luxrays/core/hardwareintersectiondevice.h"

#if defined(__APPLE__) && !defined(LUXRAYS_DISABLE_METAL)

namespace luxrays {

//------------------------------------------------------------------------------
// MetalIntersectionDevice
//------------------------------------------------------------------------------

class MetalIntersectionDevice : public MetalDevice, public HardwareIntersectionDevice {
public:
	MetalIntersectionDevice(const Context & context,
		MetalDeviceDescriptionConstRef desc, const size_t devIndex);
	virtual ~MetalIntersectionDevice();

	virtual void SetDataSet(DataSetSPtr newDataSet);
	virtual void Start();
	virtual void Stop();

	//--------------------------------------------------------------------------
	// Data parallel interface: to trace a multiple rays (i.e. on the GPU)
	//--------------------------------------------------------------------------

	virtual void EnqueueTraceRayBuffer(HardwareDeviceBuffer *rayBuff,
			HardwareDeviceBuffer *rayHitBuff,
			const unsigned int rayCount);

	friend class Context;

protected:
	virtual void Update();

	HardwareIntersectionKernelUPtr kernel;
};

}

#endif

#endif	/* _LUXRAYS_METALINTERSECTIONDEVICE_H */
// vim: autoindent noexpandtab tabstop=4 shiftwidth=4
