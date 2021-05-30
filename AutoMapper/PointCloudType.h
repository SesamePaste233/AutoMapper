#pragma once
#include "BaseType.h"
#include <unordered_map>


namespace types {
	typedef struct Point3DRGBA{
		float X;
		float Y;
		float Z;

		uchar r, g, b, a;
	}CloudPoint;

	typedef std::unordered_map<int, CloudPoint> PointCloud;


}