#include "geometry.cuh"
#include "constants.cuh"
#include <math.h>

	geometry::geometry(int reso_in, float xdim_in, float ydim_in, float timestep)
	{
		resolution = reso_in;
		xdim = xdim_in;
		ydim = ydim_in;
		xsize = xdim * resolution;
		ysize = ydim * resolution;
		tstep = timestep;
	}

