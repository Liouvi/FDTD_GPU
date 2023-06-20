#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <vector>
#include <cmath>
#include <math.h>

using namespace std;

class geometry
{
public:
	int resolution;
	float xdim;
	float ydim;
	int xsize;
	int ysize;
	float tstep;

	
	geometry(int reso_in, float xdim_in, float ydim_in, float timestep);

	void block(float xcenter, float ycenter, float length, float width, float eps, float mu);

	void disk(float xcenter, float ycenter, float radius, float eps, float mu);
	
	void tiltedblock(float theta, float dist, float eps, float mu);
};

#endif
