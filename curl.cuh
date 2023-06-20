#ifndef CURLV2_H
#define CURLV2_H
#include<algorithm>
#include<iostream>
#include<vector>
#include "constants.cuh"
#include <cmath>
#include "geometry.cuh"

class curl
{
public:
	std::vector<float> Ex,Ey,Ez;

	std::vector<float> Dx,Dy,Dz;

	std::vector<float> Hx,Hy,Hz;

	std::vector<float> CurlEx,CurlEy,CurlEz;

	std::vector<float> CurlHx,CurlHy,CurlHz;


	float tstep;
	int resolution;
	float xstep;
	float courantnumber = 0.5;

	curl(geometry& GEO);

	friend class geometry;

};


#endif
