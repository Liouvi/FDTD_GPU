#include "curl.cuh"
#include "geometry.cuh"
#include <cmath>

std::vector<float> flatten(const std::vector<std::vector<float>> &orig)
{   
    std::vector<float> ret;
    for(const auto &v: orig)
        ret.insert(ret.end(), v.begin(), v.end());                                                                                         
    return ret;
}   


curl::curl(geometry& GEO) {

	//Ez Mode update coefficients

	resolution = GEO.resolution;
	tstep = GEO.tstep;

	
	xstep = 1e-6 / resolution;

	// Coefficients


	// Update Hx coeff

}