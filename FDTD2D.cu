#include "geometry.cuh"
#include "curl.cuh"
#include "constants.cuh"
#include "dataacc.cuh"
#include <cmath>
#include <string>
#include <iterator>
#include <cuda_fp16.h>

__global__ void update_Hz(float *Ex, float *Ey, float *Hz,float*mu,int xmax,int ymax,float xstep,float tstep, float speedoflight)
{
	// Curl
	float CurlEz;
	int i = xmax * blockIdx.y * blockDim.y + xmax * threadIdx.y + blockDim.x * blockIdx.x + threadIdx.x;
	if(i < xmax * ymax) {

			// Curl
	if (blockIdx.x * blockDim.x + threadIdx.x == xmax-1 && blockIdx.y * blockDim.y + threadIdx.y == ymax-1)
	{
		CurlEz=( (0 - Ey[i]) - (0 - Ex[i]) ) / xstep;
	}

	else if (blockIdx.x * blockDim.x + threadIdx.x == xmax-1 && blockIdx.y * blockDim.y + threadIdx.y != ymax-1) {
		CurlEz =  ( ( 0 - Ey[i] ) - (0 - Ex[i] ) ) /xstep;
	}

	else {
		CurlEz = ( (Ey[i + 1] - Ey[i]) - (Ex[i + xmax] - Ex[i]) ) / xstep;
	}

	//Update
	__syncthreads();
	
	Hz[i] = Hz[i] + (-speedoflight*tstep / mu[i]) * CurlEz; 

	}



}

__global__ void update_Dx(float *Hz, float *Dx,int xmax,int ymax,float xstep,float tstep)
{
	// Curl
	float CurlHx;
	unsigned int i = xmax * blockIdx.y * blockDim.y + xmax * threadIdx.y + blockDim.x * blockIdx.x + threadIdx.x;
	if(i < xmax * ymax) {
			if(blockIdx.y == 0 && threadIdx.y == 0) {
		CurlHx = (Hz[i] - 0) / xstep;
	}
	
	else {
		CurlHx = (Hz[i] - Hz[i - xmax]) / xstep;
	}

	// Dx Update

	__syncthreads();

	Dx[i] = Dx[i] + (speedoflight * tstep) * CurlHx;


	}


}


__global__ void update_Dy(float *Hz, float *Dy,int xmax,int ymax,float xstep,float tstep)
{
	// Curl
	float CurlHy;
	unsigned int i = xmax * blockIdx.y * blockDim.y + xmax * threadIdx.y + blockDim.x * blockIdx.x + threadIdx.x;

	if(i < xmax * ymax) {
		// Curl
	if(blockIdx.x == 0 && threadIdx.x == 0) {
		CurlHy = -(Hz[i] - 0);
	}




	else {
		CurlHy = -(Hz[i] - Hz[i - 1]) / xstep;
	}

	// Dy Update
	__syncthreads();
	Dy[i] = Dy[i] + (speedoflight * tstep) * CurlHy;

	}
	
}

__global__ void update_E(float *Dx, float *Dy, float *Ex, float *Ey, float* eps,int xmax,int ymax)
{
	unsigned int i = xmax * blockIdx.y * blockDim.y + xmax * threadIdx.y + blockDim.x * blockIdx.x + threadIdx.x;

	if(i < xmax*ymax) {
		Ex[i] = Dx[i] / eps[i];
		Ey[i] = Dy[i] / eps[i];
	}
	

}

__global__ void source(float *Hz,float freq, float maxtime,float timestep,int qtime,int xmax,int ymax)
{		
		if (blockIdx.x == 16 && blockIdx.y  == 16 && threadIdx.x > 11 && threadIdx.x < 21 && threadIdx.y < 21 && threadIdx.y > 11) {
			Hz[xmax * blockIdx.y * blockDim.y + xmax * threadIdx.y + blockDim.x * blockIdx.x + threadIdx.x] += 1 * sin(freq * timestep * qtime) * (1-exp(-qtime / 5));
		}
}


int main()
{
	// Geometry Defn and Other Definitions

	float xdim =16.0;
	float ydim = 16.0;
	int resolution = 64;
	float eps_si = pow(3.47,2);
	int maxtime = 2500;

	int xsize = xdim * resolution;
	float freq = 193.41449e12 * 7;
	float xstep = 1e-6 / resolution;
	float courantnumber = 0.5;
	float timestep = courantnumber * xstep / speedoflight;

	// Create Geometry Object

	geometry lattice {resolution,xdim, ydim,timestep};

	// Create E-H field grid and set the simulation functions

	curl Grid2D{lattice};
	
	// Main Loop
	
	int frame = 0;

	float *dE_x,*dE_y;
	float *dD_x,*dD_y;

	float *dH_z;

	float *dmu;
	float *deps;

	int *dxmax;
	int *dymax;

	float *dxstep;
	float *dtstep;
	float *dfreq;
	float *dmaxtime;
	float *dtimestep;
	float *mu;
	float *eps;

	float *E_x, *E_y, *H_z;
	float *D_x, *D_y;

	E_x = (float *)malloc(xsize*xsize * sizeof(float));
	E_y = (float *)malloc(xsize*xsize * sizeof(float));
	H_z = (float *)malloc(xsize*xsize * sizeof(float));
	D_x = (float *)malloc(xsize*xsize * sizeof(float));
	D_y = (float *)malloc(xsize*xsize * sizeof(float));
	mu = (float *)malloc(xsize*xsize * sizeof(float));
	eps = (float *)malloc(xsize*xsize * sizeof(float));
	D_y = (float *)malloc(xsize*xsize * sizeof(float));



	for (int i = 0; i < xsize*xsize; i++)
    	{
        E_x[i] = 0;
        E_y[i] = 0;
		H_z[i] = 0;
		eps[i] = 13;
		mu[i] = 1;
    	}

	cudaMalloc((void **)&dE_x, xsize*xsize * sizeof(float));
    cudaMalloc((void **)&dE_y, xsize*xsize * sizeof(float));

	cudaMalloc((void **)&dD_x, xsize*xsize * sizeof(float));
    cudaMalloc((void **)&dD_y, xsize*xsize * sizeof(float));

    cudaMalloc((void **)&dH_z, xsize*xsize * sizeof(float));

    cudaMalloc((void **)&dmu, xsize*xsize * sizeof(float));
    cudaMalloc((void **)&deps, xsize*xsize * sizeof(float));

	cudaMalloc((void **)&dymax,sizeof(int));
    cudaMalloc((void **)&dxmax,sizeof(int));

	cudaMalloc((void **)&dxstep,sizeof(float));
    cudaMalloc((void **)&dtstep,sizeof(float));
	cudaMalloc((void **)&dmaxtime,sizeof(float));
	cudaMalloc((void **)&dtimestep,sizeof(float));
	cudaMalloc((void **)&dfreq,sizeof(float));

	cudaMemcpy(dE_x, E_x, xsize*xsize* sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dE_y, E_y, xsize*xsize * sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(dD_x, D_x, xsize*xsize * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dD_y, D_y, xsize*xsize * sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(dH_z, H_z, xsize*xsize * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dmu, mu, xsize*xsize * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(deps, eps, xsize*xsize * sizeof(float), cudaMemcpyHostToDevice);
	


	int xmax = 1024;
	int ymax = 1024;

	cudaMemcpy(dxmax, &xmax, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dymax, &ymax,sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(dxstep, &xstep,sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dtstep, &timestep,sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dmaxtime, &maxtime, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dfreq, &freq, sizeof(float), cudaMemcpyHostToDevice);

	int BLOCK_SIZE=32;
    unsigned int grid_rows = (xmax + BLOCK_SIZE - 1) / BLOCK_SIZE + 1;
    unsigned int grid_cols = (ymax + BLOCK_SIZE - 1) / BLOCK_SIZE + 1;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);



	for (int qtime = 0; qtime < maxtime; qtime++) {

		update_Hz<<<dimGrid, dimBlock>>>(dE_x, dE_y, dH_z, dmu,xmax,ymax,xstep,timestep,speedoflight);
		cudaDeviceSynchronize();
		update_Dx<<<dimGrid, dimBlock>>>(dH_z, dD_x,xmax,ymax,xstep,timestep);
		cudaDeviceSynchronize();
		update_Dy<<<dimGrid, dimBlock>>>(dH_z, dD_y,xmax,ymax,xstep,timestep);
		cudaDeviceSynchronize();
		source<<<dimGrid, dimBlock>>>(dH_z,freq, maxtime,timestep,qtime,xmax,ymax);
		cudaDeviceSynchronize();
		update_E<<<dimGrid, dimBlock>>>(dD_x, dD_y, dE_x, dE_y, deps,xmax,ymax);
		cudaDeviceSynchronize();

		if (qtime % 10 == 0) {
			cudaDeviceSynchronize();
    		cudaMemcpy(H_z, dH_z, 1024 * 1024 * sizeof(float), cudaMemcpyDeviceToHost);

			string name = "outputs/Hz";
			name.append(to_string(frame));
			string exten = ".dat";


			std::ofstream output_file(name+exten);
			std::ostream_iterator<float> output_iterator(output_file, " ");
			std::copy(H_z, H_z +xsize*xsize, output_iterator);
			frame += 1;
			
	}
			cudaDeviceSynchronize();
		}

	}





	
