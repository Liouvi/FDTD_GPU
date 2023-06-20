# FDTD_GPU
Parallel implementation of FDTD algorithm.

You can change simulation parameters in FDTD2D.cu. 

In order to compile and run:

```
nvcc -o FDTD FDTD2D.cu dataacc.cu curl.cu geometry.cu -rdc=true --expt-relaxed-constexpr && . /FDTD
```
