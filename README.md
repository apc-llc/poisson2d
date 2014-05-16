## A solver for 2D Poisson problem with Dirichlet or Neumann boundary conditions

### Building

Two possible library backends for FFT are supported: FFTW and FFTW-compatible interfaces of Intel MKL.

#### FFTW

```
$ mkdir build_fftw
$ cd build_fftw
$ cmake -DHAVE_FFTW=ON ..
```

#### MKL

```
$ mkdir build_mkl
$ cd build_mkl
$ cmake -DHAVE_FFTW_MKL=ON -DHAVE_FFTW_MKL_HEADER=/opt/intel/composer_xe_2013.3.163/mkl/include/fftw/ -DHAVE_FFTW_MKL_THREADS_LIBRARY=/opt/intel/composer_xe_2013.3.163/mkl/lib/intel64/ ..
```

### Test run

```
$ ./poisson2d_fft 4094 4094
Solve 2D Poisson equation
with Dirichlet b.c. by X,
and Neumann b.c. by Y:
Lx = f in D, phi = g on dD.

Method: 1d fft + shutter

Solver ncores = 1

Init time = 1.032101
Solver time = 0.902502
residual min = -0.025945, max = 0.025945, sum = 0.000104
Check time = 1.060190
Output time = 0.148881
Deinit time = 0.008483
```

### Visualize with GrADS

```
$ perl ./poisson2d_fft_phi1.pl 

Grid Analysis and Display System (GrADS) Version 2.0.a9
Copyright (c) 1988-2010 by Brian Doty and the
Institute for Global Environment and Society (IGES)
GrADS comes with ABSOLUTELY NO WARRANTY
See file COPYRIGHT for more information

Config: v2.0.a9 little-endian readline printim grib2 netcdf hdf4-sds hdf5 opendap-grids geotiff shapefile
Issue 'q config' command for more detailed configuration information
GX Package Initialization: Size = 11 8.5 
Running in Batch mode
Warning:  X axis labels overridden by SET XAXIS.
   Labels may not reflect correct scaling for dimensions or data.
Warning:  Y axis labels overridden by SET YAXIS.
   Labels may not reflect correct scaling for dimensions or data.

$ perl ./poisson2d_fft_phi2.pl 

Grid Analysis and Display System (GrADS) Version 2.0.a9
Copyright (c) 1988-2010 by Brian Doty and the
Institute for Global Environment and Society (IGES)
GrADS comes with ABSOLUTELY NO WARRANTY
See file COPYRIGHT for more information

Config: v2.0.a9 little-endian readline printim grib2 netcdf hdf4-sds hdf5 opendap-grids geotiff shapefile
Issue 'q config' command for more detailed configuration information
GX Package Initialization: Size = 11 8.5 
Running in Batch mode
Warning:  X axis labels overridden by SET XAXIS.
   Labels may not reflect correct scaling for dimensions or data.
Warning:  Y axis labels overridden by SET YAXIS.
   Labels may not reflect correct scaling for dimensions or data.
```

