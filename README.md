Running A1,A2,A3,A4,A5 sequentially will give the final result of 3D-SIM reconstruction. 
However, since the resulting datasets by A1 and A2, or the PSF_256Pixels.mat and Sample1000Tubule_3D_256Pixels.mat, have already been in the folder, running A3-A5, or A4-A5 if the wide-field image is not needed for comparison, is enough.

The running time of A4 is a little long, since 3D fft is used repetitively, mimicking the z-scanning process during the 3D-SIM imaging process. 

2026.06.13
A6, A7 and A8 can be run separately, given that 256Pixels.mat and Sample1000Tubule_3D_256Pixels.mat was generated and in folder
