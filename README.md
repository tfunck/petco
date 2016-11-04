# petco
## About
### Miscallaneous programs for analyzing (mostly) PET images  
**idSURF**: voxel-wise partial-volume correction of PET images [1]
**gtm**: geometric transfer matrix partial-volume correction method [2]
**srv**: create super-resolution volumetric image based on pair of surfaces  
**tka**: perform tracer kinetic analysis on PET images (patlak[3], logan[4], suvr)  
**surf_dist**: calculate geodesic surface distance from a roi [5]

## Installation:  
mkdir build  
cd build   
ccmake -D MINC_INSTALL_DIR="*path to your minc installation*" -D CMAKE_INSTALL_PREFIX="*install path*" .. 
make  
make install  

## References:
[1] Funck, T., Paquette, C., Evans, A. & Thiel, A. Surface-based partial-volume correction for high-resolution PET. Neuroimage 102 Pt 2, 674–87 (2014).
[2] Rousset, O. G., Ma, Y. & Evans, A. C. Correction for Partial Volume Effects in PET : Principle and Validation. Journal of Nuclear Medicine. 39, 904–911 (1998).
[3] Patlak, C. S., Blasberg, R. G. & Fenstermacher, J. D. Graphical evaluation of blood-to-brain transfer constants from multiple-time uptake data. J. Cereb. Blood Flow Metab. 3, 1–7 (1983).
[4] Logan, J. et al. Distribution Volume Ratios Without Blood Sampling from Graphical Analysis of PET Data. J. Cereb. Blood Flow Metab. 16, 834–840 (1996).
[5] Funck, T. et al. Assessing neuronal density in peri-infarct cortex with PET: Effects of cortical topology and partial volume correction. Hum. Brain Mapp. (2016). doi:10.1002/hbm.23363
