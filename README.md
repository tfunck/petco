# petco
## About
### Miscallaneous programs for analyzing (mostly) PET images  
**idSURF: voxel-wise partial-volume correction of PET images  
**srv**: create super-resolution volumetric image based on pair of surfaces  
**tka**: perform tracer kinetic analysis on PET images (patlak, logan, suvr)  

## Installation:  
mkdir build  
cd build   
ccmake -D MINC_INSTALL_DIR=*"path to your minc installation"* -D CMAKE_INSTALL_PREFIX=*"path to which you want to install"* ..  
make  
make install  


