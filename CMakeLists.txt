project(PhD C)
enable_language(C)
#sudo apt-get install libnifti-dev libfftw3-dev

cmake_minimum_required(VERSION 2.8.7)
SET(PHD_PACKAGE_VERSION_MAJOR 1)
SET(PHD_TOOLKIT_PACKAGE_VERSION_MINOR 1)
SET(PHD_TOOLKIT_PACKAGE_VERSION_PATCH 1)
SET(PHD_TOOLKIT_RELEASE_DATE "20160801")

SET( CMAKE_C_FLAGS  "-g" )


SET(PHD_VERSION_FULL "${PHD_PACKAGE_VERSION_MAJOR}.${PHD_PACKAGE_VERSION_MINOR}.${PHD_PACKAGE_VERSION_PATCH}-${PHD_RELEASE_DATE}")

set(CMAKE_INSTALL_PREFIX /usr/local/ CACHE PATH "Path for install directory.")

###############
#   Details   # 
###############
SET(CPACK_PACKAGE_NAME "")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MINC PET image processing tools")
SET(CPACK_PACKAGE_VENDOR "Brain Imaging Center, Montreal Neurological Institute, McGill University")
SET(CPACK_PACKAGE_VERSION_MAJOR ${MINC_TOOLKIT_PACKAGE_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${MINC_TOOLKIT_PACKAGE_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${MINC_TOOLKIT_PACKAGE_VERSION_PATCH})
SET(CPACK_PACKAGE_CONTACT "Thomas Funck <thomas.funck@mail.mcgill.ca>")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "pet-toolkit-${MINC_TOOLKIT_VERSION_FULL}")
SET(CPACK_MONOLITHIC_INSTALL ON)
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.txt")
SET(CPACK_PACKAGE_RELOCATABLE OFF)

####################################
# Set MINC2 installation directory #
####################################
SET(MINC_INSTALL_DIR /usr/local CACHE PATH "Path for MINC prefix, prepended onto install directories.")
set(MINC_INCLUDE_DIR ${MINC_INSTALL_DIR}/include/)
set(MINC_LIB_DIR ${MINC_INSTALL_DIR}/lib/)
set(MINC2_LIB  minc2 hdf5 netcdf z znz dl rt)
set(BICPL_LIB  bicpl  minc2 hdf5 netcdf niftiio z znz dl rt)
link_directories(${MINC_LIB_DIR})

include_directories(include)

###############################
# Compile minc_helper library #
###############################
set(minc_helper_src src/minc_helper/minc_helper-1.2.c ${MINC_INCLUDE_DIR}/minc2.h)
add_library(minc_helper STATIC ${minc_helper_src})
target_link_libraries(minc_helper ${MINC2_LIB} )
target_include_directories(minc_helper PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(minc_helper PROPERTIES COMPILE_FLAGS "-std=c99")

#############################
# Compile minc_fftw library #
#############################
set(minc_fftw_src src/minc_fftw/minc_fftw-1.0.c include/minc_fftw-1.0.h ${MINC_INCLUDE_DIR}/minc2.h include/minc_helper.h)
add_library(minc_fftw STATIC ${minc_fftw_src})
target_link_libraries(minc_fftw fftw3 ${MINC2_LIB} minc_helper )
target_include_directories(minc_fftw PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(minc_fftw PROPERTIES COMPILE_FLAGS "-std=c99")

##################
# Compile idSURF #
##################
set(idSURF_SRC src/idSURF/idSURF-4.1.1.c src/idSURF/gaussianiir3d.c   ${MINC_INCLUDE_DIR}/minc2.h include/minc_helper.h)
add_executable(idSURF ${idSURF_SRC})
target_link_libraries(idSURF  minc_helper ${MINC2_LIB} m pthread )
target_include_directories(idSURF PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(idSURF PROPERTIES COMPILE_FLAGS "-static -std=c99")

##########################
# Compile mincgroupstats #
##########################
set(mincgroupstats_SRC  src/mincgroupstats/mincgroupstats-4.0.c  ${MINC_INCLUDE_DIR}/minc2.h include/minc_helper.h)
add_executable(mincgroupstats ${mincgroupstats_SRC})
target_link_libraries(mincgroupstats  minc_helper ${MINC2_LIB} m pthread )
target_include_directories(mincgroupstats PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(mincgroupstats PROPERTIES COMPILE_FLAGS "-static -std=c99")

###############
# Compile gtm #
###############
set(gtm_SRC  src/gtm/gtm-3.2.c src/gtm/gaussj.c src/nrutil/nrutil.c include/nrutil.h  ${MINC_INCLUDE_DIR}/minc2.h include/minc_fftw-1.0.h include/minc_helper.h)
add_executable(gtm ${gtm_SRC})
target_link_libraries(gtm minc_fftw minc_helper ${MINC2_LIB} m pthread fftw3 )
target_include_directories(gtm PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(gtm PROPERTIES COMPILE_FLAGS "-static -std=c99")

###############
# Compile srv #
###############
set(srv_SRC  src/srv/srv-3.0.c src/srv/rrefTest.c src/srv/my_interpolate_sphere.c   src/srv/ray_casting.c ${MINC_INCLUDE_DIR}/minc2.h include/minc_helper.h)
add_executable(srv ${srv_SRC})
target_link_libraries(srv  minc_helper ${BICPL_LIB} m pthread )
target_include_directories(srv PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(srv PROPERTIES COMPILE_FLAGS "-static -std=c99")

#####################
# Compile mincsplit #
#####################
set(mincsplit_SRC  src/mincsplit/mincsplit.c   include/minc_helper.h  ${MINC_INCLUDE_DIR}/minc2.h )
add_executable(mincsplit ${mincsplit_SRC})
target_link_libraries(mincsplit minc_helper  ${MINC2_LIB} m pthread)
target_include_directories(mincsplit PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(mincsplit PROPERTIES COMPILE_FLAGS "-static -std=c99")

#####################
# Compile surf_dist #
#####################
set(surf_dist_SRC  src/surf_dist/surf_dist-5.1.c src/srv/my_interpolate_sphere.c    include/minc_helper.h  ${MINC_INCLUDE_DIR}/minc2.h )
add_executable(surf_dist ${surf_dist_SRC})
target_link_libraries(surf_dist minc_helper  ${BICPL_LIB} m pthread)
target_include_directories(surf_dist PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(surf_dist PROPERTIES COMPILE_FLAGS "-static -std=c99")

#######################
# Compile surf_defrag #
#######################
set(surf_defrag_SRC  src/surf_defrag/surf_defrag-1.1.c src/srv/my_interpolate_sphere.c    include/minc_helper.h  ${MINC_INCLUDE_DIR}/minc2.h )
add_executable(surf_defrag ${surf_defrag_SRC})
target_link_libraries(surf_defrag minc_helper  ${BICPL_LIB} m pthread)
target_include_directories(surf_defrag PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(surf_defrag PROPERTIES COMPILE_FLAGS "-static -std=c99")

###############
# Compile tka #
###############
set(tka_SRC  src/tka/1.0/pet_input.c src/tka/1.0/pet_core.c src/tka/1.0/pet_analyze.c src/tka/1.0/llsqwt.c src/tka/1.0/tracerDatabase.c include/llsqwt.h include/mtga.h include/nrutil.h include/minc_helper.h  ${MINC_INCLUDE_DIR}/minc2.h )
add_executable(tka ${tka_SRC})
target_link_libraries(tka minc_helper  ${MINC2_LIB} m pthread)
target_include_directories(tka PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(tka PROPERTIES COMPILE_FLAGS "-static -std=c99")

######################
# Compile minc2ascii #
######################
set(minc2ascii_SRC  src/minc2ascii/minc2ascii.c  include/minc_helper.h  ${MINC_INCLUDE_DIR}/minc2.h )
add_executable(minc2ascii ${minc2ascii_SRC})
target_link_libraries(minc2ascii minc_helper  ${MINC2_LIB} m pthread)
target_include_directories(minc2ascii PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(minc2ascii PROPERTIES COMPILE_FLAGS "-static -std=c99")


##########################
# Compile crop_from_mask #
##########################
set(crop_from_mask_SRC  src/crop_from_mask/crop_from_mask.c  include/minc_helper.h  ${MINC_INCLUDE_DIR}/minc2.h )
add_executable(crop_from_mask ${crop_from_mask_SRC})
target_link_libraries(crop_from_mask minc_helper  ${MINC2_LIB} m pthread)
target_include_directories(crop_from_mask PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(crop_from_mask PROPERTIES COMPILE_FLAGS "-static -std=c99")


######################
# Compile tac_volume #
######################
set(tac_volume_SRC  src/tac_volume/tac_volume-1.0.c  include/minc_helper.h  ${MINC_INCLUDE_DIR}/minc2.h )
add_executable(tac_volume ${tac_volume_SRC})
target_link_libraries(tac_volume minc_helper  ${MINC2_LIB} m pthread)
target_include_directories(tac_volume PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(tac_volume PROPERTIES COMPILE_FLAGS "-static -std=c99")


##########################
# Compile minc_integrate #
#########################
set(minc_integrate_SRC  src/minc_integrate/minc_integrate-1.0.c  include/minc_helper.h  ${MINC_INCLUDE_DIR}/minc2.h )
add_executable(minc_integrate ${minc_integrate_SRC})
target_link_libraries(minc_integrate minc_helper  ${MINC2_LIB} m pthread)
target_include_directories(minc_integrate PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(minc_integrate PROPERTIES COMPILE_FLAGS "-static -std=c99")


############################
# Compile join_mask_labels #
############################
set(join_mask_labels_SRC  src/join_mask_labels/join_mask_labels.c  include/minc_helper.h  ${MINC_INCLUDE_DIR}/minc2.h )
add_executable(join_mask_labels ${join_mask_labels_SRC})
target_link_libraries(join_mask_labels minc_helper  ${MINC2_LIB} m pthread)
target_include_directories(join_mask_labels PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(join_mask_labels PROPERTIES COMPILE_FLAGS "-static -std=c99")

##########################
# Compile nl-anisotropic #
##########################
file(GLOB nl-anisotropic_SRC src/nl-anisotropic/anisotropic_diffusion-3.0.c   src/nl-anisotropic/nl-anisotropic.c src/nl-anisotropic/3x3-C/dsyevh3.c src/nl-anisotropic/3x3-C/dsyevd3.c src/nl-anisotropic/3x3-C/dsyevc3.c src/nl-anisotropic/3x3-C/dsyevq3.c src/nl-anisotropic/3x3-C/dsytrd3.c src/nl-anisotropic/3x3-C/dsyev2.c src/nl-anisotropic/3x3-C/slvsec3.c  include/nl-anisotropic.h  include/minc_helper.h ${MINC_INCLUDE_DIR}/minc2.h )
add_executable(nl-anisotropic ${nl-anisotropic_SRC})
message(${nl-anisotropic_SRC})
target_link_libraries(nl-anisotropic minc_helper  ${MINC2_LIB} m pthread)
target_include_directories(nl-anisotropic PUBLIC ${MINC_INCLUDE_DIR} "include/3x3-C/" )
set_target_properties(nl-anisotropic PROPERTIES COMPILE_FLAGS "-static -std=c99")

###################
# Compile wm_dist #
###################
set(wm_dist_SRC  src/wm_dist/wm_dist-3.0.c  include/minc_helper.h  ${MINC_INCLUDE_DIR}/minc2.h )
add_executable(wm_dist ${wm_dist_SRC})
target_link_libraries(wm_dist minc_helper  ${MINC2_LIB} m pthread)
target_include_directories(wm_dist PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(wm_dist PROPERTIES COMPILE_FLAGS "-static -std=c99  -O0 -ffast-math  -fauto-inc-dec   -fbranch-count-reg -fcombine-stack-adjustments    -fcompare-elim  -fcprop-registers -fdce  -fdefer-pop -fdelayed-branch -fdse -fforward-propagate -fguess-branch-probability -fif-conversion2   -fif-conversion -finline-functions-called-once -fipa-pure-const -fipa-profile   -fipa-reference -fmerge-constants   -fmove-loop-invariants -freorder-blocks  -fsplit-wide-types   -ftree-bit-ccp -ftree-ccp -ftree-ch  -ftree-copy-prop  -ftree-dce -ftree-dominator-opts -ftree-dse  -ftree-forwprop -ftree-fre  -ftree-phiprop -ftree-sink    -ftree-sra -ftree-pta -ftree-ter -funit-at-a-time -fthread-jumps  -falign-functions  -falign-jumps  -falign-loops  -falign-labels -fcaller-saves  -fcrossjumping -fcse-follow-jumps  -fcse-skip-blocks  -fdelete-null-pointer-checks -fdevirtualize    -fexpensive-optimizations  -fgcse  -fgcse-lm      -finline-small-functions   -findirect-inlining    -fipa-cp    -fipa-sra          -foptimize-sibling-calls -fpartial-inlining  -fpeephole2   -freorder-blocks-and-partition -freorder-functions  -frerun-cse-after-loop  -fsched-interblock  -fsched-spec   -fschedule-insns  -fschedule-insns2 -fstrict-aliasing -fstrict-overflow  -ftree-builtin-call-dce  -ftree-switch-conversion -ftree-pre  -ftree-vrp     -fpeel-loops  -fipa-cp-clone")
#        
###########################
# Set install directories #
###########################
set(progs gtm idSURF mincgroupstats srv mincsplit surf_dist surf_defrag tka minc2ascii crop_from_mask join_mask_labels nl-anisotropic wm_dist tac_volume minc_integrate)
foreach(prog ${progs}) 
	install(TARGETS ${prog} DESTINATION bin)
endforeach(prog)

set(libs minc_helper minc_fftw)
foreach(lib ${libs})
	install(TARGETS ${lib} DESTINATION lib)
endforeach(lib)

