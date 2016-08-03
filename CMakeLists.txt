project(PhD C)
enable_language(C)


cmake_minimum_required(VERSION 2.8.7)
SET(PHD_PACKAGE_VERSION_MAJOR 1)
SET(PHD_TOOLKIT_PACKAGE_VERSION_MINOR 1)
SET(PHD_TOOLKIT_PACKAGE_VERSION_PATCH 1)
SET(PHD_TOOLKIT_RELEASE_DATE "20160801")

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
message(${MINC_INSTALL_DIR})
set(MINC_INCLUDE_DIR ${MINC_INSTALL_DIR}/include/)
set(MINC_LIB_DIR ${MINC_INSTALL_DIR}/lib/)
set(MINC2_LIB minc2 hdf5 netcdf z dl)
link_directories(${MINC_LIB_DIR})

include_directories(include)

#######################
# Compile minc_helper #
#######################
set(minc_helper_src src/minc_helper/minc_helper-1.2.c ${MINC_INCLUDE_DIR}/minc2.h)
add_library(minc_helper STATIC ${minc_helper_src})
target_link_libraries(minc_helper ${MINC2_LIB} )
target_include_directories(minc_helper PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(minc_helper PROPERTIES COMPILE_FLAGS "-std=c99")


##################
# Compile idSURF #
##################
set(idSURF_SRC src/idSURF/idSURF-4.2.c src/idSURF/gaussianiir3d.c  ${MINC_INCLUDE_DIR}/minc2.h include/minc_helper.h)
add_executable(idSURF ${idSURF_SRC})
target_link_libraries(idSURF  minc_helper ${MINC2_LIB} m pthread )
target_include_directories(idSURF PUBLIC ${MINC_INCLUDE_DIR})
set_target_properties(idSURF PROPERTIES COMPILE_FLAGS "-static -std=c99")
#SET(CMAKE_EXE_LINKER_FLAGS "-static")
message(${MINC2_LIB})

install(TARGETS idSURF DESTINATION bin)
install(TARGETS minc_helper DESTINATION lib)

