cmake_minimum_required (VERSION 2.8.9)
project (HPMVS)

ADD_DEFINITIONS(
    -std=c++11
    -O3
)

# list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

#find_package(Loccom REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(JPEG REQUIRED)
find_package(Gflags REQUIRED)
# add_definitions(-DHPMVS_GFLAGS_NAMESPACE=${GFLAGS_NAMESPACE})
find_package(Glog REQUIRED)

ADD_SUBDIRECTORY(thirdLibs/nlopt-2.4.2)
ADD_SUBDIRECTORY(thirdLibs/stlplus3)
ADD_SUBDIRECTORY(thirdLibs/nvmtools)

set(HPMVS_SRC
	src/hpmvs/doctree.cpp
	src/hpmvs/CellProcessor.cpp
	src/hpmvs/PatchOptimizer.cpp
	src/hpmvs/Patch3d.cpp
	src/hpmvs/Scene.cpp
	src/hpmvs/Image.cpp
	src/hpmvs/Camera.cpp
	)

find_package(OpenMP)
if (OPENMP_FOUND)
	    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# include directories
set(HPMVS_INCLUDE_DIRS
	${PROJECT_SOURCE_DIR}/include
	${EIGEN3_INCLUDE_DIR}
	${JPEG_INCLUDE_DIR}
	${GFLAGS_INCLUDE_DIRS}
	${GLOG_INCLUDE_DIRS}
	${NVMTOOLS_INCLUDE_DIR}
	${PROJECT_SOURCE_DIR}/thirdLibs/nlopt-2.4.2/api
	${PROJECT_SOURCE_DIR}/thirdLibs/cimg
	${PROJECT_SOURCE_DIR}/thirdLibs/
	${PROJECT_SOURCE_DIR}/thirdLibs/stlplus3
	)

include_directories(${HPMVS_INCLUDE_DIRS})

# dependency libraries
set(HPMVS_LIBRARIES
	${NVMTOOLS_LIBRARIES}
	${JPEG_LIBRARIES} 
	${GFLAGS_LIBRARIES} 
	${GLOG_LIBRARIES} 
	nlopt 
	stlplus3_hpmvs
	)

# build static library
add_library(HPMVS_static STATIC ${HPMVS_SRC})
target_link_libraries(HPMVS_static ${HPMVS_LIBRARIES})



# build the executable
add_executable (hpmvs src/main.cpp)
target_link_libraries(hpmvs HPMVS_static ${HPMVS_LIBRARIES})

# install path
install(TARGETS hpmvs RUNTIME DESTINATION bin)
