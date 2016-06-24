cmake_minimum_required (VERSION 2.8.9)
project (HPMVS)

ADD_DEFINITIONS(
    -std=c++11
    -O3
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

#find_package(Loccom REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(JPEG REQUIRED)
find_package(Gflags REQUIRED)
add_definitions(-DHPMVS_GFLAGS_NAMESPACE=${GFLAGS_NAMESPACE})
find_package(Glog REQUIRED)

ADD_SUBDIRECTORY(thirdLibs/nlopt-2.4.2)
ADD_SUBDIRECTORY(thirdLibs/stlplus3)

set(MY_SRC
	src/hpmvs/NVMReader.cpp
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


include_directories(
	"${PROJECT_SOURCE_DIR}/include"
	${EIGEN3_INCLUDE_DIR}
	${JPEG_INCLUDE_DIR}
	${GFLAGS_INCLUDE_DIRS}
	${GLOG_INCLUDE_DIRS}
	${CMAKE_SOURCE_DIR}/thirdLibs/nlopt-2.4.2/api
	${CMAKE_SOURCE_DIR}/thirdLibs/cimg
	${CMAKE_SOURCE_DIR}/thirdLibs/
	${CMAKE_SOURCE_DIR}/thirdLibs/stlplus3
	)


add_executable (hpmvs src/main.cpp ${MY_SRC})
target_link_libraries(hpmvs ${JPEG_LIBRARIES} ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES} nlopt stlplus3)

install(TARGETS hpmvs RUNTIME DESTINATION bin)
