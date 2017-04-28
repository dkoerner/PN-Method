# the following links helped to get this to work:
# http://stackoverflow.com/questions/12266264/compiling-cuda-code-in-qt-creator-on-windows
# http://stackoverflow.com/questions/26205608/compile-cuda-file-error-runtime-library-mismatch-value-mdd-dynamicdebug-doe

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
			src/Camera.cpp \
			src/volume.cpp \
			src/shcache.cpp \
			src/util/string.cpp \
			src/util/bitmap.cpp \
			src/util/envmap.cpp \
			src/util/data.cpp \
			src/util/field.cpp \
			src/util/voxelgrid.cpp \
			src/util/cas.cpp \
			src/util/moexp.cpp \
			src/util/moexp2d.cpp

HEADERS += \
	include/camera.h \
	include/volume.h \
	include/shcache.h \
	include/cache.h \
	include/math/bbox.h \
	include/math/color.h \
	include/math/common.h \
	include/math/dpdf.h \
	include/math/frame.h \
	include/math/pcf.h \
	include/math/plf.h \
	include/math/ray.h \
	include/math/rng.h \
	include/math/transform.h \
	include/math/vector.h \
	include/util/string.h \
	include/util/bitmap.h \
	include/util/envmap.h \
	include/util/data.h \
	include/util/wedge.h \
	include/util/timer.h \
	include/util/field.h \
	include/util/voxelgrid.h \
	include/util/cas.h \
	include/util/moexp.h \
	include/util/moexp2d.h \
	include/cuda/cumath/cumath.h \
	include/cuda/cumath/Vec2.h \
	include/cuda/cumath/Vec2Algo.h \
	include/cuda/cumath/Vec3.h \
	include/cuda/cumath/Vec3Algo.h \
	include/cuda/cumath/Matrix22.h \
	include/cuda/cumath/Matrix33.h \
	include/cuda/cumath/Matrix33Algo.h \
	include/cuda/cumath/Matrix44.h \
	include/cuda/cumath/Matrix44Algo.h \
	include/cuda/cumath/Ray2.h \
	include/cuda/cumath/Ray3.h \
	include/cuda/CudaData.h \
	include/cuda/CudaPixelBuffer.h \
	include/cuda/CudaEnvMap.h \
	include/cuda/pathtracing.cu.h \
	include/cuda/CudaVoxelGrid.h \
	include/cuda/CudaVolume.h \
	include/cuda/CudaLight.h

include(deployment.pri)
qtcAddDeployment()

# Define output directories
DESTDIR = release
OBJECTS_DIR = release/obj
CUDA_OBJECTS_DIR = release/cuda

# CUDA settings <-- may change depending on your system
CUDA_SOURCES += vectorAddition.cu src/pt_2d_anisotropic.cu
CUDA_SDK = "c:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v7.5"
#CUDA_SDK = "C:/ProgramData/NVIDIA Corporation/NVIDIA GPU Computing SDK 4.2/C"   # Path to cuda SDK install
CUDA_DIR = "c:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v7.5"
#CUDA_DIR = "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v4.2"            # Path to cuda toolkit install
SYSTEM_NAME = x64         # Depending on your system either 'Win32', 'x64', or 'Win64'
SYSTEM_TYPE = 64            # '32' or '64', depending on your system
#CUDA_ARCH = compute_35           # Type of CUDA architecture, for example 'compute_10', 'compute_11', 'sm_10'
CUDA_ARCH = compute_20           # Type of CUDA architecture, for example 'compute_10', 'compute_11', 'sm_10'
NVCC_OPTIONS = --use_fast_math

# include paths
INCLUDEPATH +=  $$CUDA_DIR/include \
				$$CUDA_SDK/common/inc/ \
				$$CUDA_SDK/../shared/inc/


QMAKE_LIBDIR += "c:/libs-msvc2013/lib"
INCLUDEPATH += "c:/libs-msvc2013/include"

INCLUDEPATH += $$PWD/include
INCLUDEPATH += "c:/libs-msvc2013/include/eigen3"
INCLUDEPATH += "c:/libs-msvc2013/include/houio"
INCLUDEPATH += "c:/libs-msvc2013/include/boost-1_61"


# library directories
QMAKE_LIBDIR += $$CUDA_DIR/lib/$$SYSTEM_NAME \
				$$CUDA_SDK/common/lib/$$SYSTEM_NAME \
				$$CUDA_SDK/../shared/lib/$$SYSTEM_NAME
# Add the necessary libraries
CUDA_LIBS = -lcuda -lcudart
LIBS += $$CUDA_LIBS
LIBS += "zlibstatic.lib"

# The following library conflicts with something in Cuda
#QMAKE_LFLAGS_RELEASE = /NODEFAULTLIB:msvcrt.lib
#QMAKE_LFLAGS_DEBUG   = /NODEFAULTLIB:msvcrtd.lib

# The following makes sure all path names (which often include spaces) are put between quotation marks
CUDA_INC = $$join(INCLUDEPATH,'" -I"','-I"','"')

# MSVCRT link option (static or dynamic, it must be the same with your Qt SDK link option)
MSVCRT_LINK_FLAG_DEBUG = "/MDd"
MSVCRT_LINK_FLAG_RELEASE = "/MD"

# Configuration of the Cuda compiler
CONFIG(debug, debug|release){
	# Debug settings
	LIBS+= "houiod.lib"
	LIBS += IlmImf-2_1d.lib Iex-2_1d.lib IlmThread-2_1d.lib Imath-2_1d.lib Halfd.lib
	# Debug mode
	cuda_d.input    = CUDA_SOURCES
	cuda_d.output   = $$CUDA_OBJECTS_DIR/${QMAKE_FILE_BASE}_cuda.obj
	cuda_d.commands = $$CUDA_DIR/bin/nvcc.exe -D_DEBUG $$NVCC_OPTIONS $$CUDA_INC $$CUDA_LIBS \
					  --machine $$SYSTEM_TYPE -arch=$$CUDA_ARCH \
					  --compile -cudart static -g -DWIN32 -D_MBCS \
					  -Xcompiler "/wd4819,/EHsc,/W3,/nologo,/Od,/Zi,/FS,/RTC1" \
					  -Xcompiler $$MSVCRT_LINK_FLAG_DEBUG \
					  -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
	cuda_d.dependency_type = TYPE_C
	QMAKE_EXTRA_COMPILERS += cuda_d
}
else {
	 # Release settings
	LIBS+= "houio.lib"
	LIBS += IlmImf-2_1.lib Iex-2_1.lib IlmThread-2_1.lib Imath-2_1.lib Half.lib
	cuda.input    = CUDA_SOURCES
	cuda.output   = $$CUDA_OBJECTS_DIR/${QMAKE_FILE_BASE}_cuda.obj
	cuda.commands = $$CUDA_DIR/bin/nvcc.exe $$NVCC_OPTIONS $$CUDA_INC $$CUDA_LIBS \
					--machine $$SYSTEM_TYPE -arch=$$CUDA_ARCH \
					--compile -cudart static -DWIN32 -D_MBCS \
					-Xcompiler "/wd4819,/EHsc,/W3,/nologo,/O2,/Zi,/FS" \
					-Xcompiler $$MSVCRT_LINK_FLAG_RELEASE \
					-c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
	cuda.dependency_type = TYPE_C
	QMAKE_EXTRA_COMPILERS += cuda
}
