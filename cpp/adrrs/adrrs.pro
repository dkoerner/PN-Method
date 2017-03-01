TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    src/util/bitmap.cpp \
    src/util/data.cpp \
    src/util/field.cpp \
    src/util/string.cpp \
    src/util/threadpool.cpp \
    src/util/voxelgrid.cpp \
    src/Camera.cpp \
    src/volume.cpp \
    src/main.cpp \
    src/util/dda3d.cpp \
    src/light.cpp \
    src/scene.cpp \
    src/integrator.cpp \
    src/pncache.cpp \
    src/wedge.cpp \
    src/util/sh.cpp


include(deployment.pri)
qtcAddDeployment()


INCLUDEPATH += include


win32
{
	#eigen
	INCLUDEPATH += "c:/libs-msvc2013/include/eigen3"

	# ilmbase+openexr ----
	QMAKE_LIBDIR += "c:/libs-msvc2013/lib"
	INCLUDEPATH += "c:/libs-msvc2013/include"

	CONFIG(debug, debug|release) {
		LIBS += IlmImf-2_1d.lib Iex-2_1d.lib IlmThread-2_1d.lib Imath-2_1d.lib Halfd.lib
	}
	CONFIG(release, debug|release) {
		LIBS += IlmImf-2_1.lib Iex-2_1.lib IlmThread-2_1.lib Imath-2_1.lib Half.lib
	}

	# zlib ----
	LIBS+= "zlibstatic.lib"

	#houio ----
	INCLUDEPATH += "c:/libs-msvc2013/include/houio"
	CONFIG(debug, debug|release) {
				LIBS+= "houiod.lib"
	}
	CONFIG(release, debug|release) {
				LIBS+= "houio.lib"
	}


	#ceres-solver (for fitting test) ---
	#CONFIG(debug, debug|release) {
	#	LIBS+= "c:/libs-install-msvc2013/ceres-solver-debug/lib/ceres-debug.lib"
	#	INCLUDEPATH += "c:/libs-install-msvc2013/ceres-solver-debug/include"
	#	INCLUDEPATH += "c:/libs-install-msvc2013/ceres-solver-debug/include"
	#	INCLUDEPATH += "c:/libs-install-msvc2013/ceres-solver-debug/include/ceres/internal/miniglog"
	#}
	#CONFIG(release, debug|release) {
	#	LIBS+= "c:/libs-install-msvc2013/ceres-solver/lib/ceres.lib"
	#	INCLUDEPATH += "c:/libs-install-msvc2013/ceres-solver/include"
	#	INCLUDEPATH += "c:/libs-install-msvc2013/ceres-solver/include"
	#	INCLUDEPATH += "c:/libs-install-msvc2013/ceres-solver/include/ceres/internal/miniglog"
	#}

	#glog ---
	#INCLUDEPATH += "c:/libs-install-msvc2013/glog/include"
	#LIBS+= "c:/libs-install-msvc2013/glog/lib/glog.lib"

	#gflags ---
	#INCLUDEPATH += "c:/libs-install-msvc2013/gflags/include"
	#LIBS+= "c:/libs-install-msvc2013/gflags/lib/gflags_static.lib"
	#LIBS+= Shlwapi

	#DEFINES += MAX_LOG_LEVEL=-1
	#DEFINES += GLOG_STATIC=-1

	##dlib ------
	#INCLUDEPATH += "c:/libs-msvc2013/include"
	#QMAKE_LIBDIR += c:/libs-msvc2013/lib
	#CONFIG(release, debug|release)
	#{
	#	LIBS += dlib.lib
	#}

}

HEADERS += \
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
    include/util/bitmap.h \
    include/util/data.h \
    include/util/field.h \
    include/util/string.h \
    include/util/threadpool.h \
    include/util/timer.h \
    include/util/voxelgrid.h \
    include/camera.h \
    include/volume.h \
    include/util/dda3d.h \
    include/light.h \
    include/scene.h \
    include/integrator.h \
    include/pncache.h \
    include/cache.h \
    include/util/wedge.h \
	include/util/sh.h
