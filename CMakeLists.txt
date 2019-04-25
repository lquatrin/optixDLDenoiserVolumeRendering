if(GLUT_FOUND AND OPENGL_FOUND)
  include_directories(${GLUT_INCLUDE_DIR})
  add_definitions(-DGLUT_FOUND -DGLUT_NO_LIB_PRAGMA)

  OPTIX_add_sample_executable( optixDLDenoiserVolumeRendering
    optixdenvolrend.cpp
    pvm.cpp
    pvm.h
    rawloader.cpp
    rawloader.h
    reader.cpp
    reader.h
    transferfunction.cpp
    transferfunction.h
    transferfunction1d.cpp
    transferfunction1d.h
    volume.cpp
    volume.h
    tutorial.h
    raymarching.cu
    box.cu
	random.h
	voldatamanager.h
	voldatamanager.cpp
  )
else()
  # GLUT or OpenGL not found
  message("Disabling optixTutorial, which requires glut and opengl.")
endif()
