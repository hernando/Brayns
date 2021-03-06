# Copyright (c) 2015-2016, EPFL/Blue Brain Project
# All rights reserved. Do not distribute without permission.
# Responsible Author: Cyrille Favreau <cyrille.favreau@epfl.ch>
#
# This file is part of Brayns <https://github.com/BlueBrain/Brayns>

# Configure the output directories. To allow IMPI to do its magic we
# will put *executables* into the (same) build directory, but tag
# mic-executables with ".mic". *libraries* cannot use the
# ".mic"-suffix trick, so we'll put libraries into separate
# directories (names 'intel64' and 'mic', respectively)

set(BRAYNSOSPRAYPLUGIN_ISPC_SOURCES
  ispc/camera/ClippedPerspectiveCamera.ispc
  ispc/geometry/ExtendedCones.ispc
  ispc/geometry/ExtendedCylinders.ispc
  ispc/geometry/ExtendedSDFGeometries.ispc
  ispc/geometry/ExtendedSpheres.ispc
  ispc/render/BasicRenderer.ispc
  ispc/render/ExtendedOBJMaterial.ispc
  ispc/render/ParticleRenderer.ispc
  ispc/render/PathTracingRenderer.ispc
  ispc/render/ProximityRenderer.ispc
  ispc/render/SimulationRenderer.ispc
  ispc/render/utils/RandomGenerator.ispc
  ispc/render/utils/SkyBox.ispc
)

set(BRAYNSOSPRAYENGINEPLUGIN_SOURCES
  OSPRayCamera.cpp
  OSPRayEngine.cpp
  OSPRayFrameBuffer.cpp
  OSPRayMaterial.cpp
  OSPRayModel.cpp
  OSPRayRenderer.cpp
  OSPRayScene.cpp
  OSPRayVolume.cpp
  utils.cpp
  ispc/camera/ClippedPerspectiveCamera.cpp
  ispc/geometry/ExtendedCones.cpp
  ispc/geometry/ExtendedCylinders.cpp
  ispc/geometry/ExtendedSDFGeometries.cpp
  ispc/geometry/ExtendedSpheres.cpp
  ispc/render/BasicRenderer.cpp
  ispc/render/ExtendedOBJMaterial.cpp
  ispc/render/ParticleRenderer.cpp
  ispc/render/PathTracingRenderer.cpp
  ispc/render/ProximityRenderer.cpp
  ispc/render/SimulationRenderer.cpp
  ispc/render/utils/AbstractRenderer.cpp
)
list(APPEND BRAYNSOSPRAYENGINEPLUGIN_SOURCES ${BRAYNSOSPRAYPLUGIN_ISPC_SOURCES})

set(BRAYNSOSPRAYENGINEPLUGIN_PUBLIC_HEADERS
  OSPRayCamera.h
  OSPRayEngine.h
  OSPRayFrameBuffer.h
  OSPRayMaterial.h
  OSPRayModel.h
  OSPRayRenderer.h
  OSPRayScene.h
  OSPRayVolume.h
  ispc/camera/ClippedPerspectiveCamera.h
  ispc/geometry/ExtendedCones.h
  ispc/geometry/ExtendedCylinders.h
  ispc/geometry/ExtendedSDFGeometries.h
  ispc/geometry/ExtendedSpheres.h
  ispc/render/BasicRenderer.h
  ispc/render/ExtendedOBJMaterial.h
  ispc/render/ParticleRenderer.h
  ispc/render/PathTracingRenderer.h
  ispc/render/ProximityRenderer.h
  ispc/render/SimulationRenderer.h
  ispc/render/utils/AbstractRenderer.h
)

set(BRAYNSOSPRAYENGINEPLUGIN_HEADERS
  utils.h
)

set_source_files_properties(
    OSPRayScene.cpp
    OSPRayMaterial.cpp
    PROPERTIES COMPILE_FLAGS -Wno-deprecated-declarations)

set(BRAYNSOSPRAYENGINEPLUGIN_LINK_LIBRARIES
    PRIVATE vmmlib braynsCommon braynsParameters ${OSPRAY_LIBRARIES}
)

# reuse ispc setup and macros from ospray
list(APPEND CMAKE_MODULE_PATH ${OSPRAY_CMAKE_ROOT})
include(ispc)

# Enable/disable use of hardware randomizer
option(BRAYNS_HARDWARE_RANDOMIZER_ENABLED "Enables hardware randomizer if supported by target CPU" OFF)
if(BRAYNS_HARDWARE_RANDOMIZER_ENABLED)
  set(ISPC_EXTRA_ARGS -DBRAYNS_ISPC_USE_HARDWARE_RANDOMIZER=ON)
endif()

# Compile ispc code
include_directories_ispc(${PROJECT_SOURCE_DIR} ${OSPRAY_INCLUDE_DIRS})
ospray_ispc_compile(${BRAYNSOSPRAYPLUGIN_ISPC_SOURCES})
list(APPEND BRAYNSOSPRAYENGINEPLUGIN_SOURCES ${ISPC_OBJECTS})

# Compile C++ code
set(BRAYNSOSPRAYENGINEPLUGIN_OMIT_LIBRARY_HEADER ON)
set(BRAYNSOSPRAYENGINEPLUGIN_OMIT_VERSION_HEADERS ON)

common_library(braynsOSPRayEnginePlugin)
