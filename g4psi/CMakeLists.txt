# modified for g4PSI from Geant4 example B1
#
# to create the make files do
#
# cmake -DGeant4_DIR=/Applications/Research/geant4.10.01.p02-install/lib/Geant4-10.1.2 <path-to-file>

#

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(g4PSI)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build g4PSI with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
find_package(ROOT REQUIRED)
include_directories(${PROJECT_SOURCE_DIR}/include
	${ROOT_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Find ROOT (required package)
#

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add ROOT libraries properties
#
set(libsuffix .so)
set(ROOT_LIBRARY_PROPERTIES ${ROOT_LIBRARY_PROPERTIES} SUFFIX ${libsuffix})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 and ROOT libraries
#
add_executable(g4PSI g4PSI.cc ${sources} ${headers})
add_definitions(-std=c++11) 
target_link_libraries(g4PSI ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build g4PSI. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

configure_file(
    ${PROJECT_SOURCE_DIR}/mac/Settings.mac ${PROJECT_BINARY_DIR}/Settings.mac
    COPYONLY
)
configure_file(
    ${PROJECT_SOURCE_DIR}/mac/vis_ss.mac ${PROJECT_BINARY_DIR}/vis.mac
    COPYONLY
)
configure_file(
    ${PROJECT_SOURCE_DIR}/mac/gui_ss.mac ${PROJECT_BINARY_DIR}/gui.mac
    COPYONLY
)
configure_file(
    ${PROJECT_SOURCE_DIR}/icons.mac ${PROJECT_BINARY_DIR}/icons.mac
    COPYONLY
)
configure_file(
    ${PROJECT_SOURCE_DIR}/run.png ${PROJECT_BINARY_DIR}/run.png
    COPYONLY
)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
# install(TARGETS g4PSI DESTINATION bin)


