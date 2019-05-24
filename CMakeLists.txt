cmake_minimum_required (VERSION 3.10)
project(LibMolGrid LANGUAGES C CXX )

#version number
set (VERSION_MAJOR 0)
set (VERSION_MINOR 1)
set (VERSION_PATCH 0)

#set default build type
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

# guard against in-source builds and bad build-type strings
include(cmake/safeguards.cmake)

# place binaries and libraries according to GNU standards
include(GNUInstallDirs)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

#dependencies
find_package(Boost REQUIRED COMPONENTS system filesystem iostreams)
find_package(OpenBabel2 REQUIRED)
find_package(RDKit)
find_package( ANN REQUIRED )
find_package( Eigen3 REQUIRED )
find_package( LEMON REQUIRED )

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_STANDARD 14)

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wall -Wno-unknown-pragmas -Werror")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_FLAGS_PROFILE " -fprofile-arcs -ftest-coverage -pg")
endif()

# so that we will find libmolgrid_config.h
include_directories(${CMAKE_SOURCE_DIR}/)
include_directories(${CMAKE_SOURCE_DIR}/shapedb)
include_directories(SYSTEM ${OPENBABEL2_INCLUDE_DIR})
include_directories(SYSTEM ${Boost_INCLUDE_DIR})
include_directories(SYSTEM ${RDKIT_INCLUDE_DIR})
include_directories( ${EIGEN3_INCLUDE_DIR} )
include_directories( ${LEMON_INCLUDE_DIRS} )

set(LEADIT_SRCS
CommandLine2/CommandLine.cpp
CommandLine2/CommandLine.h
DatabaseCreator.cpp
DatabaseCreator.h
DatabaseSearcher.cpp
DatabaseSearcher.h
DatabaseStrutures.h
FragmentIndexer.cpp
FragmentIndexer.h
FragmentSearcher.cpp
FragmentSearcher.h
main.cpp
MolMatcher.cpp
MolMatcher.h
Orienter.h
pharmacophores.cpp
pharmacophores.h
QueryObject.cpp
QueryObject.h
Reaction.cpp
Reaction.h
ScaffoldIndexer.cpp
ScaffoldIndexer.h
shapedb/Cube.h
shapedb/DataViewers.h
shapedb/GSSTreeCreator.cpp
shapedb/GSSTreeCreator.h
shapedb/GSSTreeSearcher.cpp
shapedb/GSSTreeSearcher.h
shapedb/GSSTreeStructures.cpp
shapedb/GSSTreeStructures.h
shapedb/GSSTypes.h
shapedb/KSamplePartitioner.cpp
shapedb/KSamplePartitioner.h
#shapedb/main.cpp
shapedb/MappableOctTree.cpp
shapedb/MappableOctTree.h
shapedb/MemMapped.cpp
shapedb/MemMapped.h
shapedb/MGrid.cpp
shapedb/MGrid.h
shapedb/MiraObject.h
shapedb/molecules/MolSphere.cpp
shapedb/molecules/OBMoleculeAnalytic.cpp
shapedb/molecules/PMol.cpp
shapedb/molecules/RDMoleculeAnalytic.cpp
shapedb/packers/FullMergePacker.cpp
shapedb/packers/GreedyPacker.cpp
shapedb/packers/MatcherPacker.cpp
shapedb/packers/Packer.cpp
shapedb/packers/SpectralPacker.cpp
shapedb/Results.h
shapedb/ShapeDistance.cpp
shapedb/ShapeDistance.h
shapedb/Timer.h
shapedb/TopDownPartitioner.h
shapedb/WorkFile.cpp
shapedb/WorkFile.h
)

add_executable(leadit ${LEADIT_SRCS})
target_link_libraries(leadit ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES} ${RDKIT_LIBRARIES} ${ANN_LIBRARY})

