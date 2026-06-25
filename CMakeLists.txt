cmake_minimum_required(VERSION 3.21)

project(AsFem
        VERSION 0.9.0
        DESCRIPTION "Advanced Simulation kit based on Finite Element Method"
        HOMEPAGE_URL "https://github.com/MatMechLab/AsFem"
        LANGUAGES CXX)

#################################################################
### C++ standard                                              ###
#################################################################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# install layout (CMAKE_INSTALL_BINDIR etc.); needed for conda / "make install"
include(GNUInstallDirs)

#################################################################
### Build type: default to Release for single-config gens     ###
### (multi-config gens such as MSVC/Xcode pick at build time). ###
#################################################################
get_property(_is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT _is_multi_config AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()
# Backward compatibility: the old custom "Depdebug" type == Debug + AddressSanitizer.
if(CMAKE_BUILD_TYPE STREQUAL "Depdebug")
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type" FORCE)
    set(ASFEM_ENABLE_ASAN ON CACHE BOOL "" FORCE)
endif()

#################################################################
### Options (all OFF by default => portable / conda-safe)     ###
#################################################################
# -march=native bakes in the build machine's CPU features. It is great for a
# local build but produces NON-portable binaries, so it MUST stay OFF for conda
# and any redistributable package. Local users: add -DASFEM_ENABLE_NATIVE=ON.
option(ASFEM_ENABLE_NATIVE "Optimize for the build machine's CPU (-march=native); NOT portable" OFF)
option(ASFEM_ENABLE_WERROR "Treat warnings as errors (developer builds only)"                  OFF)
option(ASFEM_ENABLE_ASAN   "Build with AddressSanitizer (Debug only)"                          OFF)

#################################################################
### Dependency: MPI                                           ###
### find_package(MPI) is cross-platform: it handles a plain   ###
### compiler + system/conda MPI (MPICH/OpenMPI), MS-MPI on    ###
### Windows, or an MPI compiler wrapper (mpicxx) transparently###
#################################################################
# honour a user-provided MPI_DIR as a search hint
if(DEFINED ENV{MPI_DIR} AND NOT MPI_HOME)
    set(MPI_HOME "$ENV{MPI_DIR}")
endif()
find_package(MPI REQUIRED COMPONENTS CXX)
message(STATUS "MPI: ${MPI_CXX_COMPILER} (MPI_HOME=${MPI_HOME})")

#################################################################
### Dependency: PETSc (+ the METIS it ships with)             ###
### Primary: pkg-config (the portable & conda-native way).    ###
### Fallback: locate libpetsc/petsc.h directly, honouring     ###
### PETSC_DIR / CONDA_PREFIX. find_library/find_path resolve  ###
### the right extension (.so/.dylib/.lib) on every platform.  ###
#################################################################
# make a non-conda PETSc (pointed to by PETSC_DIR) visible to pkg-config too
if(DEFINED ENV{PETSC_DIR})
    if(DEFINED ENV{PETSC_ARCH} AND NOT "$ENV{PETSC_ARCH}" STREQUAL "")
        set(ENV{PKG_CONFIG_PATH} "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig:$ENV{PETSC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    else()
        set(ENV{PKG_CONFIG_PATH} "$ENV{PETSC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    endif()
endif()

find_package(PkgConfig)
if(PkgConfig_FOUND)
    pkg_check_modules(PETSC IMPORTED_TARGET PETSc)
endif()

if(PETSC_FOUND)
    set(ASFEM_PETSC_TARGET PkgConfig::PETSC)
    message(STATUS "PETSc: ${PETSC_VERSION} (pkg-config)")
else()
    message(STATUS "PETSc pkg-config not found, falling back to PETSC_DIR/library search")
    find_path(PETSC_INCLUDE_DIR NAMES petsc.h
              HINTS $ENV{PETSC_DIR}/include $ENV{CONDA_PREFIX}/include
              PATHS /usr/include/petsc /usr/lib/petsc/include /opt/petsc/include)
    find_library(PETSC_LIBRARY NAMES petsc
                 HINTS $ENV{PETSC_DIR}/lib $ENV{PETSC_DIR}/lib64 $ENV{CONDA_PREFIX}/lib
                 PATHS /usr/lib /usr/lib64)
    if(NOT PETSC_INCLUDE_DIR OR NOT PETSC_LIBRARY)
        message(FATAL_ERROR
            "PETSc not found.\n"
            "  - set PETSC_DIR to your PETSc install, or\n"
            "  - install petsc so that 'pkg-config --exists PETSc' succeeds (conda: 'conda install petsc').")
    endif()
    add_library(asfem_petsc INTERFACE)
    target_include_directories(asfem_petsc INTERFACE ${PETSC_INCLUDE_DIR})
    target_link_libraries(asfem_petsc INTERFACE ${PETSC_LIBRARY})
    set(ASFEM_PETSC_TARGET asfem_petsc)
    message(STATUS "PETSc: ${PETSC_LIBRARY}")
endif()

# METIS ships with PETSc (PETSc builds it); conda provides a 'metis' package.
find_library(METIS_LIBRARY NAMES metis
             HINTS $ENV{PETSC_DIR}/lib $ENV{PETSC_DIR}/lib64 $ENV{CONDA_PREFIX}/lib
             PATHS /usr/lib /usr/lib64)
if(NOT METIS_LIBRARY)
    message(FATAL_ERROR "METIS library not found (it normally ships with PETSc). "
                        "Set PETSC_DIR or install metis.")
endif()
message(STATUS "METIS: ${METIS_LIBRARY}")

#################################################################
### Sources (explicit list in SourceCode.cmake)              ###
#################################################################
include(SourceCode.cmake)

add_executable(asfem ${inc} ${src})

target_include_directories(asfem PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/external"        # nlohmann/json
    "${CMAKE_CURRENT_SOURCE_DIR}/external/eigen") # Eigen (header-only)

target_link_libraries(asfem PRIVATE
    MPI::MPI_CXX
    ${ASFEM_PETSC_TARGET}
    ${METIS_LIBRARY})

#################################################################
### Compiler flags (per-compiler, portable by default)        ###
#################################################################
if(MSVC)
    # /O2 etc. for the optimized configs are added by CMake; keep warnings sane.
    target_compile_options(asfem PRIVATE
        /W1
        /bigobj
        $<$<BOOL:${ASFEM_ENABLE_WERROR}>:/WX>)
else()
    # GCC / Clang / AppleClang (Linux + macOS)
    target_compile_options(asfem PRIVATE
        -Wall
        $<$<CONFIG:Debug>:-Wextra>
        $<$<BOOL:${ASFEM_ENABLE_NATIVE}>:-march=native -mtune=native>
        $<$<BOOL:${ASFEM_ENABLE_WERROR}>:-Werror>
        $<$<AND:$<CONFIG:Debug>,$<BOOL:${ASFEM_ENABLE_ASAN}>>:-fsanitize=address;-fno-omit-frame-pointer;-g>)
    if(ASFEM_ENABLE_ASAN)
        target_link_options(asfem PRIVATE $<$<CONFIG:Debug>:-fsanitize=address>)
    endif()
endif()
# CMake's per-config defaults already provide -O3/-O2/-g and -DNDEBUG (Release).

#################################################################
### Link-time optimization for the optimized configurations   ###
#################################################################
include(CheckIPOSupported)
check_ipo_supported(RESULT _ipo_ok OUTPUT _ipo_err)
if(_ipo_ok)
    set_property(TARGET asfem PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    set_property(TARGET asfem PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE)
    message(STATUS "IPO/LTO enabled for optimized builds")
else()
    message(STATUS "IPO/LTO not supported: ${_ipo_err}")
endif()

#################################################################
### Output location + install rules                           ###
#################################################################
# Convenience copy in <source>/bin for the in-tree run workflow (asfem -i ...).
set_target_properties(asfem PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

# Installable target (conda / package managers / "cmake --install").
install(TARGETS asfem RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

message(STATUS "AsFem ${PROJECT_VERSION} | build type: ${CMAKE_BUILD_TYPE} | native opt: ${ASFEM_ENABLE_NATIVE}")
message(STATUS "  binary -> ${CMAKE_CURRENT_SOURCE_DIR}/bin/asfem  (install -> ${CMAKE_INSTALL_FULL_BINDIR})")

#################################################################
### Tests                                                     ###
#################################################################
include(Test.cmake)
