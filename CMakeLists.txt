
cmake_minimum_required(VERSION 3.0)

SET(CMAKE_VERBOSE_MAKEFILE ON)

### TODO
SET(CMAKE_BUILD_TYPE Debug)

# =========================================================================

set(VCPKG_FLAG BOOL OFF)
if (DEFINED ENV{CMAKE_TOOLCHAIN_FILE} )
    SET (CMAKE_TOOLCHAIN_FILE $ENV{CMAKE_TOOLCHAIN_FILE} )
    MESSAGE("CMAKE_TOOLCHAIN_FILE: ${CMAKE_TOOLCHAIN_FILE}")
    set(VCPKG_FLAG BOOL ON)
elseif (DEFINED ENV{VCPKG_ROOT})
    SET(VCPKG_ROOT_VAR $ENV{VCPKG_ROOT} )
    SET(CMAKE_TOOLCHAIN_FILE ${VCPKG_ROOT_VAR}/scripts/buildsystems/vcpkg.cmake )
    MESSAGE("VCPKG_ROOT=${VCPKG_ROOT_VAR}")
    MESSAGE("CMAKE_TOOLCHAIN_FILE:= ${CMAKE_TOOLCHAIN_FILE}")
    set(VCPKG_FLAG BOOL ON)
else()
    MESSAGE("VCPKG not found..")
endif()

# =========================================================================

project(NLOPTWR)

# =========================================================================

# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

# =========================================================================

if(DEFINED ENV{VCPKG_DEFAULT_TRIPLET} AND NOT DEFINED VCPKG_TARGET_TRIPLET)
  set(VCPKG_TARGET_TRIPLET "$ENV{VCPKG_DEFAULT_TRIPLET}" CACHE STRING "")
endif()

if (VCPKG_FLAG)
    MESSAGE("VCPKG found.")
    # include(${CMAKE_TOOLCHAIN_FILE})
else()
    MESSAGE("VCPKG not found.")
endif()

# =========================================================================

# Compile with C++11 standard
SET(CMAKE_CXX_STANDARD 11)

# SET(CMAKE_BUILD_TYPE Debug)
if (NOT DEFINED CMAKE_BUILD_TYPE)
# SET(CMAKE_BUILD_TYPE Release)
# SET(CMAKE_BUILD_TYPE RelWithDebInfo)
SET(CMAKE_BUILD_TYPE Debug)
endif()

# =========================================================================

FIND_PACKAGE (OpenMP  REQUIRED)
if (OPENMP_FOUND)
    SET (HAVE_OPENMP 1)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ")
endif()

# =========================================================================

if (VCPKG_FLAG)
    # use vcpkg
	message("VCPKG_FLAG=${VCPKG_FLAG}")

    find_package(NLopt CONFIG REQUIRED)
	message("NLOPT_LIBRARIES (vcpkg):      NLopt::nlopt ")
    
    # target_link_libraries(main PRIVATE NLopt::nlopt)
    SET(NLOPTCXX_LIBS NLopt::nlopt )
    
    # MESSAGE("VCPKG_INSTALLED_DIR:  ${VCPKG_INSTALLED_DIR}")
    set(VCPKG_INCL "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include")
    MESSAGE("VCPKG_INSTALLED_DIR/VCPKG_TARGET_TRIPLET/include:  ${VCPKG_INCL}")
    
else() 
    # do not use vcpkg
    
    find_package(NLOPT  REQUIRED)
    if (NLOPT_FOUND)
	    message("NLOPT_LIBRARY (no vcpkg): ${NLOPT_LIBRARY} ")
    else(NLOPT_FOUND)
        message("NLOPT not found. ")
    endif(NLOPT_FOUND)

    find_package(NLOPTCXX REQUIRED)
    if (NLOPTCXX_FOUND)
        message("  NLOPTCXX_LIBRARIES: ${NLOPTCXX_LIBRARIES} ")
    else(NLOPTCXX_FOUND)
        message("NLOPTCXX not found. ")
    endif(NLOPTCXX_FOUND)
    
    SET(NLOPTCXX_LIBS ${NLOPT_LIBRARY} ${NLOPTCXX_LIBRARY} )

endif()

# =========================================================================

# check if Doxygen is installed
find_package(Doxygen)

SET (NLOPTWR_VERSION_MAJOR "0")
SET (NLOPTWR_VERSION_MINOR "1")

# =========================================================================

include(GNUInstallDirs)

add_subdirectory(nloptwr)
add_subdirectory(nloptwrtest)

# =========================================================================

if (CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")
    # build the docs
    if (DOXYGEN_FOUND)
        message("A technical documentation will be generated with doxygen")

        message("CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
        message("DOXYGEN_EXECUTABLE:       ${DOXYGEN_EXECUTABLE}")
        message("CMAKE_CURRENT_LIST_DIR:   ${CMAKE_CURRENT_LIST_DIR}")
        
        # SET(DOXYGEN_IN  ${CMAKE_CURRENT_LIST_DIR}/Doxyfile.in)
        SET(DOXYGEN_OUT ${CMAKE_CURRENT_LIST_DIR}/Doxyfile)

        # request to configure the file
        # configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        message("Doxygen build started")

        # note the option ALL which allows to build the docs together with the application
        add_custom_target( doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} 
            WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
    else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
    endif (DOXYGEN_FOUND)
endif()

# =========================================================================

SET (NLOPTWR_VERSION_MAJOR "0")
SET (NLOPTWR_VERSION_MINOR "1")

# ----------------------------

SET (CPACK_SOURCE_IGNORE_FILES build/*,./.kdev4,.git/* )
SET(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/readme.md")

# NEU 
SET (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

IF (UNIX)
SET (CPACK_SOURCE_GENERATOR "TGZ;ZIP;DEB;RPM")
SET (CPACK_GENERATOR "TGZ;TBZ2;DEB;RPM")
ELSE (UNIX)
SET (CPACK_SOURCE_GENERATOR "ZIP")
SET (CPACK_GENERATOR "NSIS")
ENDIF (UNIX)

# build a CPack driven installer package
SET (CPACK_RESOURCE_FILE_LICENSE "${NLOPTWR_SOURCE_DIR}/License.txt")

SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Description of NLOPTWR")
SET (CPACK_PACKAGE_VENDOR "opherden")
SET (CPACK_PACKAGE_DESCRIPTION_FILE "${NLOPTWR_SOURCE_DIR}/readme.md")
SET (CPACK_PACKAGE_VENDOR " Team")
SET (CPACK_PACKAGE_CONTACT "Knut Opherden <knut.opherden@hotmail.com>")
SET (CPACK_PACKAGE_SECTION "library")
SET (CPACK_PACKAGE_VERSION_MAJOR "${NLOPTWR_VERSION_MAJOR}")
SET (CPACK_PACKAGE_VERSION_MINOR "${NLOPTWR_VERSION_MINOR}")

SET(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})

# include the necessary files:
include (CPack)

# =========================================================================

ENABLE_TESTING()

# =========================================================================

