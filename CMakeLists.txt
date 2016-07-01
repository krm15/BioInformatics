PROJECT( ACME )

# CMake >= 3.4 MUST be installed
CMAKE_MINIMUM_REQUIRED( VERSION 3.4 )

#---------------------------------------------------------------------------
# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
    "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

ADD_EXECUTABLE(readFastQ ReadFastQ.cxx)
TARGET_LINK_LIBRARIES( readFastQ ${ITK_LIBRARIES} )

