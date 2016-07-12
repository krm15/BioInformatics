PROJECT( BioInformatics )

# CMake >= 3.4 MUST be installed
CMAKE_MINIMUM_REQUIRED( VERSION 3.4 )

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY
 ${BioInformatics_BINARY_DIR}/bin
 CACHE PATH "Single output directory for building all executables."
)

SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY
 ${BioInformatics_BINARY_DIR}/lib
 CACHE PATH "Single output directory for building all libraries."
)

ADD_EXECUTABLE(readFastQ ReadFastQ.cxx)
TARGET_LINK_LIBRARIES( readFastQ )

