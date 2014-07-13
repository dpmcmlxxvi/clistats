CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# Set the install prefix to build directory
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
    SET(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)

PROJECT(clistats)

# Add executable to project
SET(clistats_src "src/clistats.cpp")
ADD_EXECUTABLE(clistats ${clistats_src})

# Add executable to install project
INSTALL(TARGETS clistats RUNTIME DESTINATION .)

