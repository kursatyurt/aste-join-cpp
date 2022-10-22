cmake_minimum_required (VERSION 3.1)

project(ASTE-Joiner)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED YES)
set (CMAKE_CXX_EXTENSIONS NO)

find_package(Boost 1.65.1 REQUIRED COMPONENTS log log_setup system program_options filesystem unit_test_framework)

option(ASTE_USE_VTK_COMPONENTS "Find VTK with components" ON)
if(ASTE_USE_VTK_COMPONENTS)
  # This is broken on ubuntu 20.04
  find_package(VTK REQUIRED COMPONENTS CommonCore CommonDataModel IOXML IOLegacy)
else()
  # This is broken on arch
  find_package(VTK REQUIRED)
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

add_executable(aste-join main.cpp)

message(STATUS "Boost_INCLUDE_DIRS: ${Boost_LIBRARIES}")

target_link_libraries(aste-join
  ${Boost_LIBRARIES}
  ${VTK_LIBRARIES}
)

if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
else ()
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS aste-join
    MODULES ${VTK_LIBRARIES}
    )
endif()
