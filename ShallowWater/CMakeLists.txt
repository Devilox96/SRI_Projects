cmake_minimum_required(VERSION 3.12)
project(ShallowWater)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(Libs/dMath)

#---Qt---#
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

set(QtPath /opt/Qt/5.10.1/gcc_64/lib/cmake)
set(QtPackages
        Qt5::Core
        Qt5::Widgets
        Qt5::DataVisualization)

find_package(Qt5 COMPONENTS Core Widgets DataVisualization PATHS ${QtPath} NO_DEFAULT_PATH)
#---Qt---#

add_executable(ShallowWater
        main.cpp
        Calculation/Calculation.cpp
        UI/MainWindow.cpp)

target_link_libraries(ShallowWater dMathLib ${QtPackages})