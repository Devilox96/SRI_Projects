cmake_minimum_required(VERSION 3.10)
project(ShallowWater)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(Libs/dMath)

#---Qt---#
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

set(QtPath /home/devilox/DeviloxDocuments/ProgTools/Qt/5.11.1/gcc_64/lib/cmake)

find_package(Qt5 COMPONENTS Core Widgets DataVisualization PATHS ${QtPath})

set(QtPackages
        Qt5::Core
        Qt5::Widgets
        Qt5::DataVisualization)
#---Qt---#

add_executable(ShallowWater
        main.cpp
        Calculation/Calculation.cpp
        UI/MainWindow.cpp)

target_link_libraries(ShallowWater dMathLib ${QtPackages})