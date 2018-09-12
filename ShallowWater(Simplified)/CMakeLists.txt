cmake_minimum_required(VERSION 3.12)
project(ShallowWater_Simplified_)

set(CMAKE_CXX_STANDARD 11)

#---Qt---#
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(QtPath /opt/Qt/5.10.1/gcc_64/lib/cmake)
set(QtPackages
        Qt5::Core
        Qt5::Widgets
        Qt5::DataVisualization)

find_package(Qt5 COMPONENTS Core Widgets DataVisualization PATHS ${QtPath} NO_DEFAULT_PATH)
#---Qt---#

add_executable(ShallowWater_Simplified_
        main.cpp
        include/dvectors.cpp
        include/mainwindow.cpp)
target_link_libraries(ShallowWater_Simplified_ ${QtPackages})