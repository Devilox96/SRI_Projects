cmake_minimum_required(VERSION 3.12)
project(ShallowWater_Simplified_)

set(CMAKE_CXX_STANDARD 11)

#---Qt---#
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(QtPath /home/devilox/DeviloxDocuments/ProgTools/Qt/5.11.1/gcc_64/lib/cmake)
set(QtPackages
        Qt5::Core
        Qt5::Widgets
        Qt5::DataVisualization)

find_package(Qt5 COMPONENTS Core Widgets DataVisualization PATHS ${QtPath} NO_DEFAULT_PATH)
#---Qt---#

add_subdirectory(Libs/dMath)

add_executable(ShallowWater_Simplified_
        main.cpp
        include/mainwindow.cpp)
<<<<<<< HEAD
target_link_libraries(ShallowWater_Simplified_ ${QtPackages} ${dMath})
=======
target_link_libraries(ShallowWater_Simplified_ ${QtPackages})
>>>>>>> f10d94e86212ecdf5334690c6c08d22fb72cd80c
