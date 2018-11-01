cmake_minimum_required(VERSION 3.12)
project(SphereTessellation)

set(CMAKE_CXX_STANDARD 14)

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

find_package(OpenGL REQUIRED)

add_executable(SphereTessellation
        main.cpp
        MainWindow.cpp
        Include/DivisionGL.cpp)
target_link_libraries(SphereTessellation
        ${QtPackages}
        ${OPENGL_LIBRARIES})