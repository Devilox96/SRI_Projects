cmake_minimum_required(VERSION 3.19)
project(MagneticWaves)

set(CMAKE_CXX_STANDARD 17)

find_package(Vulkan REQUIRED FATAL_ERROR)
find_package(SDL2 REQUIRED FATAL_ERROR)
find_package(OpenMP REQUIRED)

add_executable(MagneticWaves
        main.cpp
        lib/dMath/Core/dVectors.cpp
        lib/dMath/NumerCalc/dWENO2D.cpp
        lib/dPlot/Renderer.cpp
        lib/dPlot/dPlotMeshBase.cpp
        lib/dPlot/dVulkanMesh.cpp
        lib/dPlot/dDensityPlot2D.cpp)
target_include_directories(MagneticWaves PRIVATE Vulkan::Vulkan)
target_include_directories(MagneticWaves PRIVATE ${SDL2_INCLUDE_DIRS})
target_link_libraries(MagneticWaves
        ${QtPackages}
        Vulkan::Vulkan
        ${SDL2_LIBRARIES}
        OpenMP::OpenMP_CXX)