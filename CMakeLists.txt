cmake_minimum_required(VERSION 3.22)

project(n2d VERSION 0.1.0 LANGUAGES CXX)

option(N2D_USE_AMGCL "Enable amgcl" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(assimp CONFIG REQUIRED)
find_package(CLI11 CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)

add_library(n2d
    src/io/inspect.cpp
)

target_include_directories(n2d
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_features(n2d PUBLIC cxx_std_20)

target_link_libraries(n2d
    PUBLIC
        assimp::assimp
        nlohmann_json::nlohmann_json
)

add_executable(n2d_cli
    apps/n2d_cli.cpp
)

target_link_libraries(n2d_cli
    PRIVATE
        n2d
        CLI11::CLI11
        spdlog::spdlog
)

add_subdirectory(tests)
