#pragma once

#include <cstddef>
#include <filesystem>
#include <set>
#include <string>
#include <vector>

#include <tl/expected.hpp>

#include "n2d/core/error.hpp"

namespace n2d::io {

struct MaterialInfo {
    std::string name;
    std::vector<unsigned int> uv_sets;
};

struct UVSetInfo {
    unsigned int index = 0;
    unsigned int component_count = 0;
    std::set<int> udim_tiles;
    std::vector<std::string> meshes_using;
};

struct InspectionResult {
    std::filesystem::path mesh_path;
    std::size_t mesh_count = 0;
    std::size_t vertex_count = 0;
    std::size_t face_count = 0;
    std::vector<std::string> meshes;
    std::vector<MaterialInfo> materials;
    std::vector<UVSetInfo> uv_sets;
    std::set<int> udim_tiles;
};

tl::expected<InspectionResult, N2DError> inspect_mesh(const std::filesystem::path& mesh_path);

} // namespace n2d::io

