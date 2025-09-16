#include "n2d/io/inspect.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <string>
#include <vector>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>

namespace n2d::io {

namespace {
constexpr unsigned int kMaxUVSets = AI_MAX_NUMBER_OF_TEXTURECOORDS;

std::string make_mesh_name(const aiMesh& mesh, std::size_t index) {
    if (mesh.mName.length > 0) {
        return std::string{mesh.mName.C_Str()};
    }
    return "Mesh_" + std::to_string(index);
}

int compute_udim_tile(const aiVector3D& uvw) {
    const int tile_u = static_cast<int>(std::floor(static_cast<double>(uvw.x)));
    const int tile_v = static_cast<int>(std::floor(static_cast<double>(uvw.y)));
    return 1001 + tile_u + (10 * tile_v);
}

} // namespace

tl::expected<InspectionResult, N2DError> inspect_mesh(const std::filesystem::path& mesh_path) {
    Assimp::Importer importer;
    const aiScene* scene = importer.ReadFile(mesh_path.string(), 0u);
    if (scene == nullptr) {
        return tl::unexpected(make_io_error(importer.GetErrorString() ? importer.GetErrorString() : "Failed to load mesh"));
    }

    InspectionResult result;
    result.mesh_path = mesh_path;
    result.mesh_count = scene->mNumMeshes;

    std::array<UVSetInfo, kMaxUVSets> uv_data{};
    std::array<bool, kMaxUVSets> uv_present{};
    std::array<std::set<std::string>, kMaxUVSets> uv_mesh_usage;
    for (unsigned int uv_index = 0; uv_index < kMaxUVSets; ++uv_index) {
        uv_data[uv_index].index = uv_index;
    }

    std::vector<std::set<unsigned int>> material_uv_sets(scene->mNumMaterials);

    for (unsigned int mesh_index = 0; mesh_index < scene->mNumMeshes; ++mesh_index) {
        const aiMesh* mesh = scene->mMeshes[mesh_index];
        if (mesh == nullptr) {
            continue;
        }

        result.vertex_count += mesh->mNumVertices;
        result.face_count += mesh->mNumFaces;

        std::string mesh_name = make_mesh_name(*mesh, mesh_index);
        result.meshes.push_back(mesh_name);

        for (unsigned int uv_index = 0; uv_index < kMaxUVSets; ++uv_index) {
            if (!mesh->HasTextureCoords(uv_index) || mesh->mTextureCoords[uv_index] == nullptr) {
                continue;
            }

            uv_present[uv_index] = true;
            uv_data[uv_index].component_count = std::max(uv_data[uv_index].component_count, mesh->mNumUVComponents[uv_index]);
            uv_mesh_usage[uv_index].insert(mesh_name);

            for (unsigned int vertex_index = 0; vertex_index < mesh->mNumVertices; ++vertex_index) {
                const aiVector3D& uvw = mesh->mTextureCoords[uv_index][vertex_index];
                const int udim_tile = compute_udim_tile(uvw);
                uv_data[uv_index].udim_tiles.insert(udim_tile);
                result.udim_tiles.insert(udim_tile);
            }

            if (mesh->mMaterialIndex < material_uv_sets.size()) {
                material_uv_sets[mesh->mMaterialIndex].insert(uv_index);
            }
        }
    }

    result.materials.reserve(scene->mNumMaterials);
    for (unsigned int material_index = 0; material_index < scene->mNumMaterials; ++material_index) {
        const aiMaterial* material = scene->mMaterials[material_index];
        MaterialInfo material_info{};
        if (material != nullptr) {
            aiString name;
            if (material->Get(AI_MATKEY_NAME, name) == AI_SUCCESS && name.length > 0) {
                material_info.name = name.C_Str();
            }
        }

        if (material_info.name.empty()) {
            material_info.name = "Material_" + std::to_string(material_index);
        }

        const auto& uv_sets = material_uv_sets[material_index];
        material_info.uv_sets.assign(uv_sets.begin(), uv_sets.end());
        result.materials.push_back(std::move(material_info));
    }

    for (unsigned int uv_index = 0; uv_index < kMaxUVSets; ++uv_index) {
        if (!uv_present[uv_index]) {
            continue;
        }

        UVSetInfo info = uv_data[uv_index];
        const auto& meshes_using = uv_mesh_usage[uv_index];
        info.meshes_using.assign(meshes_using.begin(), meshes_using.end());
        result.uv_sets.push_back(std::move(info));
    }

    return result;
}

} // namespace n2d::io

