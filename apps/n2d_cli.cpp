#include <CLI/CLI.hpp>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <fstream>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "n2d/io/inspect.hpp"

namespace {

std::string join_integers(const std::set<int>& values) {
    if (values.empty()) {
        return "(none)";
    }
    std::ostringstream stream;
    bool first = true;
    for (int value : values) {
        if (!first) {
            stream << ", ";
        }
        first = false;
        stream << value;
    }
    return stream.str();
}

std::string join_indices(const std::vector<unsigned int>& indices) {
    if (indices.empty()) {
        return "(none)";
    }
    std::ostringstream stream;
    bool first = true;
    for (unsigned int idx : indices) {
        if (!first) {
            stream << ", ";
        }
        first = false;
        stream << idx;
    }
    return stream.str();
}

std::string join_strings(const std::vector<std::string>& values) {
    if (values.empty()) {
        return "(none)";
    }
    std::ostringstream stream;
    bool first = true;
    for (const auto& value : values) {
        if (!first) {
            stream << ", ";
        }
        first = false;
        stream << value;
    }
    return stream.str();
}

} // namespace

int main(int argc, char** argv) {
    CLI::App app{"n2d (normal2disp): convert tangent-space normals to displacement maps"};
    app.set_version_flag("--version", "n2d version 0.1.0");

    std::filesystem::path inspect_mesh_path;
    std::optional<std::filesystem::path> inspect_json_path;

    auto* inspect_cmd = app.add_subcommand("inspect", "Inspect mesh materials, UV sets, and UDIM tiles");
    inspect_cmd->add_option("--mesh", inspect_mesh_path, "Mesh file to inspect")
        ->required()
        ->check(CLI::ExistingFile);
    inspect_cmd->add_option("--inspect-json", inspect_json_path, "Write inspection report to JSON file");

    CLI11_PARSE(app, argc, argv);

    if (inspect_cmd->parsed()) {
        auto inspection = n2d::io::inspect_mesh(inspect_mesh_path);
        if (!inspection.has_value()) {
            const auto& error = inspection.error();
            spdlog::error("Inspection failed for '{}': {}", inspect_mesh_path.string(), error.message);
            return 1;
        }

        const auto& result = inspection.value();
        spdlog::info("Mesh inspection report for {}", result.mesh_path.string());
        spdlog::info("  Meshes: {} ({} vertices, {} faces)", result.mesh_count, result.vertex_count, result.face_count);
        for (std::size_t mesh_index = 0; mesh_index < result.meshes.size(); ++mesh_index) {
            spdlog::info("    [{}] {}", mesh_index, result.meshes[mesh_index]);
        }

        spdlog::info("  Materials ({}):", result.materials.size());
        for (std::size_t material_index = 0; material_index < result.materials.size(); ++material_index) {
            const auto& material = result.materials[material_index];
            spdlog::info("    [{}] {} | UV Sets: {}", material_index, material.name, join_indices(material.uv_sets));
        }

        spdlog::info("  UV Sets ({}):", result.uv_sets.size());
        for (const auto& uv_set : result.uv_sets) {
            spdlog::info("    [{}] components={} | meshes={} | UDIMs={}",
                         uv_set.index,
                         uv_set.component_count,
                         join_strings(uv_set.meshes_using),
                         join_integers(uv_set.udim_tiles));
        }

        if (!result.udim_tiles.empty()) {
            spdlog::info("  UDIM tiles detected: {}", join_integers(result.udim_tiles));
        } else {
            spdlog::info("  UDIM tiles detected: (none)");
        }

        if (inspect_json_path.has_value()) {
            nlohmann::json json;
            json["mesh_path"] = result.mesh_path.generic_string();
            json["mesh_count"] = result.mesh_count;
            json["vertex_count"] = result.vertex_count;
            json["face_count"] = result.face_count;
            json["meshes"] = result.meshes;

            nlohmann::json materials_json = nlohmann::json::array();
            for (const auto& material : result.materials) {
                materials_json.push_back({{"name", material.name}, {"uv_sets", material.uv_sets}});
            }
            json["materials"] = std::move(materials_json);

            nlohmann::json uv_sets_json = nlohmann::json::array();
            for (const auto& uv_set : result.uv_sets) {
                uv_sets_json.push_back({
                    {"index", uv_set.index},
                    {"component_count", uv_set.component_count},
                    {"meshes", uv_set.meshes_using},
                    {"udim_tiles", std::vector<int>(uv_set.udim_tiles.begin(), uv_set.udim_tiles.end())},
                });
            }
            json["uv_sets"] = std::move(uv_sets_json);
            json["udim_tiles"] = std::vector<int>(result.udim_tiles.begin(), result.udim_tiles.end());

            std::ofstream output(*inspect_json_path);
            if (!output) {
                spdlog::error("Failed to open '{}' for writing", inspect_json_path->string());
                return 1;
            }

            output << json.dump(2);
            output.close();
            if (!output) {
                spdlog::error("Failed to write full inspection report to '{}'", inspect_json_path->string());
                return 1;
            }

            spdlog::info("Inspection report written to {}", inspect_json_path->string());
        }

        return 0;
    }

    spdlog::info("n2d starting up");
    return 0;
}

