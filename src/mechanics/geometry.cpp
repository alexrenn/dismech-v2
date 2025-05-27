#include "geometry.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <stdexcept>
#include <iostream>
#include <filesystem>
#include <unordered_set>
#include "../eigenIncludes.h"

Geometry::Geometry(const std::string& fname)
{   
    std::cout << "Reading mesh from file: " << fname << std::endl;
    from_txt(fname);
    std::cout << "Finished reading mesh.\n";

    // Separate the joint edges from rod edges
    auto [rod_shell_joint_edges, rod_edges] = separate_joint_edges(face_nodes_, edges_);
    rod_shell_joint_edges_ = rod_shell_joint_edges;
    rod_edges_ = rod_edges;

    std::cout << "Separated joint edges and rod edges.\n";
    std::cout << "Rod-shell joint edges: " << rod_shell_joint_edges_.size() << " entries, each with 2 ints\n";
    std::cout << "Rod edges: " << rod_edges_.size() << " entries, each with 2 ints\n";

    // General counting (for internal use, if needed)
    n_nodes = nodes_.size();
    n_rod_edges = rod_edges_.size();
    n_rod_shell_joints = rod_shell_joint_edges_.size();
    n_faces = face_nodes_.size();

    std::cout << "Counts:\n";
    std::cout << "  Nodes: " << n_nodes << "\n";
    std::cout << "  Rod edges: " << n_rod_edges << "\n";
    std::cout << "  Rod-shell joints: " << n_rod_shell_joints << "\n";
    std::cout << "  Faces: " << n_faces << "\n";


    // Initialize shell-related arrays
    n_edges = 3 * n_faces; 
    shell_edges_.resize(n_edges, {0, 0}); // n_edges entries, each initialized to {0, 0}
    third_node_.resize(n_edges, -1);
    face_shell_edges_.resize(n_faces, std::array<int, 3>{-1, -1, -1}); // Initialize with -1
    sign_faces_.resize(n_faces, std::array<int, 3>{0, 0, 0});          // Initialize with 0
    edge_faces_.resize(n_edges, std::array<int, 2>{-1, -1});           // Initialize with -1
    As_.resize(n_faces, 0.0); // n_faces e    edge_faces_.resize(n_edges, std::array<int, 2>{-1, -1});   
    face_unit_norms_.resize(n_faces, Eigen::Vector3d::Zero());
    face_edges_.resize(n_faces, std::array<int, 3>{-1, -1, -1}); // Initialize with -1
    hinges_.resize(n_edges); // Resize hinges_ to the maximum possible size



    std::cout << "Shell-related arrays initialized.\n";
    std::cout << "  Total edges (expected 3 * n_faces): " << n_edges << "\n";
    std::cout << "  third_node_ size: " << third_node_.size() << "\n";
    std::cout << "  As_ size: " << As_.size() << "\n";
    std::cout << "  face_unit_norms_ size: " << face_unit_norms_.size() << "\n";

    // Calculate shell and hinge edges
    calculate_shell_and_hinge_edges(n_faces);

    // Trim unused values
    trim_unused_values();

    // Ghost edges for rod-shell joint bent-twist springs
    calculate_ghost_edges(n_rod_shell_joints, n_faces);

    // Sequence edges
    sequence_edges(n_faces, face_nodes_);

    create_stretch_springs();

    calculate_face_edges(face_nodes_);

    calculate_twist_angles();
}

void Geometry::calculate_shell_and_hinge_edges(int n_faces) {
    // Iterate over faces

    std::cout << "Starting calculate_shell_and_hinge_edges with n_faces = " << n_faces << std::endl;

    for (int i = 0; i < n_faces; ++i) {
        int n1 = face_nodes_[i][0];
        int n2 = face_nodes_[i][1];
        int n3 = face_nodes_[i][2];
        
        Eigen::Vector3d p1 = nodes_[n1];
        Eigen::Vector3d p2 = nodes_[n2];
        Eigen::Vector3d p3 = nodes_[n3];
        
        // Calculate face normal and area (cross product and norm)
        Eigen::Vector3d face_cross_prod = (p2 - p1).cross(p3 - p1);
        double face_norm = face_cross_prod.norm();
        As_[i] = face_norm / 2;
        face_unit_norms_[i] = face_cross_prod / face_norm;

        std::cout << std::setprecision(17);  // match Python's float precision
        const auto& n = face_unit_norms_[i];
        std::cout << "Face " << i << " area: " << As_[i]
                << ", normal: [" << n(0) << "." << n(1) << "." << n(2) << ".]\n";

        // Iterate over edge pairs (permutations of edges)
        std::vector<std::array<int, 3>> permutations = {
            {n2, n3, n1}, {n3, n1, n2}, {n1, n2, n3}
        };

        for (int j = 0; j < permutations.size(); ++j) {
            int n1_perm = permutations[j][0];
            int n2_perm = permutations[j][1];
            int n3_perm = permutations[j][2];

            std::array<int, 2> edge = {n1_perm, n2_perm};
            std::array<int, 2> edge_neg = {n2_perm, n1_perm};

            std::cout << "Checking edge: {" << edge[0] << ", " << edge[1] << "} (negative: {"
            << edge_neg[0] << ", " << edge_neg[1] << "})" << std::endl;

            // Check if edge already exists
            bool edge_exists = false;
            int exist_id = -1;
            for (int k = 0; k < s_i_; ++k) {
                if (shell_edges_[k] == edge_neg) {
                    edge_exists = true;
                    exist_id = k;
                    break;
                }
                if (shell_edges_[k] == edge) {
                    edge_exists = true;
                    exist_id = k;
                    break;
                }
            }

            if (!edge_exists) {
                // Free-standing edge
                std::cout << "Edge {" << edge[0] << ", " << edge[1] << "} is new. Adding to shell_edges_." << std::endl;
                shell_edges_[s_i_] = edge; 
                third_node_[s_i_] = n3_perm;
                face_shell_edges_[i][j] = s_i_;
                sign_faces_[i][j] = 1;
                edge_faces_[s_i_][0]= i;
                edge_faces_[s_i_][1] = i;
                ++s_i_;
            } else {
                // Existing hinge edge
                std::cout << "Edge {" << edge[0] << ", " << edge[1] << "} already exists at index " << exist_id << "." << std::endl;
                std::array<int, 2> existing_edge = shell_edges_[exist_id];
                int existing_n3 = third_node_[exist_id];
                hinges_[h_i_] = {n1_perm, n2_perm, existing_n3, n3_perm};
                face_shell_edges_[i][j] = exist_id;

                // Determine sign based on direction of existing edge
                if (existing_edge == edge) { 
                    sign_faces_[i][j] = 1;
                } else {
                    sign_faces_[i][j] = -1;
                }

                edge_faces_[exist_id][1] = i;
                ++h_i_;
            }
        }
    }
    std::cout << "Finished calculate_shell_and_hinge_edges. Total edges: " << s_i_
    << ", Total hinges: " << h_i_ << std::endl;
}

void Geometry::trim_unused_values() {
    shell_edges_.resize(s_i_); // Keep only the first s_i_ rows
    hinges_.resize(h_i_); // Keep only the first h_i_ rows
}

void Geometry::calculate_ghost_edges(int n_rod_shell_joints, int n_faces) {
    std::cout << "Starting calculate_ghost_edges with n_rod_shell_joints = " << n_rod_shell_joints
    << ", n_faces = " << n_faces << std::endl;

    // Ghost edges for rod-shell joint bent-twist springs
    ghost_rod_shell_joint_edges_.push_back({0, 0});
    std::cout << "Initialized ghost_rod_shell_joint_edges_ with a dummy edge {0, 0}." << std::endl;


    for (int i = 0; i < n_rod_shell_joints; ++i) {
        int s_node = rod_shell_joint_edges_[i][1];
        std::vector<int> s_faces;
        std::cout << "Processing rod-shell joint " << i << " with s_node = " << s_node << std::endl;

        // Find faces with s_node
        for (int j = 0; j < n_faces; ++j) {
            if (std::find(face_nodes_[j].begin(), face_nodes_[j].end(), s_node) != face_nodes_[j].end()) {
                s_faces.push_back(j);
                std::cout << "  Found face " << j << " containing s_node = " << s_node << std::endl;
            }
        }

        // Add any edges that are not already considered
        std::vector<int> s_edges;
        for (int j = 0; j < s_faces.size(); ++j) {
            std::array<int, 3> temp_edges = face_shell_edges_[s_faces[j]];
            std::cout << "  Processing face " << s_faces[j] << " with edges: {"
                << temp_edges[0] << ", " << temp_edges[1] << ", " << temp_edges[2] << "}" << std::endl;
            for (int k = 0; k < 3; ++k) {
                std::array<int, 2> edge = shell_edges_[temp_edges[k]];
                std::cout << "    Checking edge {" << edge[0] << ", " << edge[1] << "}" << std::endl;
                if (std::find(ghost_rod_shell_joint_edges_.begin(), ghost_rod_shell_joint_edges_.end(), edge) == ghost_rod_shell_joint_edges_.end()) {
                    if (std::find(s_edges.begin(), s_edges.end(), temp_edges[k]) == s_edges.end()) {
                        s_edges.push_back(temp_edges[k]);
                        std::cout << "      Added edge index " << temp_edges[k] << " to s_edges." << std::endl;
                    }
                }
            }
        }

        // Add the new edges to the ghost edges list
        for (int s : s_edges) {
            ghost_rod_shell_joint_edges_.push_back(shell_edges_[s]);
            std::cout << "  Added ghost edge {" << shell_edges_[s][0] << ", " << shell_edges_[s][1] << "} to ghost_rod_shell_joint_edges_." << std::endl;
        }
    }
    std::cout << "Finished calculate_ghost_edges. Total ghost edges: " << ghost_rod_shell_joint_edges_.size() << std::endl;
}

void Geometry::calculate_bend_twist_springs(int n_nodes) {
    // Remove jugaad and concatenate rod_shell_joint_edges_
    std::vector<std::array<int, 2>> rod_shell_joint_edges_total = rod_shell_joint_edges_;
    // Appends a range of elements ([begin()+1, end())) from ghost_rod_shell_joint_edges_ into rod_shell_joint_edges_total.
    rod_shell_joint_edges_total.insert(
        rod_shell_joint_edges_total.end(), 
        ghost_rod_shell_joint_edges_.begin() + 1, 
        ghost_rod_shell_joint_edges_.end() 
    );

    std::vector<std::array<int, 2>> rod_edges_modified = rod_edges_;
    rod_edges_modified.insert(
        rod_edges_modified.end(), 
        rod_shell_joint_edges_total.begin(), 
        rod_shell_joint_edges_total.end()
    );

    for (int i = 0; i < n_nodes; ++i) {
        // Find edges that point into/out of the center node
        std::vector<int> into, outof;
        
        for (int j = 0; j < rod_edges_modified.size(); ++j) {
            if (rod_edges_modified[j][1] == i) {
                into.push_back(j);
            } else if (rod_edges_modified[j][0] == i) {
                outof.push_back(j);
            }
        }

        // Pairs to be created
        struct SpringPair {
            std::array<int, 2> sign1;
            std::array<int, 2> sign2;
            std::vector<std::array<int, 2>> spring_edges;
        };

        std::vector<SpringPair> pairs;

        // 1. All combinations of two edges that point into the node
        if (into.size() >= 2) {
            std::vector<std::array<int, 2>> combinations = get_combinations(into);
            pairs.push_back({{0, 0}, {1, -1}, combinations});
        }

        // 2. All combinations of two edges that point out of the node
        if (outof.size() >= 2) {
            std::vector<std::array<int, 2>> combinations = get_combinations(outof);
            pairs.push_back({{1, 1}, {-1, 1}, combinations});
        }

        // 3. All combinations of an edge into/out of the node
        if (!outof.empty() && !into.empty()) {
            std::vector<std::array<int, 2>> combinations = get_combinations_outof_into(into, outof);
            pairs.push_back({{0, 1}, {1, 1}, combinations});
        }

        // By pairing the combinations with indices, we can use one loop
        for (auto& spring_pair : pairs) {
            int n1 = spring_pair.sign1[0];
            int n2 = spring_pair.sign1[1];
            int s1 = spring_pair.sign2[0];
            int s2 = spring_pair.sign2[1];
            // Eigen::MatrixXd spring_nodes(spring_edges.size(), 3);
            for (const auto& edge_pair : spring_pair.spring_edges) {
                int node1 = rod_edges_modified[edge_pair[0]][n1];
                int node2 = i;
                int node3 = rod_edges_modified[edge_pair[1]][n2];

                bend_twist_springs_.push_back({node1, edge_pair[0], node2, edge_pair[1], node3});
                bend_twist_signs_.push_back({s1, s2});

            }    
        }
    }
}

// TODO can make cleaner
// Helper function to get combinations
std::vector<std::array<int, 2>> Geometry::get_combinations(const std::vector<int>& vec) {
    std::vector<std::array<int, 2>> combinations;
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j = i + 1; j < vec.size(); ++j) {
            combinations.push_back({vec[i], vec[j]});
        }
    }
    return combinations;
}

// Helper function to get combinations between into and outof
std::vector<std::array<int, 2>> Geometry::get_combinations_outof_into(const std::vector<int>& into, const std::vector<int>& outof) {
    std::vector<std::array<int, 2>> combinations;
    for (int i : into) {
        for (int j : outof) {
            combinations.push_back({i, j});
        }
    }
    return combinations;
}

void Geometry::sequence_edges(int n_faces, std::vector<std::array<int, 3>>& face_nodes_) {
    // Sequence edges by concatenating rod edges with rod shell joint edges
    edges_ = rod_edges_;
    edges_.insert(edges_.end(), rod_shell_joint_edges_total_.begin(), rod_shell_joint_edges_total_.end());

    // Only add unique shell_edges
    if (!edges_.empty()) {
        for (const auto& shell_edge : shell_edges_) {
            bool exists = false;
            for (const auto& edge : edges_) {
                if (edge == shell_edge) {
                    exists = true;
                    break;
                }
            }
            if (!exists) {
                edges_.push_back(shell_edge);
            } else {
                edges_ = shell_edges_;
            }
        } 
    }
    // Extract shell_edges_ slice from edges_
    int n_rod_edges = static_cast<int>(rod_edges_.size());
    int n_rod_shell_joints = static_cast<int>(rod_shell_joint_edges_total_.size());
    int start_idx = n_rod_edges + n_rod_shell_joints;
    shell_edges_.clear();
    shell_edges_.insert(shell_edges_.end(), edges_.begin() + start_idx, edges_.end());
}

void Geometry::create_stretch_springs() {
    // Clear and fill rod_stretch_springs_
    rod_stretch_springs_.clear();
    rod_stretch_springs_.insert(rod_stretch_springs_.end(), rod_edges_.begin(), rod_edges_.end());
    rod_stretch_springs_.insert(rod_stretch_springs_.end(), rod_shell_joint_edges_.begin(), rod_shell_joint_edges_.end());

    // Copy shell_edges_ directly into shell_stretch_springs_
    shell_stretch_springs_ = shell_edges_;
}

void Geometry::calculate_face_edges(std::vector<std::array<int, 3>>& face_nodes_) {
    // Calculate face edges
    for (int i = 0; i < n_faces; ++i) {
        int n1 = face_nodes_[i][0];
        int n2 = face_nodes_[i][1];
        int n3 = face_nodes_[i][2];

        std::vector<std::array<int, 2>> permutations = {{n2, n3}, {n3, n1}, {n1, n2}};

        // Loop through permutations and assign the edge indices
        for (int j = 0; j < 3; ++j) {
            const auto& edge = permutations[j];
            if (sign_faces_[i][j] > 0) {
                face_edges_[i][j] = find_edge_index(edge[0], edge[1]); // Edge index for face i, edge j
            } else {
                face_edges_[i][j]= find_edge_index(edge[1], edge[0]);
            }
        }
    }
    // Python TODO: make this mutable
}

int Geometry::find_edge_index(int n1, int n2) {
    // Find the index of an edge in the __edges_ array
    for (int i = 0; i < edges_.size(); ++i) {
        if ((edges_[i][0] == n1 && edges_[i][1] == n2) || (edges_[i][0] == n2 && edges_[i][1] == n1)) {
            return i;
        }
    }
    return -1;  // Return -1 if not found, though it should always be found
}

void Geometry::calculate_twist_angles() {
    // Placeholder method to handle twist angle calculation (if needed later)
}

void Geometry::from_txt(const std::string& fname) {
    if (!std::filesystem::exists(fname)) {
        throw std::runtime_error(fname + " is not a valid path");
    }
    
    // Constants
    enum Header { NODES = 0, EDGES = 1, TRIANGLES = 2 };
    const std::unordered_map<std::string, Header> header_map = {
        {"*nodes", NODES}, {"*edges", EDGES}, {"*triangles", TRIANGLES}
    };

    std::array<size_t, 3> expected_values = {3, 2, 3};  // per line
    std::array<bool, 3> h_flag = {false, false, false};  // seen headers

    // Flags & parameters
    int cur_h = -1;  // Tracks current header
    // std::vector<Eigen::MatrixXd> params(valid_headers.size());  // For storing matrices

    std::vector<std::vector<double>> temp;  // Temporary storage for values

    std::ifstream file(fname);
    std::string line;
    size_t line_num = 0;

    while (std::getline(file, line)) {
        ++line_num;

        // Trim
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;

        // Handle headers
        if (line[0] == '*') {
            std::string header = line;
            std::transform(header.begin(), header.end(), header.begin(), ::tolower); // Take every character in the string header, convert it to lowercase, and overwrite the string.

            auto it = header_map.find(header);
            if (it == header_map.end()) {
                throw std::runtime_error("Unknown header: " + line + " at line " + std::to_string(line_num));
            }

            int h_id = it->second;
            if (h_flag[h_id]) {
                throw std::runtime_error("Duplicate header: " + line + " at line " + std::to_string(line_num));
            }

            process_temp(cur_h, temp, nodes_, edges_, face_nodes_);
            cur_h = h_id;
            h_flag[h_id] = true;
            continue;
        }

        // Parse line of data
        std::stringstream ss(line);
        std::vector<double> values;
        std::string val;
        while (std::getline(ss, val, ',')) {
            try {
                values.push_back(std::stod(val));
            } catch (...) {
                throw std::runtime_error("Invalid number at line " + std::to_string(line_num) + ": " + val);
            }
        }

        if (cur_h < 0 || values.size() != expected_values[cur_h]) {
            throw std::runtime_error(
                "Incorrect number of values at line " + std::to_string(line_num) +
                " (got " + std::to_string(values.size()) +
                ", expected " + std::to_string(expected_values[cur_h]) + ")"
            );
        }

        temp.push_back(values);
    }

    // Process any remaining data
    process_temp(cur_h, temp, nodes_, edges_, face_nodes_);
}

void Geometry::process_temp(int header,
                            std::vector<std::vector<double>>& temp,
                            std::vector<Eigen::Vector3d>& nodes,
                            std::vector<std::array<int, 2>>& edges,
                            std::vector<std::array<int, 3>>& triangles) {
    if (header < 0) return;

    if (header == 0) {  // NODES
        for (const auto& vals : temp) {
            nodes.emplace_back(Eigen::Vector3d(vals[0], vals[1], vals[2]));
        }
    } else if (header == 1) {  // EDGES
        for (const auto& vals : temp) {
            edges.push_back({static_cast<int>(vals[0]) - 1,
                             static_cast<int>(vals[1]) - 1});
        }
    } else if (header == 2) {  // TRIANGLES
        for (const auto& vals : temp) {
            triangles.push_back({static_cast<int>(vals[0]) - 1,
                                 static_cast<int>(vals[1]) - 1,
                                 static_cast<int>(vals[2]) - 1});
        }
    }

    temp.clear();
}

std::pair<std::vector<std::array<int, 2>>, std::vector<std::array<int, 2>>>
Geometry::separate_joint_edges(const std::vector<std::array<int, 3>>& triangles,
                     const std::vector<std::array<int, 2>>& edges) {
                        
    cout << "all edges: " << std::endl;
    for (const auto& edge : edges) {
        cout << "Edge: " << edge[0] << ", " << edge[1] << std::endl;
    }
    // Return empty vectors if no edges
    if (edges.empty()) {
        return {{}, {}};
    }
    std::vector<std::array<int, 2>> rod_shell_joint_edges;
    std::vector<std::array<int, 2>> rod_edges;

    // Step 1: collect all unique shell node indices
    std::unordered_set<int> shell_nodes;
    for (const auto& tri : triangles) {
        shell_nodes.insert(tri[0]);
        shell_nodes.insert(tri[1]);
        shell_nodes.insert(tri[2]);
    }

    // Step 2: classify edges
    for (const auto& edge : edges) {
        bool n0_is_shell = shell_nodes.count(edge[0]);
        bool n1_is_shell = shell_nodes.count(edge[1]);

        if (n0_is_shell || n1_is_shell) {
            // It’s a joint edge — reverse if needed
            if (n0_is_shell && !n1_is_shell) {
                rod_shell_joint_edges.push_back({edge[1], edge[0]}); // reversed
            } else {
                rod_shell_joint_edges.push_back(edge); // original order
            }
        } else {
            rod_edges.push_back(edge); // neither node is in a triangle → rod-only edge
        }
    }

    return std::make_pair(rod_shell_joint_edges, rod_edges);
}

// FOR TESTING
// int main() {
//     Geometry geom("hex_parachute_n6.txt");
//     // Geometry geom("horizontal_rod_n21.txt");
//     // Geometry geom("mesh.txt");
//     return 0;
// }