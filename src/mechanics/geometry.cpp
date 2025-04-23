#include "geometry.h"

Geometry(const std::vector<Eigen::MatrixXd>& params) {
    // Initialize your internal data members with params
    nodes_ = params[0];
    edges_ = params[1];
    face_nodes_ = params[2]; 

    // Separate the joint edges from rod edges
    auto [rod_shell_joint_edges, rod_edges] = separate_joint_edges(face_nodes, edges);
    rod_shell_joint_edges_ = rod_shell_joint_edges;
    rod_edges_ = rod_edges;

    // General counting (for internal use, if needed)
    n_nodes = nodes_.rows();
    n_rod_edges = rod_edges_.rows();
    n_rod_shell_joints = rod_shell_joint_edges_.rows();
    n_faces = face_nodes_.rows();

    // Initialize shell-related arrays
    n_edges = 3 * n_faces;
    Eigen::MatrixXi shell_edges(n_edges, 2);
    Eigen::MatrixXi hinges(n_edges, 4);
    Eigen::VectorXi third_node(n_edges);

    Eigen::MatrixXi edge_faces(n_edges, 2);
    Eigen::VectorXd As(n_faces);

    face_shell_edges_ = Eigen::MatrixXi(n_faces, 3);
    sign_faces_ = Eigen::MatrixXi(n_faces, 3);
    face_unit_norms_ = Eigen::MatrixXd(3, n_faces);

    // Calculate shell and hinge edges
    calculate_shell_and_hinge_edges(n_faces);

    // Trim unused values
    trim_unused_values();

    // Ghost edges for rod-shell joint bent-twist springs
    calculate_ghost_edges(n_rod_shell_joints, n_faces);

    // Sequence edges
    sequence_edges(n_faces, face_nodes);

    create_stretch_springs();

    calculate_face_edges(face_nodes);

    find_edge_index(face_nodes);
    calculate_twist_angles();
}

Geometry::Geometry(const Eigen::MatrixXd& nodes,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& face_nodes)
    : nodes_(nodes), face_nodes_(face_nodes)
{
    // Separate the joint edges from rod edges
    auto [rod_shell_joint_edges, rod_edges] = separate_joint_edges(face_nodes, edges);
    rod_shell_joint_edges_ = rod_shell_joint_edges;
    rod_edges_ = rod_edges;

    // General counting (for internal use, if needed)
    n_nodes = nodes_.rows();
    n_rod_edges = rod_edges_.rows();
    n_rod_shell_joints = rod_shell_joint_edges_.rows();
    n_faces = face_nodes_.rows();

    // Initialize shell-related arrays
    n_edges = 3 * n_faces;
    Eigen::MatrixXi shell_edges(n_edges, 2);
    Eigen::MatrixXi hinges(n_edges, 4);
    Eigen::VectorXi third_node(n_edges);

    Eigen::MatrixXi edge_faces(n_edges, 2);
    Eigen::VectorXd As(n_faces);

    face_shell_edges_ = Eigen::MatrixXi(n_faces, 3);
    sign_faces_ = Eigen::MatrixXi(n_faces, 3);
    face_unit_norms_ = Eigen::MatrixXd(3, n_faces);

    // Calculate shell and hinge edges
    calculate_shell_and_hinge_edges(n_faces);

    // Trim unused values
    trim_unused_values();

    // Ghost edges for rod-shell joint bent-twist springs
    calculate_ghost_edges(n_rod_shell_joints, n_faces);

    // Sequence edges
    sequence_edges(n_faces, face_nodes);

    create_stretch_springs();

    calculate_face_edges(face_nodes);

    find_edge_index(face_nodes);
    calculate_twist_angles();
}

void Geometry::calculate_shell_and_hinge_edges(int n_faces) {
    Eigen::VectorXd As(n_faces);
    Eigen::MatrixXd face_cross(3, n_faces);
    
    // Iterate over faces
    for (int i = 0; i < n_faces; ++i) {
        int n1 = face_nodes_(i, 0);
        int n2 = face_nodes_(i, 1);
        int n3 = face_nodes_(i, 2);
        
        Eigen::VectorXd p1 = nodes_.row(n1);
        Eigen::VectorXd p2 = nodes_.row(n2);
        Eigen::VectorXd p3 = nodes_.row(n3);
        
        // Calculate face normal and area (cross product and norm)
        Eigen::VectorXd face_cross_prod = (p2 - p1).cross(p3 - p1);
        double face_norm = face_cross_prod.norm();
        As(i) = face_norm / 2;
        face_unit_norms_.col(i) = face_cross_prod / face_norm;

        // Iterate over edge pairs (permutations of edges)
        std::vector<std::tuple<int, int, int>> permutations = {
            {n2, n3, n1}, {n3, n1, n2}, {n1, n2, n3}
        };

        for (int j = 0; j < permutations.size(); ++j) {
            int n1_perm = std::get<0>(permutations[j]);
            int n2_perm = std::get<1>(permutations[j]);
            int n3_perm = std::get<2>(permutations[j]);

            Eigen::Vector2i edge(n1_perm, n2_perm);
            Eigen::Vector2i edge_neg(n2_perm, n1_perm);

            // Check if edge already exists
            bool edge_exists = false;
            int exist_id = -1;
            for (int k = 0; k < s_i_; ++k) {
                if ((shell_edges_.row(k) == edge_neg).all()) {
                    edge_exists = true;
                    exist_id = k;
                    break;
                }
                if ((shell_edges_.row(k) == edge).all()) {
                    edge_exists = true;
                    exist_id = k;
                    break;
                }
            }

            if (!edge_exists) {
                // Free-standing edge
                shell_edges_.row(s_i_) = edge;
                third_node_(s_i_) = n3_perm;
                face_shell_edges_(i, j) = s_i_;
                sign_faces_(i, j) = 1;
                edge_faces_(s_i_, 0) = i;
                edge_faces_(s_i_, 1) = i;
                ++s_i_;
            } else {
                // Existing hinge edge
                Eigen::Vector2i existing_edge = shell_edges_.row(exist_id);
                int existing_n3 = third_node_(exist_id);
                hinges_.row(h_i_) << n1_perm, n2_perm, existing_n3, n3_perm;
                face_shell_edges_(i, j) = exist_id;

                // Determine sign based on direction of existing edge
                if (existing_edge == edge) {
                    sign_faces_(i, j) = 1;
                } else {
                    sign_faces_(i, j) = -1;
                }

                edge_faces_(exist_id, 1) = i;
                ++h_i_;
            }
        }
    }
}

void Geometry::trim_unused_values() {
    shell_edges.conservativeResize(s_i, 2);
    hinges.conservativeResize(h_i, 4);
}

void Geometry::calculate_ghost_edges(int n_rod_shell_joints, int n_faces) {
    // Ghost edges for rod-shell joint bent-twist springs
    std::vector<Eigen::Vector2i> ghost_rod_shell_joint_edges;
    ghost_rod_shell_joint_edges.push_back(Eigen::Vector2i(0, 0));  // Using Eigen for vector

    for (int i = 0; i < n_rod_shell_joints; ++i) {
        int s_node = rod_shell_joint_edges_(i, 1);
        std::vector<int> s_faces;

        // Find faces with s_node
        for (int j = 0; j < n_faces; ++j) {
            if (std::find(face_nodes_.row(j).data(), face_nodes_.row(j).data() + 3, s_node) != face_nodes_.row(j).data() + 3) {
                s_faces.push_back(j);
            }
        }

        // Add any edges that are not already considered
        std::vector<int> s_edges;
        for (int j = 0; j < s_faces.size(); ++j) {
            Eigen::Vector3i temp_edges = face_shell_edges_.row(s_faces[j]);
            for (int k = 0; k < 3; ++k) {
                Eigen::Vector2i edge = shell_edges_.row(temp_edges(k));
                if (std::find(ghost_rod_shell_joint_edges.begin(), ghost_rod_shell_joint_edges.end(), edge) == ghost_rod_shell_joint_edges.end()) {
                    if (std::find(s_edges.begin(), s_edges.end(), temp_edges(k)) == s_edges.end()) {
                        s_edges.push_back(temp_edges(k));
                    }
                }
            }
        }

        // Add the new edges to the ghost edges list
        for (int s : s_edges) {
            ghost_rod_shell_joint_edges.push_back(shell_edges_.row(s));
        }
    }
}

void Geometry::calculate_bend_twist_springs(int n_nodes) {
    // Remove jugaad and concatenate rod_shell_joint_edges
    std::vector<Eigen::Vector2i> rod_shell_joint_edges_total = rod_shell_joint_edges_;
    rod_shell_joint_edges_total.insert(
        rod_shell_joint_edges_total.end(), ghost_rod_shell_joint_edges_.begin() + 1, ghost_rod_shell_joint_edges_.end()
    );

    std::vector<Eigen::Vector2i> rod_edges_modified = safe_concat(rod_edges_, rod_shell_joint_edges_total);

    // Bend-twist springs
    std::vector<Eigen::MatrixXd> bend_twist_springs;
    std::vector<Eigen::MatrixXd> bend_twist_signs;

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
        std::vector<std::tuple<std::pair<int, int>, std::pair<int, int>, std::vector<std::pair<int, int>>>> pairs;

        // 1. All combinations of two edges that point into the node
        if (into.size() >= 2) {
            std::vector<std::pair<int, int>> combinations = get_combinations(into);
            pairs.push_back({{0, 0}, {1, -1}, combinations});
        }

        // 2. All combinations of two edges that point out of the node
        if (outof.size() >= 2) {
            std::vector<std::pair<int, int>> combinations = get_combinations(outof);
            pairs.push_back({{1, 1}, {-1, 1}, combinations});
        }

        // 3. All combinations of an edge into/out of the node
        if (!outof.empty() && !into.empty()) {
            std::vector<std::pair<int, int>> combinations = get_combinations_outof_into(into, outof);
            pairs.push_back({{0, 1}, {1, 1}, combinations});
        }

        // By pairing the combinations with indices, we can use one loop
        for (auto& [sign1_sign2, sign2_sign2, spring_edges] : pairs) {
            Eigen::MatrixXd spring_nodes(spring_edges.size(), 3);

            for (int idx = 0; idx < spring_edges.size(); ++idx) {
                spring_nodes(idx, 0) = rod_edges_modified[spring_edges[idx].first][0];
                spring_nodes(idx, 1) = i;
                spring_nodes(idx, 2) = rod_edges_modified[spring_edges[idx].second][1];
            }

            bend_twist_springs.push_back(spring_nodes);

            Eigen::MatrixXd spring_signs(spring_edges.size(), 2);
            for (int idx = 0; idx < spring_edges.size(); ++idx) {
                spring_signs(idx, 0) = sign1_sign2.first;
                spring_signs(idx, 1) = sign2_sign2.second;
            }

            bend_twist_signs.push_back(spring_signs);
        }
    }

    // Concatenate all bend twist springs
    if (!bend_twist_springs.empty()) {
        this->bend_twist_springs_ = concatenate(bend_twist_springs);
        this->bend_twist_signs_ = concatenate(bend_twist_signs);
    } else {
        this->bend_twist_springs_ = Eigen::MatrixXd(0, 0);
        this->bend_twist_signs_ = Eigen::MatrixXd(0, 0);
    }
}

// Helper function to get combinations
std::vector<std::pair<int, int>> Geometry::get_combinations(const std::vector<int>& vec) {
    std::vector<std::pair<int, int>> combinations;
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j = i + 1; j < vec.size(); ++j) {
            combinations.push_back({vec[i], vec[j]});
        }
    }
    return combinations;
}

// Helper function to get combinations between into and outof
std::vector<std::pair<int, int>> Geometry::get_combinations_outof_into(const std::vector<int>& into, const std::vector<int>& outof) {
    std::vector<std::pair<int, int>> combinations;
    for (int i : into) {
        for (int j : outof) {
            combinations.push_back({i, j});
        }
    }
    return combinations;
}

void Geometry::sequence_edges(int n_faces, const Eigen::MatrixXi& face_nodes) {
    // Sequence edges by concatenating rod edges with rod shell joint edges
    this->edges_ = safe_concat(this->rod_edges_, this->rod_shell_joint_edges_total_);

    // Only add unique shell_edges
    if (this->edges_.rows() > 0) {
        for (int i = 0; i < shell_edges_.rows(); ++i) {
            Eigen::RowVector2i edge = shell_edges_.row(i);
            bool exists = false;
            for (int j = 0; j < this->edges_.rows(); ++j) {
                if ((this->edges_.row(j) == edge).all()) {
                    exists = true;
                    break;
                }
            }
            if (!exists) {edges_
                this->edges_.conservativeResize(this->edges_.rows() + 1, Eigen::NoChange);
                this->edges_.row(this->edges_.rows() - 1) = edge;
            }
        }
    } else {
        this->edges_ = shell_edges_;
    }

    this->shell_edges_ = this->edges_.block(n_rod_edges + n_rod_shell_joints, 0, this->edges_.rows() - (n_rod_edges + n_rod_shell_joints), 2);
}

void Geometry::create_stretch_springs() {
    // Create rod and shell stretch springs by concatenating their respective edges
    rod_stretch_springs_ = safe_concat(rod_edges_, rod_shell_joint_edges_);
    shell_stretch_springs_ = shell_edges_;
}

void Geometry::calculate_face_edges(const Eigen::MatrixXi& face_nodes) {
    // Initialize face edges
    __face_edges_ = Eigen::MatrixXi(n_faces_, 3);

    // Calculate face edges
    for (int i = 0; i < n_faces_; ++i) {
        int n1 = face_nodes(i, 0);
        int n2 = face_nodes(i, 1);
        int n3 = face_nodes(i, 2);

        std::array<std::pair<int, int>, 3> permutations = {{{n2, n3}, {n3, n1}, {n1, n2}}};

        // Loop through permutations and assign the edge indices
        for (int j = 0; j < 3; ++j) {
            const auto& edge = permutations[j];
            if (__sign_faces_[i, j] > 0) {
                __face_edges_(i, j) = find_edge_index(edge.first, edge.second);
            } else {
                __face_edges_(i, j) = find_edge_index(edge.second, edge.first);
            }
        }
    }
}

int Geometry::find_edge_index(int n1, int n2) {
    // Find the index of an edge in the __edges_ array
    for (int i = 0; i < edges_.size(); ++i) {
        if ((edges_  == n1 && edges_  == n2) || (edges_  == n2 && edges_  == n1)) {
            return i;
        }
    }
    return -1;  // Return -1 if not found, though it should always be found
}

void Geometry::calculate_twist_angles() {
    // Placeholder method to handle twist angle calculation (if needed later)
    twist_angles_ = Eigen::MatrixXd::Zero(n_rod_edges + this->rod_shell_joint_edges_total_.rows());
}

static Eigen::MatrixXd safe_concat(const Eigen::MatrixXd& arr1, const Eigen::MatrixXd& arr2) {
    // Check if both arrays are empty
    if (arr1.rows() > 0 && arr2.rows() > 0) {
        // Concatenate along the rows (axis 0)
        Eigen::MatrixXd result(arr1.rows() + arr2.rows(), arr1.cols());
        result << arr1, arr2;  // Concatenating along the vertical axis
        return result;
    } else if (arr1.rows() > 0) {
        // If only arr1 is non-empty, return it
        return arr1;
    } else if (arr2.rows() > 0) {
        // If only arr2 is non-empty, return it
        return arr2;
    } else {
        // If both arrays are empty, return an empty matrix
        return Eigen::MatrixXd(0, 0);
    }
}

static Geometry from_txt(const std::string& fname) {
    // Constants
    const std::unordered_map<std::string, int> valid_headers = {
        {"*nodes", 0}, {"*edges", 1}, {"*triangles", 2}
    };
    const std::vector<int> h_len = {3, 2, 3};  // Expected number of values per line
    const std::vector<int> h_dtype = {0, 1, 1};  // Representing GEOMETRY_FLOAT, GEOMETRY_INT
    const int GEOMETRY_INT = 1;
    const int GEOMETRY_FLOAT = 0;

    // Flags & parameters
    std::vector<bool> h_flag(valid_headers.size(), false);
    int cur_h = -1;  // Tracks current header
    std::vector<Eigen::MatrixXd> params(valid_headers.size());  // For storing matrices

    std::vector<std::string> temp_array;  // Temporary storage for values

    // Open file
    std::ifstream file(fname);
    if (!file.is_open()) {
        throw std::invalid_argument(fname + " is not a valid path");
    }

    std::string line;
    while (std::getline(file, line)) {
        // Trim whitespace and newlines
        line = trim(line);

        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Check for header lines (start with "*")
        if (line[0] == '*') {
            auto h_it = valid_headers.find(line);
            if (h_it == valid_headers.end()) {
                throw std::invalid_argument("Unknown header: " + line);
            }

            int h_id = h_it->second;

            if (h_flag[h_id]) {
                throw std::invalid_argument(line + " header used twice");
            }

            process_temp_array(cur_h, temp_array, params, h_dtype);
            h_flag[h_id] = true;
            cur_h = h_id;
        } else {  // Data line
            std::vector<std::string> vals = split(line, ",");
            if (vals.size() != h_len[cur_h]) {
                throw std::invalid_argument("Line " + line + " should have " + std::to_string(h_len[cur_h]) + " values");
            }

            temp_array.clear();
            for (const auto& val : vals) {
                temp_array.push_back(val);
            }
        }
    }

    // Process last collected data
    process_temp_array(cur_h, temp_array, params, h_dtype);

    return Geometry(params);
}

static void process_temp_array(int cur_h, const std::vector<std::string>& temp_array, std::vector<Eigen::MatrixXd>& params, const std::vector<int>& h_dtype) {
    if (!temp_array.empty()) {
        // Converting to Eigen matrix
        Eigen::MatrixXd temp_matrix(temp_array.size(), h_dtype[cur_h]);

        for (int i = 0; i < temp_array.size(); ++i) {
            temp_matrix(i) = std::stod(temp_array[i]);
        }

        if (h_dtype[cur_h] == 1) {  // GEOMETRY_INT
            temp_matrix.array() -= 1;  // Adjust for 0-based indexing
        }

        params[cur_h] = temp_matrix;
    }
}

// Utility function to trim whitespace from a string
static std::string trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\n\r");
    size_t end = str.find_last_not_of(" \t\n\r");
    return str.substr(start, end - start + 1);
}

// Utility function to split a string based on a delimiter
static std::vector<std::string> split(const std::string& str, const std::string& delimiter) {
    std::vector<std::string> tokens;
    size_t pos = 0;
    std::string token;
    while ((pos = str.find(delimiter)) != std::string::npos) {
        token = str.substr(0, pos);
        tokens.push_back(token);
        str.erase(0, pos + delimiter.length());
    }
    tokens.push_back(str);
    return tokens;
}

// TODO check function
std::pair<Eigen::MatrixXi, Eigen::MatrixXi> Geometry::separate_joint_edges(
    const Eigen::MatrixXi& triangles, const Eigen::MatrixXi& edges)
{
    // If no edges, return empty matrices
    if (edges.size() == 0) {
        return std::make_pair(Eigen::MatrixXi(), Eigen::MatrixXi());
    }

    // Create a set of unique shell nodes
    std::unordered_set<int> shell_nodes_set;
    for (int i = 0; i < triangles.rows(); ++i) {
        shell_nodes_set.insert(triangles(i, 0));
        shell_nodes_set.insert(triangles(i, 1));
        shell_nodes_set.insert(triangles(i, 2));
    }

    // Find joint edges (those that contain shell nodes)
    std::vector<Eigen::Vector2i> joint_edges;
    std::vector<Eigen::Vector2i> rod_edges;

    for (int i = 0; i < edges.rows(); ++i) {
        int node1 = edges(i, 0);
        int node2 = edges(i, 1);

        // If either node is in the shell, it's a joint edge
        if (shell_nodes_set.count(node1) > 0 || shell_nodes_set.count(node2) > 0) {
            joint_edges.push_back(Eigen::Vector2i(node1, node2));
        } else {
            rod_edges.push_back(Eigen::Vector2i(node1, node2));
        }
    }

    // Now handle reordering of the joint edges if necessary (if node1 is a shell node)
    for (auto& edge : joint_edges) {
        int node1 = edge(0);
        int node2 = edge(1);

        if (shell_nodes_set.count(node1) > 0) {
            // Reverse the edge if node1 is a shell node
            edge(0) = node2;
            edge(1) = node1;
        }
    }

    // Convert the vectors back to Eigen matrices
    Eigen::MatrixXi joint_edges_matrix(joint_edges.size(), 2);
    Eigen::MatrixXi rod_edges_matrix(rod_edges.size(), 2);

    for (int i = 0; i < joint_edges.size(); ++i) {
        joint_edges_matrix(i, 0) = joint_edges ;
        joint_edges_matrix(i, 1) = joint_edges ;
    }

    for (int i = 0; i < rod_edges.size(); ++i) {
        rod_edges_matrix(i, 0) = rod_edges ;
        rod_edges_matrix(i, 1) = rod_edges ;
    }

    return std::make_pair(joint_edges_matrix, rod_edges_matrix);
}



// Getter methods for Geometry Class
const Eigen::MatrixXd& nodes() const { return nodes_; }
const Eigen::MatrixXi& edges() const { return edges_; }
const Eigen::MatrixXi& rod_edges() const { return rod_edges_; }
const Eigen::MatrixXi& rod_shell_joint_edges() const { return rod_shell_joint_edges_; }
const Eigen::MatrixXi& face_nodes() const { return face_nodes_; }
