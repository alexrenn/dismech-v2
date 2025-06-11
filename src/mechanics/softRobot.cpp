#include "softRobot.h"
#include <unordered_set>
#include <iostream>

SoftRobot::SoftRobot(const GeomParams& geom, const Material& material, const Geometry& geo, 
                            const SimParams& sim_params, const Environment& env, const RobotState& state)
        : sim_params_(sim_params), env_(env), state_(state){

        _init_geometry(geo);
        _init_stiffness(geom, material);
        _init_state(geo);
        this->_init_springs(geo);
        _get_mass_matrix(geom, material);
        _get_voronoi_area();
    }

void SoftRobot::_init_geometry(const Geometry& geo)
 {
    // Initialize Geometry Properties
    n_nodes_ = geo.getNodes().size();
    n_edges_ = geo.getEdges().size();    
    n_edges_rod_only = geo.getRodEdges().size();
    n_edges_shell_only = geo.getShellEdges().size();
    n_edges_dof_ = n_edges_rod_only + geo.getRodShellJointEdgesTotal().size();
    n_faces_ = geo.getFaceNodes().size();
    face_nodes_shell_ = geo.getFaceNodes();
    face_edges_ = geo.getFaceEdges();
    nodes_ = geo.getNodes();
    edges_ = geo.getEdges();
    face_shell_edges_ = geo.getFaceShellEdges();
    twist_angles_ = geo.getTwistAngles();
    sign_ = geo.getBendTwistSigns();
    bend_twist_springs_ = geo.getBendTwistSprings();
    bend_twist_signs_ = geo.getBendTwistSigns();
    rod_stretch_springs_ = geo.getRodStretchSprings();
    shell_stretch_springs_ = geo.getShellStretchSprings();
    hinges_ = geo.getHinges();


    // Initialize DOF vector
    n_dof_ = 3 * n_nodes_ + n_edges_dof_;
    q0_ = Eigen::VectorXd::Zero(n_dof_);

    // Fill node positions
    for (int i = 0; i < n_nodes_; ++i) {
        q0_.segment<3>(3 * i) = nodes_[i];
    }

    // Fill twist angles
    if(twist_angles_.size() != 0) {
        for (int i = 0; i < n_edges_dof_; ++i) {
            q0_(3 * n_nodes_ + i) = twist_angles_[i];
        }
    }


// Midedge bending has more DOF
    if (sim_params_.use_mid_edge) {
        n_dof_ += n_edges_shell_only;
    
        Eigen::VectorXd extra = Eigen::VectorXd::Zero(n_edges_shell_only); // Create extra zeros
        Eigen::VectorXd q0_extended(n_dof_); // Concatenate the two vectors
        q0_extended << q0_, extra;
        q0_ = q0_extended;
    
        std::tie(init_ts_, init_fs_, init_cs_, init_xis_) = _init_curvature_midedge(geo);
    } else {
        tau0_.clear();
        init_ts_.clear();
        init_fs_.clear();
        init_cs_.clear();
        init_xis_.clear();
    }
    // Now compute reference lengths and face areas
    _get_ref_len();
    // TODO TEST
    _get_voronoi_ref_len();
    _get_face_area();

    // Testing
    std::cout << "\n=== Geometry Initialization ===" << std::endl;
    std::cout << "Number of nodes: " << n_nodes_ << std::endl;
    std::cout << "Number of edges (rod only): " << n_edges_rod_only << std::endl;
    std::cout << "Number of edges (shell only): " << n_edges_shell_only << std::endl;
    std::cout << "Number of edges (joint): " << (n_edges_dof_ - n_edges_rod_only) << std::endl;
    std::cout << "Number of edges: " << n_edges_ << std::endl;
    std::cout << "DOF: " << n_dof_ << std::endl;
    std::cout << "Number of faces: " << n_faces_ << std::endl;
 }

// Rod Stiffness and Shell Stiffness Initialization
void SoftRobot::_init_stiffness(const GeomParams& geom, const Material& material) {
    RodStiffness rod = computeRodStiffness(geom, material);
    this->EA = rod.EA;
    this->EI1 = rod.EI1;
    this->EI2 = rod.EI2;
    this->GJ = rod.GJ;

    ShellStiffness shell = computeShellStiffness(geom, material, this->ref_len, this->sim_params_.use_mid_edge);
    this->ks = shell.ks;
    this->kb = shell.kb;
    this->nu = material.poisson_shell;

    std::cout << "\n=== Stiffness Initialization ===" << std::endl;
    std::cout << "Rod stiffness (EA, EI1, EI2, GJ): " << EA << ", " << EI1 << ", " << EI2 << ", " << GJ << std::endl;
    std::cout << "Shell stiffness (ks size, kb): " << ks.size() << ", " << kb << std::endl;
    std::cout << "Poisson ratio: " << nu << std::endl;
}

void SoftRobot::_init_state(const Geometry& geo) {
    // Initialize RobotState state for q0
    auto [a1, a2] = this->compute_space_parallel(); 

    auto [m1, m2] = this->compute_material_directors(q0_, a1, a2);

    compute_time_parallel(a1, q0_, q0_);

    int n = static_cast<int>(geo.getBendTwistSprings().size());
    for (int i = 0; i < n; ++i) {
        edges_[i][0] = geo.getBendTwistSprings()[i][1];
        edges_[i][1] = geo.getBendTwistSprings()[i][3];
    }

    std::cout << edges_.size() << std::endl;
    std::cout << sign_.size() << std::endl;
    std::cout << a1.size() << std::endl;
    std::cout << tangent_.size() << std::endl;


    Eigen::VectorXi undef_ref_twist; 
    if (edges_.size() > 0) {
        Eigen::VectorXi zero_twist = Eigen::VectorXi::Zero(sign_.size());
        undef_ref_twist = compute_reference_twist_util(edges_, sign_, a1, tangent_, zero_twist);
    } else {
        undef_ref_twist.resize(0, 0); // Empty matrix
    }

    // Initialize robot state
    this->state_ = RobotState (
        this->q0_,
        Eigen::VectorXd::Zero(q0_.size()), // Empty vector for 'u'
        std::vector<Eigen::Vector3d>(), // Empty vector for 'a'
        a1, // computed above
        a2, // computed above
        m1, // computed above
        m2, // computed above
        undef_ref_twist, // need compute_reference_twist (springs)
        tau0_, // empty for now tau0_
        free_dof_);

    std::cout << "\n=== State Initialization ===" << std::endl;
    std::cout << "q0 size: " << q0_.size() << std::endl;
    std::cout << "Material directors (m1, m2) sizes: " << m1.size() << ", " << m2.size() << std::endl;
    std::cout << "Space parallel (a1, a2) sizes: " << a1.size() << ", " << a2.size() << std::endl;
    std::cout << "Reference twist size: " << undef_ref_twist.size() << std::endl;
}

void SoftRobot::_init_springs(const Geometry& geo) {
    // create lambda functions
    auto map_node_to_dof_vec = [this](const std::vector<int>& node_indices) -> std::vector<std::array<int, 3>> {
        return this->map_node_to_dof(node_indices);
    };
    auto map_node_to_dof_single = [this](int node_index) -> std::array<int, 3> {
        return this->map_node_to_dof({node_index})[0];
    };
    auto map_edge_to_dof_single = [this](size_t index) -> int {
        return 3 * this->n_nodes_ + static_cast<int>(index); };

    // Stretch springs - rod
    int n_rod = rod_stretch_springs_.size();
    for (int i = 0; i < n_rod; ++i) {
        rod_stretch_springs_object_.emplace_back(
            rod_stretch_springs_[i],
            this->ref_len[i],
            this->EA,
            map_node_to_dof_single
        );
    }

    // Stretch springs - shell
    int n_shell = shell_stretch_springs_.size();
    for (int i = 0; i < n_shell; ++i) {
        shell_stretch_springs_object_.emplace_back(
            shell_stretch_springs_[i],
            this->ref_len[i + n_rod],
            this->ks[i + n_rod],
             map_node_to_dof_single
        );
    }

    // Combine

    this->stretch_springs_.clear();
    this->stretch_springs_.reserve(rod_stretch_springs_object_.size() + shell_stretch_springs_object_.size());
    this->stretch_springs_.insert(this->stretch_springs_.end(), rod_stretch_springs_object_.begin(), rod_stretch_springs_object_.end());
    this->stretch_springs_.insert(this->stretch_springs_.end(), shell_stretch_springs_object_.begin(), shell_stretch_springs_object_.end());

    // Bend/twist springs
    this->bend_twist_springs_.clear();

    for (size_t i = 0; i < bend_twist_springs_.size(); ++i) {
        int edge_dof0 = map_edge_to_dof_single(i * 2);
        int edge_dof1 = map_edge_to_dof_single(i * 2 + 1);

        bend_twist_springs_object_.emplace_back(
            bend_twist_springs_[i],
            bend_twist_signs_[i],
            this->ref_len,
            std::vector<double>{this->EI1, this->EI2},
            this->GJ,
            map_node_to_dof_vec,
            [edge_dof0, edge_dof1](int which) -> int {
            return (which == 0) ? edge_dof0 : edge_dof1;
            }
        );
    }

    auto map_face_edge_to_dof_single = [this](int face_edge_index) -> int {
        int offset = 3 * this->n_nodes_ + this->n_edges_dof_;
        return offset + face_edge_index;
    };

    // Conditional: triangle or hinge springs
    if (this->sim_params_.use_mid_edge) {
        for (size_t i = 0; i < face_nodes_shell_.size(); ++i) {
            triangle_springs_object_.emplace_back(
                face_nodes_shell_[i],
                face_shell_edges_[i],
                face_edges_[i],
                sign_faces_[i],
                this->ref_len,
                std::vector<double>{this->face_area[i]},
                this->init_ts_[i],
                this->init_fs_[i],
                this->init_cs_[i],
                this->init_xis_[i],
                this->kb,
                this->nu,
                map_node_to_dof_single,
                map_face_edge_to_dof_single
            );
        }
        this->hinge_springs_.clear(); // No hinge springs
    } else {
        // Only hinge springs
        for (const auto& hinge : hinges_) {
            this->hinge_springs_object_.emplace_back(
                hinge,
                this->kb,
                map_node_to_dof_single
            );
        }
        this->triangle_springs_.clear(); // No triangle springs
    }
}
void SoftRobot::_get_mass_matrix(const GeomParams& geom, const Material& material) {
    mass_matrix = Eigen::VectorXd::Zero(this->n_dof_);  // using the member variable

    // Shell face contributions
    if (this->n_faces_ > 0) {
        std::vector<std::array<int, 3>>faces = this->face_nodes_shell_; // (n_faces, 3)

        // First compute all face areas and masses
        std::vector<double> face_areas(this->n_faces_);
        std::vector<double> face_masses(this->n_faces_);
        
        for (size_t i = 0; i < faces.size(); ++i) {
            const auto& face = faces[i];
            Eigen::Vector3d v1 = this->nodes_[face[1]] - this->nodes_[face[0]];
            Eigen::Vector3d v2 = this->nodes_[face[2]] - this->nodes_[face[1]];
            face_areas[i] = 0.5 * v1.cross(v2).norm();
            face_masses[i] = material.density * face_areas[i] * geom.shell_h;

            double node_mass = face_masses[i] / 3.0;  // Split mass equally among 3 nodes
            
            // Add mass to each node's DOFs
            for (int j = 0; j < 3; ++j) {
                int node_id = face[j];
                auto dof_indices = map_node_to_dof({node_id})[0];
                for (int dof : dof_indices) {
                    mass_matrix(dof) += node_mass;
                }
            }
        }
    }

    // Node contributions
    if (this->n_nodes_ > 0) {
        std::vector<double> dm_nodes(this->n_nodes_);
        double scale = 0.0; 
        if (geom.axs) {
            scale = *geom.axs * material.density;
        } else {
            double area = 3 * geom.rod_r0 * geom.rod_r0;
            scale = area * material.density;
        }

        // Manually fill in scaled values
        for (int i = 0; i < this->n_nodes_; ++i) {
            dm_nodes[i] = this->voronoi_ref_len[i] * scale;
            auto dof_indices = map_node_to_dof({i})[0]; // returns std::array<int, 3>
            for (int dof : dof_indices) {
                mass_matrix(dof) += dm_nodes[i];
            }
        }
    }

    // Edge contributions
    if (this->n_edges_dof_ > 0) {
        std::vector<double> dm_edges(this->n_edges_dof_);
        std::vector<double> edge_mass(this->n_edges_dof_);
        if (geom.axs) {
            for (int i = 0; i < this->n_edges_dof_; ++i) {
                dm_edges[i] = this->ref_len[i] * *geom.axs * material.density;
                edge_mass[i] = dm_edges[i] * (*geom.jxs / *geom.axs);
                int dof = map_edge_to_dof({i})(0);
                mass_matrix(dof) += edge_mass[i];
            }

        } else {
            double area = 3 * geom.rod_r0 * geom.rod_r0;
            for (int i = 0; i < this->n_edges_dof_; ++i) {
                dm_edges[i] = this->ref_len[i] * area * material.density;
                edge_mass[i] = dm_edges[i] * 0.5 * (geom.rod_r0 * geom.rod_r0);
                int dof = map_edge_to_dof({i})(0); // get the single DOF index
                mass_matrix(dof) += edge_mass[i];
            }
        }
    }

    std::cout << "\n=== Mass Matrix Initialization ===" << std::endl;
    std::cout << "Mass matrix size: " << mass_matrix.size() << std::endl;
    
    // Count non-zero elements
    int non_zero_count = 0;
    for (int i = 0; i < mass_matrix.size(); ++i) {
        if (mass_matrix(i) != 0) {
            non_zero_count++;
        }
    }
    std::cout << "Number of non-zero elements: " << non_zero_count << std::endl;

       std::cout << "Zero-mass DOF indices: ";
   for (int i = 0; i < mass_matrix.size(); ++i) {
       if (mass_matrix(i) == 0) std::cout << i << " ";
   }
   std::cout << std::endl;
    
    // Find min and max values
    double min_val = mass_matrix.minCoeff();
    double max_val = mass_matrix.maxCoeff();
    std::cout << "Mass matrix range: [" << min_val << ", " << max_val << "]" << std::endl;
    
    if (this->n_faces_ > 0) {
        std::cout << "Number of shell faces: " << this->n_faces_ << std::endl;
    }
    std::cout << "Number of nodes: " << this->n_nodes_ << std::endl;
    std::cout << "Number of edge DOFs: " << this->n_edges_dof_ << std::endl;
}

void SoftRobot::_get_ref_len() {
    this->ref_len = std::vector<double>();  // Allocate output vector
    this->ref_len.resize(n_edges_); // Allocate memory for n_edges_ elements

    for (size_t i = 0; i < n_edges_; ++i) {
        const Eigen::Vector3d& node1 = nodes_[edges_[i][0]];
        const Eigen::Vector3d& node2 = nodes_[edges_[i][1]];
        this->ref_len[i] = (node2 - node1).norm();
    }
}

void SoftRobot::_get_voronoi_ref_len() {
    std::vector<std::array<int, 2>> local_edges(edges_.begin(), edges_.begin() + n_edges_dof_);
    std::vector<double> weights(this->n_edges_dof_);
    for (int i = 0; i < this->n_edges_dof_; ++i) {
        weights[i] = 0.5 * this->ref_len[i];
    }
    
    std::vector<double> contributions(n_nodes_, 0.0); // Initialize with zeros

    for (size_t i = 0; i < local_edges.size(); ++i) {
    const auto& edge = local_edges[i];
            contributions[edge[0]] += weights[i]; // Add half length to first node
            contributions[edge[1]] += weights[i]; // Add half length to second node
    }
   
    this->voronoi_ref_len = contributions;
}

void SoftRobot::_get_voronoi_area() {
    if (this->face_nodes_shell_.size() == 0) {
        this->voronoi_area.clear(); 
        return; // Return empty vector if no faces
    }

    std::vector<std::array<int, 3>> faces = this->face_nodes_shell_;
    n_faces_ = faces.size();

    std::vector<double> areas(n_faces_);
    for (int i = 0; i < n_faces_; ++i) {
        Eigen::Vector3d a = this->nodes_[faces[i][0]]; 
        Eigen::Vector3d b = this->nodes_[faces[i][1]];
        Eigen::Vector3d c = this->nodes_[faces[i][2]];

        Eigen::Vector3d v1 = b - a;
        Eigen::Vector3d v2 = c - b;
        areas[i] = 0.5 * v1.cross(v2).norm();  // triangle area
    }
    std::vector<double> node_areas(this->n_nodes_, 0.0);
    for (int i = 0; i < n_faces_; ++i) {
        double share = areas[i] / 3.0;
        node_areas[faces[i][0]] += share;
        node_areas[faces[i][1]] += share;
        node_areas[faces[i][2]] += share;
    }

    this->voronoi_area = node_areas;
}

void SoftRobot::_get_face_area() {
    if (this->n_faces_ == 0) {
        this->face_area.clear();  // Clear if no faces
        return;
    }

    const std::vector<std::array<int, 3>>& faces = this->face_nodes_shell_;
    std::vector<double> areas(faces.size(), 0.0); 

    for (size_t i = 0; i < faces.size(); ++i) {
        const Eigen::Vector3d& a = this->nodes_[faces[i][0]];
        const Eigen::Vector3d& b = this->nodes_[faces[i][1]];
        const Eigen::Vector3d& c = this->nodes_[faces[i][2]];

        Eigen::Vector3d v1 = b - a;
        Eigen::Vector3d v2 = c - b;
        areas[i] = 0.5 * v1.cross(v2).norm();  // triangle area
    }

    this->face_area = areas;  // Save entire vector
}

void SoftRobot::scale_mass_matrix(const std::vector<Eigen::Vector3d>& nodes, double scale) {
    std::vector<int> node_indices;
    for (const auto& pos : nodes) {
        // Find the index of the node with the given position
        for (int i = 0; i < n_nodes_; ++i) {
            if (nodes_[i].isApprox(pos)) { // Use isApprox for floating-point comparison
                node_indices.push_back(i);
                break; // Found the node, move to the next position
            }
        }
    }

    std::vector<std::array<int, 3>> dof_indices = this->map_node_to_dof(node_indices);
    for (int i = 0; i < static_cast<int>(dof_indices.size()); ++i) {
        for (int j = 0; j < 3; ++j) {
            this->mass_matrix(dof_indices[i][j]) *= scale;
        }
    }
}

// ask radha types
std::tuple<std::vector<Eigen::Matrix3d>, std::vector<std::array<double, 3>>, std::vector<std::array<double, 3>>, std::vector<std::array<double, 3>>>
SoftRobot::_init_curvature_midedge(const Geometry& geo) {
    this->sign_faces_ = geo.getSignFaces();

    // Extract positions for each face
    std::vector<Eigen::Matrix3d> all_p_is;
    for (const auto& face : this->face_nodes_shell_) {
        Eigen::Matrix3d p_i;
        for (int j = 0; j < 3; j++) {
            p_i.col(j) = this->q0_.segment<3>(3 * face[j]);
        }
        all_p_is.push_back(p_i);
    }

    // Extract edge vectors
    std::vector<std::array<double, 3>> all_xi_is;
    for (const auto& edge : this->face_edges_) {
        std::array<double, 3> xi_i;
        for (int j = 0; j < 3; j++) {
            xi_i[j] = this->q0_(3 * edge[j] + j);  // Get x,y,z components
        }
        all_xi_is.push_back(xi_i);
    }

    // Compute initial tau0
    this->tau0_ = this->update_pre_comp_shell(this->q0_);
    
    // Convert tau0_ to Matrix3d format for each face edge
    std::vector<Eigen::Matrix3d> all_tau0_is;
    for (const auto& edge : this->face_edges_) {
        Eigen::Matrix3d tau0_i;
        for (int j = 0; j < 3; j++) {
            tau0_i.col(j) = this->tau0_[edge[j]];
        }
        all_tau0_is.push_back(tau0_i);
    }

    // Compute t, f, c for all faces simultaneously
    std::vector<Eigen::Matrix3d> t;
    std::vector<std::array<double, 3>> f;
    std::vector<std::array<double, 3>> c;

    compute_tfc_midedge(all_p_is, all_tau0_is, this->sign_faces_, t, f, c);

    return std::make_tuple(t, f, c, all_xi_is);
}

std::vector<Eigen::Vector3d> SoftRobot::update_pre_comp_shell(const Eigen::VectorXd& q) {
    if (!this->sim_params_.use_mid_edge) {
        return std::vector<Eigen::Vector3d>(); // Return empty vector if not using mid-edge
    }

    // Compute face normals
    std::vector<Eigen::Vector3d> face_normals(this->face_nodes_shell_.size());
    for (size_t i = 0; i < this->face_nodes_shell_.size(); ++i) {
        const auto& face = this->face_nodes_shell_[i];
        // Extract positions from q vector for each node
        Eigen::Vector3d q0(q[3 * face[0]], q[3 * face[0] + 1], q[3 * face[0] + 2]);
        Eigen::Vector3d q1(q[3 * face[1]], q[3 * face[1] + 1], q[3 * face[1] + 2]);
        Eigen::Vector3d q2(q[3 * face[2]], q[3 * face[2] + 1], q[3 * face[2] + 2]);

        Eigen::Vector3d v1 = q1 - q0;
        Eigen::Vector3d v2 = q2 - q1;
        face_normals[i] = v1.cross(v2).normalized();
    }

    // Accumulate edge normals
    std::vector<Eigen::Vector3d> edge_normals(this->n_edges_, Eigen::Vector3d::Zero());
    std::vector<double> edge_counts(this->n_edges_, 0);

    for (int i = 0; i < this->face_edges_.size(); ++i) {
        int e0 = this->face_edges_[i][0];
        int e1 = this->face_edges_[i][1];

        edge_counts[e0] += 1;
        edge_counts[e1] += 1;

        edge_normals[e0] += face_normals[i];
        edge_normals[e1] += face_normals[i];
    }

    for(int i = 0; i < this->n_edges_; ++i) {
        if (edge_counts[i] > 0) {
            edge_normals[i] /= edge_counts[i];
        }
    }

    // Compute edge vectors and tau_0
    std::vector<Eigen::Vector3d> tau0(this->edges_.size());
    for (size_t i = 0; i < this->edges_.size(); ++i) {
        int n0 = this->edges_[i][0];
        int n1 = this->edges_[i][1];

        // Extract positions from q vector
        Eigen::Vector3d q0(q[3 * n0], q[3 * n0 + 1], q[3 * n0 + 2]);
        Eigen::Vector3d q1(q[3 * n1], q[3 * n1 + 1], q[3 * n1 + 2]);
        Eigen::Vector3d edge_vec = q1 - q0;

        tau0[i] = edge_vec.cross(edge_normals[i]).normalized();
    }

    return tau0;
}

std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>>
SoftRobot::compute_space_parallel() {
    std::vector<Eigen::Vector3d> a1(this->n_edges_dof_, Eigen::Vector3d::Zero());
    std::vector<Eigen::Vector3d> a2(this->n_edges_dof_, Eigen::Vector3d::Zero());

    this->tangent_ = this->_compute_tangent(this->q0_); // tangent_ is std::vector<Eigen::Vector3d>

    if (!this->tangent_.empty()) {
        // Initialize first a1
        Eigen::Vector3d ref_y(0, 1, 0);
        Eigen::Vector3d ref_z(0, 0, -1);

        Eigen::Vector3d a1_init = this->tangent_[0].cross(ref_y);
        if (a1_init.norm() < 1e-6) {
            a1_init = this->tangent_[0].cross(ref_z);
        }
        a1[0] = a1_init.normalized();
        a2[0] = this->tangent_[0].cross(a1[0]);

        // Iterative parallel transport
        for (int i = 1; i < this->n_edges_dof_; ++i) {
            const Eigen::Vector3d& t_prev = this->tangent_[i - 1];
            const Eigen::Vector3d& t_curr = this->tangent_[i];
            const Eigen::Vector3d& a1_prev = a1[i - 1];

            Eigen::Vector3d transported = parallel_transport(a1_prev, t_prev, t_curr);

            // Remove component along t_curr
            transported -= transported.dot(t_curr) * t_curr;
            a1[i] = transported.normalized();
            a2[i] = t_curr.cross(a1[i]);
        }
    }

    return std::make_pair(a1, a2);
}


std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>>
SoftRobot::compute_time_parallel(const std::vector<Eigen::Vector3d>& a1_old,
                                 const Eigen::VectorXd& q0,
                                 const Eigen::VectorXd& q) {
    // Get tangent vectors (as vectors of Vector3d)
    std::vector<Eigen::Vector3d> tangent0 = this->_compute_tangent(q0);
    std::vector<Eigen::Vector3d> tangent  = this->_compute_tangent(q);

    // Parallel transport a1_old from tangent0 to tangent
    auto a1_transported = parallel_transport(a1_old, tangent0, tangent);

    size_t N = tangent.size();
    std::vector<Eigen::Vector3d> a1(N);
    std::vector<Eigen::Vector3d> a2(N);

    for (size_t i = 0; i < N; ++i) {
        double t_dot = a1_transported[i].dot(tangent[i]);
        a1[i] = a1_transported[i] - t_dot * tangent[i];
        a1[i].normalize();  // Normalize in-place
        a2[i] = tangent[i].cross(a1[i]);
    }

    std::cout << a1.size() << std::endl;
    std::cout << a2.size() << std::endl;

    return std::make_pair(a1, a2);
}
std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> SoftRobot::compute_material_directors(
    const Eigen::VectorXd& q,
    const std::vector<Eigen::Vector3d>& a1,
    const std::vector<Eigen::Vector3d>& a2)
{
    Eigen::VectorXd theta = q.segment(3 * n_nodes_, n_edges_dof_);
    Eigen::VectorXd cos_theta = theta.array().cos();
    Eigen::VectorXd sin_theta = theta.array().sin();

    std::vector<Eigen::Vector3d> m1(n_edges_dof_);
    std::vector<Eigen::Vector3d> m2(n_edges_dof_);

    for (int i = 0; i < n_edges_dof_; ++i) {
        m1[i] = cos_theta(i) * a1[i] + sin_theta(i) * a2[i];
        m2[i] = -sin_theta(i) * a1[i] + cos_theta(i) * a2[i];
    }

    return std::make_pair(m1, m2);
}

Eigen::VectorXi SoftRobot::compute_reference_twist(
    const std::vector<BendTwistSpring>& springs,
    const Eigen::VectorXd& q,
    const std::vector<Eigen::Vector3d>& a1,
    const Eigen::VectorXi& ref_twist)
{
    if (springs.empty()) {
        return Eigen::VectorXi();  // empty
    }

    std::vector<std::array<int, 2>> edges;
    std::vector<std::array<int, 2>> sgn;

    edges.reserve(springs.size());
    sgn.reserve(springs.size());

    for (const auto& s : springs) {
        edges.push_back(s.edges_ind);
        sgn.push_back(s.sgn);
    }

    std::vector<Eigen::Vector3d> tangent = this->_compute_tangent(q);

    return compute_reference_twist_util(edges, sgn, a1, tangent, ref_twist);
}


std::vector<Eigen::Vector3d> SoftRobot::_compute_tangent(const Eigen::VectorXd& q) {
    std::vector<Eigen::Vector3d> tangent(this->n_edges_dof_);

    for (int i = 0; i < this->n_edges_dof_; ++i) {
        int n0 = this->edges_[i][0];
        int n1 = this->edges_[i][1];

        // Convert node indices to DOF triplets
        std::array<int, 3> dof0 = map_node_to_dof({n0})[0];
        std::array<int, 3> dof1 = map_node_to_dof({n1})[0];

        Eigen::Vector3d pos0(q[dof0[0]], q[dof0[1]], q[dof0[2]]);
        Eigen::Vector3d pos1(q[dof1[0]], q[dof1[1]], q[dof1[2]]);

        Eigen::Vector3d edge = pos1 - pos0;
        double norm = edge.norm();
        if (norm < 1e-10) norm = 1.0;

        Eigen::Vector3d unit_vec = edge / norm;

        for (int j = 0; j < 3; ++j) {
            if (std::abs(unit_vec(j)) < _TANGENT_THRESHOLD)
                unit_vec(j) = 0.0;
        }

        tangent[i] = unit_vec;
    }

    return tangent;
}

// Fix/free nodes and edges
// ask radha should i be passing anything in here
void SoftRobot::free_nodes(const std::vector<int>& nodes, std::optional<int> axis, bool fix_edges) {
    Eigen::VectorXi node_dof_mask = _get_node_dof_mask(nodes, axis);
    Eigen::VectorXi new_dof = setdiff1d(fixed_dof_, node_dof_mask);

    if (fix_edges) {
        Eigen::VectorXi edge_dof_mask = _get_intermediate_edge_dof(nodes);
        new_dof = setdiff1d(new_dof, edge_dof_mask);
    }

    _fix_dof(new_dof);
}

void SoftRobot::fix_nodes(const std::vector<int>& nodes, std::optional<int> axis, bool fix_edges) {
    Eigen::VectorXi node_dof_mask = _get_node_dof_mask(nodes, axis);
    Eigen::VectorXi new_dof = setdiff1d(get_fixed_dof(), node_dof_mask);

    if (fix_edges) {
        Eigen::VectorXi edge_dof_mask = _get_intermediate_edge_dof(nodes);
        new_dof = setdiff1d(new_dof, edge_dof_mask);
    }

    _fix_dof(new_dof);
}

void SoftRobot::free_edges(const std::vector<int>& edge_indices) {
    Eigen::VectorXi edge_dofs = map_edge_to_dof(edge_indices);
    Eigen::VectorXi new_dof = setdiff1d(fixed_dof_, edge_dofs);
    _fix_dof(new_dof);
}

void SoftRobot::fix_edges(const std::vector<int>& edge_indices) {
    Eigen::VectorXi edge_dofs = map_edge_to_dof(edge_indices);
    Eigen::VectorXi new_dof = union1d(fixed_dof_, edge_dofs);
    _fix_dof(new_dof);
}


void SoftRobot::_fix_dof(const Eigen::VectorXi& new_fixed_dof) {
    Eigen::VectorXi all_dof = Eigen::VectorXi::LinSpaced(n_dof_, 0, n_dof_ - 1);
    free_dof_ = setdiff1d(all_dof, new_fixed_dof);
    fixed_dof_ = new_fixed_dof;
}


Eigen::VectorXi SoftRobot::_get_intermediate_edge_dof(const std::vector<int>& nodes) {
    std::vector<int> edge_indices;
    for (int i = 0; i < this->n_edges_dof_; ++i) {
        int node1_index = this->edges_[i][0];
        int node2_index = this->edges_[i][1];
        bool node1_found = std::find(nodes.begin(), nodes.end(), node1_index) != nodes.end();
        bool node2_found = std::find(nodes.begin(), nodes.end(), node2_index) != nodes.end();
        if (node1_found && node2_found) {
            edge_indices.push_back(i);
        }
    }
    if (this->sim_params_.two_d_sim) {
        for (int i = 0; i < this->n_edges_dof_; ++i) {
            if (std::find(edge_indices.begin(), edge_indices.end(), i) == edge_indices.end()) {
                edge_indices.push_back(i);
            }
        }
    }
    return map_edge_to_dof(edge_indices);
}

// Utility Helpers
Eigen::VectorXi SoftRobot::setdiff1d(const Eigen::VectorXi& a, const Eigen::VectorXi& b) {
    std::unordered_set<int> b_set(b.data(), b.data() + b.size());
    std::vector<int> result;
    for (int i = 0; i < a.size(); ++i)
        if (b_set.find(a(i)) == b_set.end())
            result.push_back(a(i));
    return Eigen::Map<Eigen::VectorXi>(result.data(), result.size());
}

Eigen::VectorXi SoftRobot::union1d(const Eigen::VectorXi& a, const Eigen::VectorXi& b) {
    std::unordered_set<int> union_set(a.data(), a.data() + a.size());
    for (int i = 0; i < b.size(); ++i)
        union_set.insert(b(i));
    std::vector<int> result(union_set.begin(), union_set.end());
    std::sort(result.begin(), result.end());
    return Eigen::Map<Eigen::VectorXi>(result.data(), result.size());
}

// Check if each a is in b
std::vector<bool> isin(const Eigen::VectorXi& a, const Eigen::VectorXi& b) {
    std::unordered_set<int> b_set(b.data(), b.data() + b.size());
    std::vector<bool> result(a.size());

    for (int i = 0; i < a.size(); ++i) {
        result[i] = b_set.count(a[i]) > 0;
    }
    return result;
}

// Perturb system
void SoftRobot::move_nodes(const std::vector<int>& nodes, const Eigen::VectorXd& perturbation, std::optional<int> axis) {
    Eigen::VectorXi dof_mask = _get_node_dof_mask(nodes, axis);
    state_.q0(dof_mask) += perturbation;
}

void SoftRobot::twist_edges(const std::vector<int>& edge_indices, const Eigen::VectorXd& perturbation) {
    Eigen::VectorXi dof_mask = map_edge_to_dof(edge_indices);
    state_.q0(dof_mask) += perturbation;
}

//Utility
Eigen::VectorXi SoftRobot::_get_node_dof_mask(const std::vector<int>& nodes, std::optional<int> axis) const {
    std::vector<std::array<int, 3>> node_dof_vec = map_node_to_dof(nodes);
    Eigen::MatrixXi node_dof(node_dof_vec.size(), 3);
    for (int i = 0; i < node_dof_vec.size(); ++i) {
        node_dof(i, 0) = node_dof_vec[i][0];
        node_dof(i, 1) = node_dof_vec[i][1];
        node_dof(i, 2) = node_dof_vec[i][2];
    }

    if (axis.has_value()) {
        // Return just one column (axis-specified)
        return node_dof.col(axis.value());
    } else {
        // Flatten the entire matrix row-wise
        Eigen::VectorXi flat(node_dof.rows() * node_dof.cols());
        int index = 0;
        for (int i = 0; i < node_dof.rows(); ++i) {
            for (int j = 0; j < node_dof.cols(); ++j) {
                flat(index++) = node_dof(i, j);
            }
        }
        return flat;
    }
}
// Maps a list of node numbers to their 3 DOFs each: [x, y, z]
std::vector<std::array<int, 3>> SoftRobot::map_node_to_dof(const std::vector<int>& node_indices) {
    int num_nodes = node_indices.size();
    std::vector<std::array<int, 3>> dof_map(num_nodes);
    for (int i = 0; i < num_nodes; ++i) {
        int base = 3 * node_indices[i];
        dof_map[i][0] = base;
        dof_map[i][1] = base + 1;
        dof_map[i][2] = base + 2;
    }

    return dof_map;
}

Eigen::VectorXi SoftRobot::map_edge_to_dof(const std::vector<int>& edge_indices) const {
    int num_edges = edge_indices.size();
    Eigen::VectorXi dof_indices(num_edges);
    int offset = 3 * this->n_nodes_;
    for (int i = 0; i < num_edges; ++i) {
        dof_indices(i) = offset + edge_indices[i];
    }
    return dof_indices;
}

// Maps face-edge indices to global DOF indices
Eigen::VectorXi SoftRobot::map_face_edge_to_dof(const std::vector<int>& edge_nums) const {
    int num_edges = edge_nums.size();
    Eigen::VectorXi dof_indices(num_edges);
    int offset = 3 * n_nodes_ + n_edges_dof_;

    for (int i = 0; i < num_edges; ++i) {
        dof_indices(i) = offset + edge_nums[i];
    }

    return dof_indices;
}

    // Getter for node_dof_indices (Node DOF indices matrix)
    std::vector<std::array<int, 3>> SoftRobot::get_node_dof_indices() const {
        std::vector<std::array<int, 3>> indices(n_nodes_);
        for (int i = 0; i < n_nodes_; ++i) {
            indices[i] = {3 * i, 
                          3 * i + 1, 
                          3 * i + 2};
        }
        return indices;
    }

    // Getter for fixed_dof (indices of constrained degrees of freedom)
Eigen::VectorXi SoftRobot::get_fixed_dof() const {
    Eigen::VectorXi fixed_dof(fixed_dof_.size());
    for (int i = 0; i < fixed_dof_.size(); ++i) {
        fixed_dof(i) = fixed_dof_[i];
    }
    return fixed_dof;
}

// Debug function to print all member variables and test all public member functions
void SoftRobot::debug() {
    std::cout << "\n===== SoftRobot Debug Info =====" << std::endl;
    std::cout << "n_nodes_: " << n_nodes_ << std::endl;
    std::cout << "n_edges_rod_only: " << n_edges_rod_only << std::endl;
    std::cout << "n_edges_shell_only: " << n_edges_shell_only << std::endl;
    std::cout << "n_edges_: " << n_edges_ << std::endl;
    std::cout << "n_edges_dof_: " << n_edges_dof_ << std::endl;
    std::cout << "n_faces_: " << n_faces_ << std::endl;
    std::cout << "n_dof_: " << n_dof_ << std::endl;
    std::cout << "EA: " << EA << ", EI1: " << EI1 << ", EI2: " << EI2 << ", GJ: " << GJ << std::endl;
    std::cout << "kb: " << kb << ", nu: " << nu << std::endl;
    std::cout << "nodes_.size(): " << nodes_.size() << std::endl;
    std::cout << "edges_.size(): " << edges_.size() << std::endl;
    std::cout << "face_nodes_shell_.size(): " << face_nodes_shell_.size() << std::endl;
    std::cout << "face_edges_.size(): " << face_edges_.size() << std::endl;
    std::cout << "sign_faces_.size(): " << sign_faces_.size() << std::endl;
    std::cout << "twist_angles_.size(): " << twist_angles_.size() << std::endl;
    std::cout << "sign_.size(): " << sign_.size() << std::endl;
    std::cout << "hinges_.size(): " << hinges_.size() << std::endl;
    std::cout << "ref_len.size(): " << ref_len.size() << std::endl;
    std::cout << "voronoi_ref_len.size(): " << voronoi_ref_len.size() << std::endl;
    std::cout << "voronoi_area.size(): " << voronoi_area.size() << std::endl;
    std::cout << "face_area.size(): " << face_area.size() << std::endl;
    std::cout << "mass_matrix size: " << mass_matrix.size() << std::endl;
    std::cout << "rod_stretch_springs_.size(): " << rod_stretch_springs_.size() << std::endl;
    std::cout << "shell_stretch_springs_.size(): " << shell_stretch_springs_.size() << std::endl;
    std::cout << "bend_twist_springs_.size(): " << bend_twist_springs_.size() << std::endl;
    std::cout << "triangle_springs_.size(): " << triangle_springs_.size() << std::endl;
    std::cout << "hinge_springs_.size(): " << hinge_springs_.size() << std::endl;
    std::cout << "bend_twist_signs_.size(): " << bend_twist_signs_.size() << std::endl;
    std::cout << "rod_stretch_springs_object_.size(): " << rod_stretch_springs_object_.size() << std::endl;
    std::cout << "shell_stretch_springs_object_.size(): " << shell_stretch_springs_object_.size() << std::endl;
    std::cout << "bend_twist_springs_object_.size(): " << bend_twist_springs_object_.size() << std::endl;
    std::cout << "triangle_springs_object_.size(): " << triangle_springs_object_.size() << std::endl;
    std::cout << "hinge_springs_object_.size(): " << hinge_springs_object_.size() << std::endl;
    std::cout << "stretch_springs_.size(): " << stretch_springs_.size() << std::endl;
    std::cout << "tau0_.size(): " << tau0_.size() << std::endl;
    std::cout << "init_ts_.size(): " << init_ts_.size() << std::endl;
    std::cout << "init_fs_.size(): " << init_fs_.size() << std::endl;
    std::cout << "init_cs_.size(): " << init_cs_.size() << std::endl;
    std::cout << "init_xis_.size(): " << init_xis_.size() << std::endl;
    std::cout << "tangent_.size(): " << tangent_.size() << std::endl;
    std::cout << "free_dof_.size(): " << free_dof_.size() << std::endl;
    std::cout << "fixed_dof_.size(): " << fixed_dof_.size() << std::endl;
    std::cout << "q0_.size(): " << q0_.size() << std::endl;
    // Test public getters
    std::cout << "get_ref_len().size(): " << get_ref_len().size() << std::endl;
    std::cout << "get_voronoi_ref_len().size(): " << get_voronoi_ref_len().size() << std::endl;
    std::cout << "get_voronoi_area().size(): " << get_voronoi_area().size() << std::endl;
    std::cout << "get_face_area().size(): " << get_face_area().size() << std::endl;
    std::cout << "get_mass_matrix().size(): " << get_mass_matrix().size() << std::endl;
    std::cout << "get_node_dof_indices().size(): " << get_node_dof_indices().size() << std::endl;
    std::cout << "get_fixed_dof().size(): " << get_fixed_dof().size() << std::endl;
    std::cout << "get_end_node_dof_index(): " << get_end_node_dof_index() << std::endl;
    std::cout << "get_q0().size(): " << get_q0().size() << std::endl;
    std::cout << "get_n_dof(): " << get_n_dof() << std::endl;
    // Test a few utility functions
    if (!nodes_.empty()) {
        std::vector<int> test_nodes = {0};
        auto dof_mask = _get_node_dof_mask(test_nodes, std::nullopt);
        std::cout << "_get_node_dof_mask({0}): size " << dof_mask.size() << std::endl;
        auto dof_map = map_node_to_dof(test_nodes);
        std::cout << "map_node_to_dof({0}): {" << dof_map[0][0] << ", " << dof_map[0][1] << ", " << dof_map[0][2] << "}" << std::endl;
    }
    if (!edges_.empty()) {
        std::vector<int> test_edge_indices = {0};
        auto edge_dof = map_edge_to_dof(test_edge_indices);
        std::cout << "map_edge_to_dof({first edge}): size " << edge_dof.size() << std::endl;
    }
    if (!face_edges_.empty()) {
        std::vector<int> test_face_edges = {0};
        auto face_edge_dof = map_face_edge_to_dof(test_face_edges);
        std::cout << "map_face_edge_to_dof({0}): size " << face_edge_dof.size() << std::endl;
    }
    std::cout << "===== End SoftRobot Debug Info =====\n" << std::endl;
}

// // main for testing
// int main() {
//     std::cout << "SoftRobot test initialization.\n";

//     // Create GeomParams with same values
//     GeomParams geom_params(1e-3,    // rod_r0
//                           1e-3    // shell_h
//                           );    // axs (no area scaling)

//     // Create Material with same values
//     Material material(1500,    // density
//                      10e6,    // youngs_rod
//                      10e8,    // youngs_shell
//                      0.5,     // poisson_rod
//                      0.30);    // poisson_shell

//     // Create SimParams with same values
//     SimParams dynamic_3d_sim(false,  // static_sim
//                            false,    // two_d_sim
//                            false,    // use_mid_edge
//                            false,    // use_line_search
//                            false,    // show_floor
//                            true,     // log_data
//                            1,        // log_step
//                            1e-2,     // dt
//                            25,       // max_iter
//                            3.0,      // total_time
//                            10,       // plot_step
//                            1e-4,     // tol
//                            1e-4,     // ftol
//                            1e-2);    // dtol

//     // Create Environment and add forces
//     Environment env;
//     // env.add_force("gravity", Eigen::Vector3d(0.0, 0.0, -9.81));
//     // env.add_force("aerodynamics", 1.0, 10.0);  // rho=1, cd=10

//     // Load geometry from file
//     Geometry geo("src/mechanics/hex_parachute_n6.txt");

//     // Create initial RobotState with required parameters
//     Eigen::VectorXd q0 = Eigen::VectorXd::Zero(1);  // Will be resized by SoftRobot
//     Eigen::VectorXd u = Eigen::VectorXd::Zero(1);   // Will be resized by SoftRobot
//     std::vector<Eigen::Vector3d> a;  // Empty vector
//     std::vector<Eigen::Vector3d> a1; // Empty vector
//     std::vector<Eigen::Vector3d> a2; // Empty vector
//     std::vector<Eigen::Vector3d> m1; // Empty vector
//     std::vector<Eigen::Vector3d> m2; // Empty vector
//     Eigen::VectorXi undef_ref_twist = Eigen::VectorXi::Zero(1); // Will be resized
//     std::vector<Eigen::Vector3d> tau0; // Empty vector
//     Eigen::VectorXi free_dof = Eigen::VectorXi::Zero(1); // Will be resized

//     RobotState initial_state(q0, u, a, a1, a2, m1, m2, undef_ref_twist, tau0, free_dof);

//     // Create SoftRobot instance
//     SoftRobot robot(geom_params, material, geo, dynamic_3d_sim, env, initial_state);

//     std::cout << "SoftRobot initialized successfully.\n";
//     robot.debug();
//     return 0;
// }