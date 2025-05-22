#include "softRobot.h"
#include <unordered_set>

SoftRobot::SoftRobot(const GeomParams& geom, const Material& material, const Geometry& geo, 
                            const SimParams& sim_params, const Environment& env, const RobotState& state)
        : sim_params_(sim_params), env_(env), state_(state){

        this->_init_geometry(geo);
        this->_init_stiffness(geom, material);
        this->_init_state(geo);
        //this->_init_springs(geo);
        _get_mass_matrix(geom, material);
    }

void SoftRobot::_init_geometry(const Geometry& geo)
 {
    // Initialize Geometry Properties

    n_nodes_ = geo.getNodes().size();
    n_edges_ = geo.getEdges().size();    
    n_edges_rod_only = geo.getRodEdges().size();
    n_edges_shell_only = geo.getShellEdges().size();
    int n_edges_joint = geo.getRodShellJointEdgesTotal().size();
    n_edges_dof_ = n_edges_rod_only + n_edges_joint;
    n_faces_ = geo.getFaceNodes().size();
    face_edges_ = geo.getFaceEdges();
    nodes_ = geo.getNodes();
    edges_ = geo.getEdges();
    face_nodes_shell_ = geo.getFaceNodes();
    twist_angles_ = geo.getTwistAngles();

    // Initialize DOF vector
    n_dof_ = 3 * n_nodes_ + n_edges_dof_;
    q0_ = Eigen::VectorXd::Zero(n_dof_);
    
    // Fill node positions
    for (int i = 0; i < n_nodes_; ++i) {
        q0_.segment<3>(3 * i) = nodes_[i];
    }
    // Fill twist angles
    for (int i = 0; i < n_edges_dof_; ++i) {
        q0_(3 * n_nodes_ + i) = twist_angles_[i];
    }


    // TODO: temporary until compute_tfc_midedge is implemented
    Eigen::VectorXi temp = Eigen::VectorXi::Zero(n_dof_); 
    // tau0_.resize(n_dof_, temp); // Initialize tau0_ with zeros

   /* // Midedge bending has more DOF
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
    */
    // Now compute reference lengths and face areas
    _get_ref_len();
    _get_voronoi_ref_len();
    _get_face_area();
 }
/*
 // TOOD: rod sitffness
void SoftRobot::_init_stiffness(const GeomParams& geom, const Material& material) {
    RodStiffness rod = compute_rod_stiffness(geom, material);
    this->EA = rod.EA;
    this->EI1 = rod.EI1;
    this->EI2 = rod.EI2;`
    this->GJ = rod.GJ;

    RodStiffness rod = compute_shell_stiffness(geom, material, ref_len, sim_params.use_mid_edge);
    this->ks = shell.ks;
    this->kb = shell.kb;
}
*/
void SoftRobot::_init_state(const Geometry& geo) {
    // Initialize RobotState state for q0
    auto [a1, a2] = this->compute_space_parallel(); // TODO parallel_transport

    auto [m1, m2] = this->compute_material_directors(q0_, a1, a2);
/*
    int n = static_cast<int>(geo.getBendTwistSprings().size());
    for (int i = 0; i < n; ++i) {
        edges_[i][0] = geo.getBendTwistSprings()[i][1];
        edges_[i][1] = geo.getBendTwistSprings()[i][3];
    }



    std::vector<int> sign(n, 0);
    for (int i = 0; i < n; ++i) {
        sign[i] = bend_twist_signs_[i];  // Assuming you have bend_twist_signs_
    }
        */

    /* std::vector<double> undef_ref_twist; 
    if (edges_.size() > 0) {
        Eigen::VectorXd zero_twist = Eigen::VectorXd::Zero(sign.size());
        undef_ref_twist = compute_reference_twist(
            edges_, sign, a1, this->tangent_, zero_twist);
    } else {
        undef_ref_twist.resize(0, 0); // Empty matrix
    } */


    // TODO temp until springs are implemented
    Eigen::VectorXi temp = Eigen::VectorXi::Zero(edges_.size());  // unsure if this is right size
    const Eigen::VectorXi& undef_ref_twist = temp;
    

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
        Eigen::VectorXd::Zero(n_dof_), // empty for now tau0_
        free_dof_);
}
/*
void SoftRobot::_init_springs(const Geometry& geo) {
    int n_rod = geo.rod_stretch_springs.size();

    // Stretch springs - rod
    std::vector<StretchSpring> rod_springs;
    for (int i = 0; i < n_rod; ++i) {
        rod_springs.emplace_back(
            geo.rod_stretch_springs[i],
            this->ref_len(i),
            this->EA,
            this->map_node_to_dof
        );
    }

    // Stretch springs - shell
    std::vector<StretchSpring> shell_springs;
    int n_shell = geo.shell_stretch_springs.size();
    for (int i = 0; i < n_shell; ++i) {
        shell_springs.emplace_back(
            geo.shell_stretch_springs[i],
            this->ref_len(i + n_rod),
            this->ks(i + n_rod),
            this->map_node_to_dof
        );
    }

    // Combine
    this->stretch_springs_.reserve(rod_springs.size() + shell_springs.size());
    this->stretch_springs_.insert(this->stretch_springs_.end(),
                                 rod_springs.begin(), rod_springs.end());
    this->stretch_springs_.insert(this->stretch_springs_.end(),
                                 shell_springs.begin(), shell_springs.end());

    // Bend/twist springs
    this->bend_twist_springs_.clear();
    for (size_t i = 0; i < geo.bend_twist_springs.size(); ++i) {
        this->bend_twist_springs_.emplace_back(
            geo.bend_twist_springs[i],
            geo.bend_twist_signs[i],
            this->ref_len,
            Eigen::Vector2d(this->EI1, this->EI2),
            this->GJ,
            this->map_node_to_dof,
            this->map_edge_to_dof
        );
    }

    // Conditional: triangle or hinge springs
    if (this->sim_params.use_mid_edge) {
        this->triangle_springs_.clear();
        for (size_t i = 0; i < geo.face_nodes.size(); ++i) {
            this->triangle_springs_.emplace_back(
                geo.face_nodes[i],
                geo.face_shell_edges[i],
                geo.face_edges_[i],
                geo.sign_faces[i],
                this->ref_len,
                this->face_area_[i],
                this->init_ts_[i]
                this->init_fs_[i]
                this->init_cs_[i]
                this->init_xis_[i]
                this->kb,
                this->nu,
                this->map_node_to_dof,
                this->map_face_edge_to_dof
            );
        }
        this->shell_hinge_springs_.clear(); // No hinge springs
    } else {
        // Only hinge springs
        this->shell_hinge_springs_.clear();
        for (const auto& hinge : geo.hinges) {
            this->shell_hinge_springs_.emplace_back(
                hinge,
                this->kb,
                this->map_node_to_dof
            );
        }
        this->triangle_springs_.clear(); // No triangle springs
    }
} */
void SoftRobot::_get_mass_matrix(const GeomParams& geom, const Material& material) {
    mass_matrix = Eigen::VectorXd::Zero(this->n_dof_);  // using the member variable

    // Shell face contributions
    if (this->n_faces_ > 0) {
        std::vector<std::array<int, 3>>faces = this->face_nodes_shell_; // (n_faces, 3)

        std::vector<double> areas(this->n_faces_);
        std::vector<double> m_shell(this->n_faces_);; 

        for (size_t i = 0; i < faces.size(); ++i) {
            const auto& face = faces[i];
            Eigen::Vector3d v1 = this->nodes_[face[1]] - this->nodes_[face[0]];
            Eigen::Vector3d v2 = this->nodes_[face[2]] - this->nodes_[face[1]];

            double area = 0.5 * v1.cross(v2).norm();
            for (int j = 0; j < 3; ++j) {
                int node_id = face[j];
                for (int k = 0; k < 3; ++k) {
                    int dof = 3 * node_id + k;
                    mass_matrix(dof) +=material.density * area * geom.shell_h / 3.0;
                }
            }
        }
    }

    // Node contributions
    if (this->n_nodes_ > 0) {
        std::vector<double> dm_nodes(this->n_nodes_);
        double scale = 0.0; 
        if (geom.axs >= 0.0) {
            scale = geom.axs * material.density;
        } else {
            double area = M_PI * geom.rod_r0 * geom.rod_r0;
            scale = area * material.density;
        }

        // Manually fill in scaled values
        for (int i = 0; i < this->n_nodes_; ++i) {
            dm_nodes[i] = this->voronoi_ref_len[i] * scale;
        }

        for (int i = 0; i < this->n_nodes_; ++i) {
            for (int j = 0; j < 3; ++j) {
                int dof = 3 * i + j;
                mass_matrix(dof) += dm_nodes[i];
            }
        }
    }

    // Edge contributions
    // TODO match python code
    if (this->n_edges_dof > 0) {
        std::vector<double> dm_edges(this->n_edges_dof);
        std::vector<double> edge_mass(this->n_edges_dof);
        if (geom.axs >= 0.0) {
            for (int i = 0; i < this->n_edges_dof; ++i) {
                dm_edges[i] = this->ref_len[i] * geom.axs * material.density;
            }
            for (int i = 0; i < this->n_edges_dof; ++i) {
                edge_mass[i] = dm_edges[i] * 0.5 * (geom.rod_r0 * geom.rod_r0);
            }

        } else {
            double area = M_PI * geom.rod_r0 * geom.rod_r0;
            for (int i = 0; i < this->n_edges_dof; ++i) {
                dm_edges[i] = this->ref_len[i] * area * material.density;
            }
            for (int i = 0; i < this->n_edges_dof; ++i) {
                edge_mass[i] = dm_edges[i] * 0.5 * (geom.rod_r0 * geom.rod_r0);
            }
        }

        for (int i = 0; i < this->n_edges_dof; ++i) {
            int dof = 3 * this->n_nodes_ + i;
            mass_matrix(dof) = edge_mass[i]; 
        }
    }
}

//TODO: check vectorization
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
    std::vector<std::array<int, 2>> local_edges(edges_.begin(), edges_.begin() + n_edges_dof);
    std::vector<double> weights(this->n_edges_dof);
    for (int i = 0; i < this->n_edges_dof; ++i) {
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
/* std::tuple<Eigen::MatrixXd, 
std::vector<std::vector<double>>,
std::vector<std::vector<double>>,
std::vector<std::vector<double>> d>
SoftRobot::_init_curvature_midedge(const Geometry& geo) {
    this->face_edges_ = geo.getFaceEdges();
    this->sign_faces = geo.getSignFaces();

    // Vectorize face processing
    // TODO  std::vector<<Eigen::Matrix3d> all_p_is
    std::vector<std::array<double, 3>> all_p_is = this->extract_block(this->q0, this->face_nodes_shell);
    std::vector<std::array<double, 3>> all_xi_is = this->extract_block(this->q0, this->face_edges_);

    // Compute initial tau0
    this->tau0_ = this->update_pre_comp_shell(this->q0);
    // TODO std::vector<<Eigen::Matrix3d> all_tau0_is
    Eigen::MatrixXd all_tau0_is = this->tau0.block(0, this->face_edges_, this->tau0.rows(), this->face_edges_.size()).transpose();

    // Compute t, f, c for all faces simultaneously
    Eigen::MatrixXd t, f, c;
    compute_tfc_midedge(all_p_is, all_tau0_is, this->sign_faces, t, f, c);

    return std::make_tuple(t, f, c, all_xi_is);
    // TODO std::vector<<Eigen::Matrix3d> t
    //TODO std::vector<<array<double, 3>> f
    // TODO std::vector<<array<double, 3>> c
    // TODO std::vector<<array<double, 3>> xi_is

}
    */

std::vector<Eigen::VectorXd> SoftRobot::update_pre_comp_shell(const Eigen::MatrixXd& q) {
    if (!this->sim_params_.use_mid_edge) {
        return std::vector<Eigen::VectorXd>(); // Return empty vector if not using mid-edge
    }

    // Compute face normals
    std::vector<Eigen::Vector3d> face_normals(this->face_nodes_shell_.size());
    for (size_t i = 0; i < this->face_nodes_shell_.size(); ++i) {
        const auto& face = this->face_nodes_shell_[i];
        Eigen::Vector3d q0 = q.col(face[0]);
        Eigen::Vector3d q1 = q.col(face[1]);
        Eigen::Vector3d q2 = q.col(face[2]);

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
        if (edge_counts[i] > 0) { // Use square brackets
            edge_normals[i] /= edge_counts[i];
        }
    }
    // Compute edge vectors and tau_0
    tau0_.resize(this->edges_.size());
    for (size_t i = 0; i < this->edges_.size(); ++i) {
        int n0 = this->edges_[i][0];
        int n1 = this->edges_[i][1];

        Eigen::Vector3d q0 = q.col(n0);
        Eigen::Vector3d q1 = q.col(n1);
        Eigen::Vector3d edge_vec = q1 - q0;

        Eigen::Vector3d tau_0 = edge_vec.cross(edge_normals[i]).normalized();
        tau0_[i] = tau_0;
    }

    return tau0_;

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
    std::vector<Eigen::Vector3d> a1_transported = parallel_transport(a1_old, tangent0, tangent);

    size_t N = tangent.size();
    std::vector<Eigen::Vector3d> a1(N);
    std::vector<Eigen::Vector3d> a2(N);

    for (size_t i = 0; i < N; ++i) {
        double t_dot = a1_transported[i].dot(tangent[i]);
        a1[i] = a1_transported[i] - t_dot * tangent[i];
        a1[i].normalize();  // Normalize in-place
        a2[i] = tangent[i].cross(a1[i]);
    }

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
/*
Eigen::VectorXd SoftRobot::compute_reference_twist(
    const std::vector<BendTwistSpring>& springs,
    const Eigen::VectorXd& q,
    const Eigen::MatrixXd& a1,
    const Eigen::VectorXd& ref_twist)
{
    if (springs.empty()) {
        return Eigen::VectorXd(); // Empty vector
    }

    Eigen::MatrixXi edges(springs.size(), 2);
    Eigen::VectorXi sgn(springs.size());

    for (size_t i = 0; i < springs.size(); ++i) {
        edges.row(i) = springs[i].edges_ind;
        sgn(i) = springs[i].sgn;
    }

    Eigen::MatrixXd tangent = _compute_tangent(q);
    return compute_reference_twist(edges, sgn, a1, tangent, ref_twist);
} */

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

void SoftRobot::free_edges(const std::vector<std::array<int, 2>>& edges) {
    Eigen::VectorXi edge_dofs = map_edge_to_dof(edges);
    Eigen::VectorXi new_dof = setdiff1d(fixed_dof_, edge_dofs);
    _fix_dof(new_dof);
}

void SoftRobot::fix_edges(const std::vector<std::array<int, 2>>& edges) {
    Eigen::VectorXi edge_dofs = map_edge_to_dof(edges);
    Eigen::VectorXi new_dof = union1d(fixed_dof_, edge_dofs);
    return _fix_dof(new_dof);
}

// TODO FIX
void SoftRobot::_fix_dof(const Eigen::VectorXi& new_fixed_dof) {
    Eigen::VectorXi all_dof = Eigen::VectorXi::LinSpaced(n_dof_, 0, n_dof_ - 1);
    free_dof_ = setdiff1d(all_dof, new_fixed_dof);
    fixed_dof_ = new_fixed_dof;
}


Eigen::VectorXi SoftRobot::_get_intermediate_edge_dof(const std::vector<int>& nodes) {
    std::vector<std::array<int, 2>> matching_edges;
    for (int i = 0; i < this->n_edges_; ++i) {
        // Extract node indices for the current edge
        int node1_index = this->edges_[i][0];
        int node2_index = this->edges_[i][1];

        // Check if both nodes of the edge are in the 'nodes' vector
        bool node1_found = false;
        bool node2_found = false;
        for (int j = 0; j < nodes.size(); ++j) {
            if (node1_index == nodes[j]) node1_found = true;
            if (node2_index == nodes[j]) node2_found = true;
        }

        if (node1_found && node2_found) {
            matching_edges.push_back(this->edges_[i]);
        }
    }

    std::vector<std::array<int, 2>> edges_vec = matching_edges;

    // map_edge_to_dof expects edge *pairs*, not edge *indices*
    return map_edge_to_dof(edges_vec);
}

// Utility Helpers
Eigen::VectorXi setdiff1d(const Eigen::VectorXi& a, const Eigen::VectorXi& b) {
    std::unordered_set<int> b_set(b.data(), b.data() + b.size());
    std::vector<int> result;
    for (int i = 0; i < a.size(); ++i)
        if (b_set.find(a(i)) == b_set.end())
            result.push_back(a(i));
    return Eigen::Map<Eigen::VectorXi>(result.data(), result.size());
}

Eigen::VectorXi union1d(const Eigen::VectorXi& a, const Eigen::VectorXi& b) {
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

void SoftRobot::twist_edges(const std::vector<std::array<int, 2>>& edges, const Eigen::VectorXd& perturbation) {
    Eigen::VectorXi dof_mask = map_edge_to_dof(edges);
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
std::vector<std::array<int, 3>> SoftRobot::map_node_to_dof(const std::vector<int>& node_indices) const {
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

Eigen::VectorXi SoftRobot::map_edge_to_dof(const std::vector<std::array<int, 2>>& edges) const {
    int num_edges = edges.size();
    Eigen::VectorXi dof_indices(num_edges);
    int offset = 3 * this->n_nodes_;

    for (int i = 0; i < num_edges; ++i) {
        dof_indices(i) = offset + i; // Assuming edges are indexed sequentially
    }

    return dof_indices;
}

// Maps face-edge indices to global DOF indices
Eigen::VectorXi SoftRobot::map_face_edge_to_dof(const std::vector<int>& edge_nums) const {
    int num_edges = edge_nums.size();
    Eigen::VectorXi dof_indices(num_edges);
    int offset = 3 * n_nodes_ + n_edges_dof;

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

// std::shared_ptr<SoftRobot> SoftRobot::update(
//     std::optional<Eigen::VectorXd> q,
//     std::optional<Eigen::VectorXd> u,
//     std::optional<Eigen::VectorXd> a,
//     std::optional<Eigen::VectorXd> a1,
//     std::optional<Eigen::VectorXd> a2,
//     std::optional<Eigen::VectorXd> m1,
//     std::optional<Eigen::VectorXd> m2,
//     std::optional<Eigen::VectorXd> ref_twist,
//     std::optional<Eigen::VectorXi> free_dof_
// ) const {
//     RobotState new_state = state_; // copy

//     if (q) new_state.q = *q;
//     if (u) new_state.u = *u;
//     if (a) new_state.a = *a;
//     if (a1) new_state.a1 = *a1;
//     if (a2) new_state.a2 = *a2;
//     if (m1) new_state.m1 = *m1;
//     if (m2) new_state.m2 = *m2;
//     if (ref_twist) new_state.ref_twist = *ref_twist;
//     if (free_dof_) new_state.free_dof_ = *free_dof_;

//     return std::make_unique<SoftRobot>(new_state);
// }


