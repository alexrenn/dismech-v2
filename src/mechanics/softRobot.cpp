#include "softRobot.h"
#include "robotState.h"
#include "geometry.h"

void SoftRobot::SoftRobot(const GeomParams& geom, const Material& material, const Geometry& geo, 
                            const SimParams& sim_params, const Environment& env, const RobotState& state)
        : sim_params_(sim_params), env_(env), state_(state){

        this->_init_geometry(geo);
        this->_init_stiffness(geom, material);
        this->_init_state(geo);
        this->_init_springs(geo);
        mass_matrix = get_mass_matrix(geom, material);
    }

void SoftRobot::_init_geometry(const Geometry& geo)
 {
    // Initialize Geometry Properties
    n_nodes_ = geo.nodes.size();
    n_edges_rod_only_ = geo.rod_edges.size();
    n_edges_shell_only_ = geo.shell_edges.size();
    int n_edges_joint = geo.rod_shell_joint_edges_total.size();
    n_edges_ = geo.edges.size();
    n_edges_dof_ = n_edges_rod_only + n_edges_joint;
    n_faces_ = geo.face_nodes.size(); // TODO come back to it

    // Assign references to the original data, avoiding duplicates
    nodes_ = geo.nodes;
    edges_ = geo.edges;
    face_nodes_shell_ = geo.face_nodes;
    twist_angles_ = geo.twist_angles;

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

    // Midedge bending has more DOF
    if (sim_params_.use_mid_edges) {
        n_dof_ += n_edges_shell_only_;
    
        Eigen::VectorXd extra = Eigen::VectorXd::Zero(n_edges_shell_only_); // Create extra zeros
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
    ref_len = _get_ref_len();
    voronoi_ref_len = _get_voronoi_ref_len();
    face_area = get_face_area();
 }

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

void SoftRobot::_init_state(const Geometry& geo) {
    // Initialize RobotState state for q0
    auto [a1, a2] = this->compute_space_parallel();

    auto [m1, m2] = this->compute_material_directors(this->q0, a1, a2);

    int n = static_cast<int>(geo.bend_twist_springs.size());
    for (int i = 0; i < n; ++i) {
        edges_[i][0] = geo.bend_twist_springs[i][1];
        edges_[i][1] = geo.bend_twist_springs[i][3];
    }

    std::vector<int> sign(n, 0);
    for (int i = 0; i < n; ++i) {
        sign[i] = bend_twist_signs_[i];  // Assuming you have bend_twist_signs_
    }

    std::vector<double> undef_ref_twist; // RADHA is this the correct type?
    if (edges.size() > 0) {
        Eigen::VectorXd zero_twist = Eigen::VectorXd::Zero(sign.size());
        undef_ref_twist = compute_reference_twist(
            edges, sign, a1, this->tangent, zero_twist);
    } else {
        undef_ref_twist.resize(0, 0); // Empty matrix
    }

    // Initialize robot state
    this->state = RobotState::init(
        this->q0, a1, a2, m1, m2, undef_ref_twist, this->tau0);
}

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
                geo.face_edges[i],
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
}

void SoftRobot::_get_mass_matrix(const GeomParams& geom, const Material& material) {
    mass_matrix = Eigen::VectorXd::Zero(this->n_dof);  // Reuse the member variable

    // Shell face contributions
    if (this->n_faces > 0) {
        std::vector<std::array<int, 3>>faces = this->face_nodes_shell; // (n_faces, 3)

        std::vector<double> areas(this->n_faces_); // TODO std vector of doubles (area of each triangle)
        std::vector<double> m_shell(this->n_faces_);; // TODO std vector of doubles (area of each triangle)

        // TODO combine into one for loop
        for (size_t i = 0; i < faces.size(); ++i) {
            const auto& face = faces[i];
            Eigen::Vector3d v1 = this->nodes_[face[1]] - this->nodes_[face[0]];
            Eigen::Vector3d v2 = this->nodes_[face[2]] - this->nodes_[face[1]];

            double area = 0.5 * v1.cross(v2).norm();
            areas[i] = area;
            m_shell(i) = material.density * area * geom.shell_h;
        }
        for (size_t i = 0; i < faces.size(); ++i) {
            const auto& face = faces[i];
            for (int j = 0; j < 3; ++j) {
                int node_id = face[j];
                for (int k = 0; k < 3; ++k) {
                    int dof = 3 * node_id + k;
                    mass_matrix(dof) += m_shell(i) / 3.0;
                }
            }
        }
    }

    // Node contributions
    if (this->n_nodes > 0) {
        std::vector<double> dm_nodes; // RADHA types?
        if (geom.axs >= 0.0) {
            dm_nodes = this->voronoi_ref_len * geom.axs * material.density; // TODO change
        } else {
            double area = M_PI * geom.rod_r0 * geom.rod_r0;
            dm_nodes = this->voronoi_ref_len * area * material.density;
        }

        // TODO: optimize??
        for (int i = 0; i < this->n_nodes; ++i) {
            for (int j = 0; j < 3; ++j) {
                int dof = 3 * i + j;
                mass_matrix(dof) += dm_nodes(i);
            }
        }
    }

    // Edge contributions
    // TODO match python code
    if (this->n_edges_dof > 0) {
        std::vector<double> dm_edges; // RADHA types?
        if (geom.axs >= 0.0) {
            dm_edges = this->ref_len.head(this->n_edges_dof) * geom.axs * material.density; // TODO change
        } else {
            double area = M_PI * geom.rod_r0 * geom.rod_r0;
            dm_edges = this->ref_len.head(this->n_edges_dof) * area * material.density; // head is first n entries
        }

        std::vector<double> edge_mass = dm_edges * 0.5 * (geom.rod_r0 * geom.rod_r0);
        for (int i = 0; i < this->n_edges_dof; ++i) {
            int dof = 3 * this->n_nodes + i;
            mass_matrix(dof) = edge_mass(i); 
        }
    }
}

//TODO: check vectorization
void SoftRobot::_get_ref_len() const {
    std::vector<Eigen::VectorXd> node1 = nodes(edges.col(0), Eigen::all);  // Start points
    std::vector<Eigen::VectorXd> node2 = nodes(edges.col(1), Eigen::all);  // End points
    std::vector<Eigen::VectorXd> edge_vectors = node2 - node1;
    this->ref_len = edge_vectors.rowwise().norm();  // Lengths of vectors
    return ref_len;
}

void SoftRobot::_get_voronoi_ref_len() const {
    Eigen::MatrixXi edges = this->edges.topRows(this->n_edges_dof); // Select first n_edges_dof edges // TODO std vector array int 2
    int n_nodes = this->n_nodes;
    Eigen::VectorXd weights = 0.5 * this->ref_len.head(this->n_edges_dof); // Half edge lengths // TODO std vector double

    Eigen::VectorXd contributions = Eigen::VectorXd::Zero(n_nodes); // Initialize with zeros // TODO std vector double

    for (int i = 0; i < edges.rows(); ++i) {
        contributions(edges(i, 0)) += weights(i); // Add half length to first node
        contributions(edges(i, 1)) += weights(i); // Add half length to second node
    }

    this->voronoi_ref_len = contributions;
}

void SoftRobot::_get_voronoi_area() const {
    if (this->face_nodes_shell.size() == 0) {
        return Eigen::VectorXd(); // Return empty vector if no faces
        // TODO std vector double
    }

    const Eigen::MatrixXi& faces = this->face_nodes_shell;
    Eigen::MatrixXd v1 = this->nodes(faces.col(1), Eigen::all) - this->nodes(faces.col(0), Eigen::all);
    Eigen::MatrixXd v2 = this->nodes(faces.col(2), Eigen::all) - this->nodes(faces.col(1), Eigen::all);
    Eigen::MatrixXd cross = v1.rowwise().cross(v2);
    Eigen::VectorXd areas = 0.5 * cross.rowwise().norm(); // TODO std vector double

    Eigen::VectorXd node_areas = Eigen::VectorXd::Zero(this->n_nodes); // TODO std vector double
    for (int i = 0; i < faces.rows(); ++i) {
        node_areas(faces(i, 0)) += areas(i) / 3.0;
        node_areas(faces(i, 1)) += areas(i) / 3.0;
        node_areas(faces(i, 2)) += areas(i) / 3.0;
    }

    this->voronoi_area = node_areas
}

void SoftRobot::_get_face_area() const {
    if (this->n_faces == 0) {
        return Eigen::VectorXd(); // Empty vector if no faces // TODO std vector double
    }

    const Eigen::MatrixXi& faces = this->face_nodes_shell; // TODO std vector array int 3
    Eigen::VectorXd areas(faces.rows()); // TODO std vector double

    for (int i = 0; i < faces.rows(); ++i) {
        Eigen::Vector3d a = this->nodes.row(faces(i, 0));
        Eigen::Vector3d b = this->nodes.row(faces(i, 1));
        Eigen::Vector3d c = this->nodes.row(faces(i, 2));

        Eigen::Vector3d v1 = b - a;
        Eigen::Vector3d v2 = c - b;
        double area = 0.5 * (v1.cross(v2)).norm();

        areas(i) = area;
    }

    this->face_area = area;
}

void SoftRobot::scale_mass_matrix(const Eigen::Vector3d& nodes, double scale) {
    Eigen::VectorXi dof_indices = this->map_node_to_dof(nodes); // TODO std vector int
    for (int i = 0; i < dof_indices.size(); ++i) {
        this->mass_matrix(dof_indices(i)) *= scale;
    }
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> SoftRobot::_init_curvature_midedge(const Geometry& geo) {
    this->face_edges = geo.face_edges;
    this->sign_faces = geo.sign_faces;

    // Vectorize face processing
    Eigen::MatrixXd all_p_is = this->q0.block(0, this->map_node_to_dof(this->face_nodes_shell), this->q0.rows(), this->face_nodes_shell.size());
    Eigen::MatrixXd all_xi_is = this->q0.block(0, this->map_edge_to_dof(this->face_edges), this->q0.rows(), this->face_edges.size());

    // Compute initial tau0
    this->tau0 = this->update_pre_comp_shell(this->q0);
    Eigen::MatrixXd all_tau0_is = this->tau0.block(0, this->face_edges, this->tau0.rows(), this->face_edges.size()).transpose();

    // Compute t, f, c for all faces simultaneously
    Eigen::MatrixXd t, f, c;
    compute_tfc_midedge(all_p_is, all_tau0_is, this->sign_faces, t, f, c);

    return std::make_tuple(t, f, c, all_xi_is);
}

Eigen::MatrixXd SoftRobot::update_pre_comp_shell(const Eigen::MatrixXd& q) {
    if (!this->sim_params.use_mid_edge) {
        return Eigen::MatrixXd::Zero(0, 0);
    }

    // Compute face normals

    // TODO std vector EigenVecXd (for all 3)
    Eigen::MatrixXd v1 = q.block(0, this->map_node_to_dof(this->face_nodes_shell.col(1)), q.rows(), this->face_nodes_shell.cols()) -
                         q.block(0, this->map_node_to_dof(this->face_nodes_shell.col(0)), q.rows(), this->face_nodes_shell.cols());
    Eigen::MatrixXd v2 = q.block(0, this->map_node_to_dof(this->face_nodes_shell.col(2)), q.rows(), this->face_nodes_shell.cols()) -
                         q.block(0, this->map_node_to_dof(this->face_nodes_shell.col(1)), q.rows(), this->face_nodes_shell.cols());
    Eigen::MatrixXd face_normals = v1.cross(v2).rowwise().normalized();

    // Accumulate edge normals
        // TODO std vector EigenVecXd
    Eigen::MatrixXd edge_normals = Eigen::MatrixXd::Zero(this->n_edges, 3);
    for (int i = 0; i < this->face_edges.rows(); ++i) {
        edge_normals.row(this->face_edges(i, 0)) += face_normals.row(i);
        edge_normals.row(this->face_edges(i, 1)) += face_normals.row(i);
    }

    // Normalize edge normals
    // combine these with the previous loop
    Eigen::VectorXd edge_counts = Eigen::VectorXd::Zero(this->n_edges);
    for (int i = 0; i < this->face_edges.rows(); ++i) {
        edge_counts(this->face_edges(i, 0)) += 1;
        edge_counts(this->face_edges(i, 1)) += 1;
    }
    for (int i = 0; i < this->n_edges; ++i) {
        if (edge_counts(i) > 0) {
            edge_normals.row(i) /= edge_counts(i);
        }
    }

    // Compute edge vectors and tau_0
    // TODO std vector EigenVecXd
    Eigen::MatrixXd edge_vecs = q.block(0, this->map_node_to_dof(this->edges.col(1)), q.rows(), this->edges.cols()) -
                                q.block(0, this->map_node_to_dof(this->edges.col(0)), q.rows(), this->edges.cols());
    // TODO std vector EigenVecXd
    Eigen::MatrixXd tau_0 = edge_vecs.cross(edge_normals);
    tau_0 = tau_0.rowwise().normalized();

    return tau_0.transpose();
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> SoftRobot::compute_space_parallel() {
    Eigen::MatrixXd a1 = Eigen::MatrixXd::Zero(this->n_edges_dof, 3);
    Eigen::MatrixXd a2 = Eigen::MatrixXd::Zero(this->n_edges_dof, 3);

    this->__tangent = this->_compute_tangent(this->__q0);

    if (this->__tangent.size() > 0) {
        // Initialize first a1
        Eigen::VectorXd a1_init = this->__tangent.row(0).cross(Eigen::Vector3d(0, 1, 0));
        if (a1_init.norm() < 1e-6) {
            a1_init = this->__tangent.row(0).cross(Eigen::Vector3d(0, 0, -1));
        }
        a1.row(0) = a1_init.normalized();
        a2.row(0) = this->__tangent.row(0).cross(a1.row(0));

        // Iterative parallel transport (depends on previous a1)
        for (int i = 1; i < this->n_edges_dof; ++i) {
            Eigen::Vector3d t_prev = this->__tangent.row(i - 1);
            Eigen::Vector3d t_curr = this->__tangent.row(i);
            Eigen::Vector3d a1_prev = a1.row(i - 1);

            a1.row(i) = parallel_transport(a1_prev, t_prev, t_curr);
            a1.row(i) -= a1.row(i).dot(t_curr) * t_curr;
            a1.row(i) = a1.row(i).normalized();
            a2.row(i) = t_curr.cross(a1.row(i));
        }
    }

    return std::make_tuple(a1, a2);
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> SoftRobot::compute_time_parallel(
    const Eigen::MatrixXd& a1_old,
    const Eigen::MatrixXd& q0,
    const Eigen::MatrixXd& q) 
{
    Eigen::MatrixXd tangent0 = this->_compute_tangent(q0);
    Eigen::MatrixXd tangent = this->_compute_tangent(q);

    Eigen::MatrixXd a1_transported = parallel_transport(a1_old, tangent0, tangent);

    // Orthonormalization
    Eigen::VectorXd t_dot = (a1_transported.array() * tangent.array()).rowwise().sum();
    Eigen::MatrixXd a1 = a1_transported - tangent.array().colwise() * t_dot.array();
    
    a1 = a1.array().rowwise() / a1.rowwise().norm().array(); // Normalize each row
    Eigen::MatrixXd a2(a1.rows(), 3);
    for (int i = 0; i < a1.rows(); ++i) {
        a2.row(i) = tangent.row(i).cross(a1.row(i));
    }

    return std::make_tuple(a1, a2);
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> SoftRobot::compute_material_directors(
    const Eigen::VectorXd& q,
    const Eigen::MatrixXd& a1,
    const Eigen::MatrixXd& a2)
{
    Eigen::VectorXd theta = q.segment(3 * n_nodes, n_edges_dof);
    Eigen::MatrixXd cos_theta = theta.array().cos().matrix();
    Eigen::MatrixXd sin_theta = theta.array().sin().matrix();

    Eigen::MatrixXd m1 = cos_theta.asDiagonal() * a1 + sin_theta.asDiagonal() * a2;
    Eigen::MatrixXd m2 = -sin_theta.asDiagonal() * a1 + cos_theta.asDiagonal() * a2;

    return std::make_pair(m1, m2);
}

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
    return ::compute_reference_twist(edges, sgn, a1, tangent, ref_twist);
}

Eigen::MatrixXd SoftRobot::_compute_tangent(const Eigen::VectorXd& q)
{
    Eigen::MatrixXi edges = this->edges.topRows(n_edges_dof);
    Eigen::VectorXi n0 = edges.col(0);
    Eigen::VectorXi n1 = edges.col(1);

    Eigen::MatrixXd pos0(n_edges_dof, 3), pos1(n_edges_dof, 3);
    for (int i = 0; i < n_edges_dof; ++i) {
        for (int j = 0; j < 3; ++j) {
            pos0(i, j) = q(map_node_to_dof(n0(i)) + j);
            pos1(i, j) = q(map_node_to_dof(n1(i)) + j);
        }
    }

    Eigen::MatrixXd vecs = pos1 - pos0;
    Eigen::VectorXd norms = vecs.rowwise().norm();
    for (int i = 0; i < norms.size(); ++i)
        if (norms(i) < 1e-10) norms(i) = 1.0;

    Eigen::MatrixXd tangent = vecs.array().colwise() / norms.array();

    for (int i = 0; i < tangent.rows(); ++i)
        for (int j = 0; j < tangent.cols(); ++j)
            if (std::abs(tangent(i, j)) < _TANGENT_THRESHOLD)
                tangent(i, j) = 0;

    return tangent;
}

// Fix/free nodes and edges

std::shared_ptr<SoftRobot> SoftRobot::free_nodes(const std::vector<int>& nodes, std::optional<int> axis, bool fix_edges) const {
    Eigen::VectorXi node_dof_mask = _get_node_dof_mask(nodes, axis);
    Eigen::VectorXi new_dof = setdiff1d(fixed_dof, node_dof_mask);

    if (fix_edges) {
        Eigen::VectorXi edge_dof_mask = _get_intermediate_edge_dof(nodes);
        new_dof = setdiff1d(new_dof, edge_dof_mask);
    }

    return std::make_shared<SoftRobot>(_fix_dof(new_dof));
}

std::shared_ptr<SoftRobot> SoftRobot::fix_nodes(const std::vector<int>& nodes, std::optional<int> axis, bool fix_edges) const {
    Eigen::VectorXi node_dof_mask = _get_node_dof_mask(nodes, axis);
    Eigen::VectorXi new_dof = union1d(fixed_dof, node_dof_mask);

    if (fix_edges) {
        Eigen::VectorXi edge_dof_mask = _get_intermediate_edge_dof(nodes);
        new_dof = union1d(new_dof, edge_dof_mask);
    }

    return std::make_shared<SoftRobot>(_fix_dof(new_dof));
}

std::shared_ptr<SoftRobot> SoftRobot::free_edges(const std::vector<int>& edges) const {
    Eigen::VectorXi edge_dofs = map_edge_to_dof(edges);
    Eigen::VectorXi new_dof = setdiff1d(fixed_dof, edge_dofs);
    return std::make_shared<SoftRobot>(_fix_dof(new_dof));
}

std::shared_ptr<SoftRobot> SoftRobot::fix_edges(const std::vector<int>& edges) const {
    Eigen::VectorXi edge_dofs = map_edge_to_dof(edges);
    Eigen::VectorXi new_dof = union1d(fixed_dof, edge_dofs);
    return std::make_shared<SoftRobot>(_fix_dof(new_dof));
}

std::shared_ptr<SoftRobot> SoftRobot::_fix_dof(const Eigen::VectorXi& new_fixed_dof) const {
    SoftRobot updated = *this; // shallow copy
    Eigen::VectorXi all_dof = Eigen::VectorXi::LinSpaced(n_dof, 0, n_dof - 1);
    Eigen::VectorXi free_dof = setdiff1d(all_dof, new_fixed_dof);
    updated.update(free_dof);
    return updated;
}


Eigen::VectorXi SoftRobot::_get_intermediate_edge_dof(const Eigen::VectorXi& nodes) const {
    Eigen::VectorXi edge_mask = (
        isin(edges.col(0), nodes) &&
        isin(edges.col(1), nodes)
    );

    std::vector<int> matching_edges;
    for (int i = 0; i < edge_mask.size(); ++i) {
        if (edge_mask(i)) matching_edges.push_back(i);
    }

    Eigen::VectorXi edges_vec = Eigen::Map<Eigen::VectorXi>(matching_edges.data(), matching_edges.size());

    if (sim_params.two_d_sim) {
        Eigen::VectorXi all_edges = Eigen::VectorXi::LinSpaced(n_edges_dof, 0, n_edges_dof - 1);
        edges_vec = union1d(edges_vec, all_edges);
    }

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
std::shared_ptr<SoftRobot> SoftRobot::move_nodes(const std::vector<int>& nodes, const Eigen::VectorXd& perturbation, std::optional<int> axis) const {
    Eigen::VectorXd q = state.q; // Copy current state
    Eigen::VectorXi dof_mask = _get_node_dof_mask(nodes, axis);
    q(dof_mask) += perturbation;
    return std::make_shared<SoftRobot> (update(q));
}

std::shared_ptr<SoftRobot> SoftRobot::twist_edges(const std::vector<int>& edges, const Eigen::VectorXd& perturbation) const {
    Eigen::VectorXd q = state.q; // Copy current state
    Eigen::VectorXi dof_mask = map_edge_to_dof(edges);
    q(dof_mask) += perturbation;
    return std::make_shared<SoftRobot> (update(q));
}

//Utility
static Eigen::VectorXi get_node_dof_mask(const std::vector<int>& nodes, std::optional<int> axis = std::nullopt) {
    Eigen::MatrixXi node_dof = map_node_to_dof(nodes);
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
static Eigen::MatrixXi map_node_to_dof(const std::vector<int>& node_nums) {
    int num_nodes = node_nums.size();
    Eigen::MatrixXi dof_map(num_nodes, 3);

    for (int i = 0; i < num_nodes; ++i) {
        int base = 3 * node_nums[i];
        dof_map(i, 0) = base;
        dof_map(i, 1) = base + 1;
        dof_map(i, 2) = base + 2;
    }
    return dof_map;
}

Eigen::VectorXi SoftRobot::map_edge_to_dof(const std::vector<int>& edge_nums) const {
    Eigen::VectorXi dof_indices(edge_nums.size());
    int offset = 3 * this->n_nodes;

    for (int i = 0; i < edge_nums.size(); ++i) {
        dof_indices(i) = offset + edge_nums[i];
    }

    return dof_indices;
}

// Maps face-edge indices to global DOF indices
Eigen::VectorXi SoftRobot::map_face_edge_to_dof(const std::vector<int>& edge_nums) const {
    int num_edges = edge_nums.size();
    Eigen::VectorXi dof_indices(num_edges);
    int offset = 3 * n_nodes + n_edges_dof;

    for (int i = 0; i < num_edges; ++i) {
        dof_indices(i) = offset + edge_nums[i];
    }

    return dof_indices;
}

std::shared_ptr<SoftRobot> SoftRobot::update(
    std::optional<Eigen::VectorXd> q,
    std::optional<Eigen::VectorXd> u,
    std::optional<Eigen::VectorXd> a,
    std::optional<Eigen::VectorXd> a1,
    std::optional<Eigen::VectorXd> a2,
    std::optional<Eigen::VectorXd> m1,
    std::optional<Eigen::VectorXd> m2,
    std::optional<Eigen::VectorXd> ref_twist,
    std::optional<Eigen::VectorXi> free_dof
) const {
    RobotState new_state = state_; // copy

    if (q) new_state.q = *q;
    if (u) new_state.u = *u;
    if (a) new_state.a = *a;
    if (a1) new_state.a1 = *a1;
    if (a2) new_state.a2 = *a2;
    if (m1) new_state.m1 = *m1;
    if (m2) new_state.m2 = *m2;
    if (ref_twist) new_state.ref_twist = *ref_twist;
    if (free_dof) new_state.free_dof = *free_dof;

    return std::make_unique<SoftRobot>(new_state);
}


