#include "bend_energy.h"

#include <Eigen/Dense>
#include <numeric>

// Batch outer product: a, b are N×3 → returns N×3×3
std::vector<Eigen::Matrix3d> batch_outer(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b) {
    const int N = a.rows();
    std::vector<Eigen::Matrix3d> result(N);
    for (int i = 0; i < N; ++i) {
        result[i] = a.row(i).transpose() * b.row(i);
    }
    return result;
}

// Cross product matrix for each row of v: N×3 → returns N×3×3
std::vector<Eigen::Matrix3d> cross_mat_batch(const Eigen::MatrixXd& v) {
    const int N = v.rows();
    std::vector<Eigen::Matrix3d> result(N);
    for (int i = 0; i < N; ++i) {
        const auto& vi = v.row(i);
        Eigen::Matrix3d m;
        m <<     0, -vi.z(),  vi.y(),
              vi.z(),     0, -vi.x(),
             -vi.y(),  vi.x(),     0;
        result[i] = m;
    }
    return result;
}

// Helper: outer product of a row vector with itself
std::vector<Eigen::Matrix3d> batch_outer_self(const Eigen::MatrixXd& a) {
    return batch_outer(a, a);
}

// Helper: identity matrix broadcast
std::vector<Eigen::Matrix3d> broadcast_identity(int N) {
    std::vector<Eigen::Matrix3d> out(N, Eigen::Matrix3d::Identity());
    return out;
}

void assign_block(Eigen::MatrixXd& H, int i, int j, const Eigen::Matrix3d& block) {
    H.block<3,3>(i,j) = block;
}

void assign_twist_blocks(Eigen::MatrixXd& H, int i, int j, const Eigen::Vector3d& v) {
    H.block<3,1>(i,j) = v;
    H.block<1,3>(j,i) = v.transpose();
}


// Per-spring second derivative tensors and Hessian construction
void assign_second_derivatives(
    const Eigen::MatrixXd& D2De2,
    const Eigen::MatrixXd& D2DeDf,
    const Eigen::MatrixXd& D2DfDe,
    const Eigen::MatrixXd& D2Df2,
    double D2t1, double D2t2,
    const Eigen::MatrixXd& D2ct0_e, const Eigen::MatrixXd& D2ct1_e,
    const Eigen::MatrixXd& D2ct0_f, const Eigen::MatrixXd& D2ct1_f,
    Eigen::MatrixXd& H) {

    // Position-position blocks
    H.block(0, 0, 3, 3) = D2De2;
    H.block(0, 3, 3, 3) = -D2De2 + D2DeDf;
    H.block(0, 6, 3, 3) = -D2DeDf;

    H.block(3, 0, 3, 3) = -D2De2 + D2DfDe;
    H.block(3, 3, 3, 3) = D2De2 - D2DeDf - D2DfDe + D2Df2;
    H.block(3, 6, 3, 3) = D2DeDf - D2Df2;

    H.block(6, 0, 3, 3) = -D2DfDe;
    H.block(6, 3, 3, 3) = D2DfDe - D2Df2;
    H.block(6, 6, 3, 3) = D2Df2;

    // Twist scalar terms
    H(9, 9) = D2t1;
    H(10, 10) = D2t2;

    // Position-twist coupling
    H.block(0, 9, 3, 1) = -D2ct0_e;
    H.block(3, 9, 3, 1) = D2ct0_e - D2ct1_e;
    H.block(6, 9, 3, 1) = D2ct1_e;

    H.block(0, 10, 3, 1) = -D2ct0_f;
    H.block(3, 10, 3, 1) = D2ct0_f - D2ct1_f;
    H.block(6, 10, 3, 1) = D2ct1_f;

    // Twist-position transpose
    H.block(9, 0, 1, 3) = -D2ct0_e.transpose();
    H.block(9, 3, 1, 3) = (D2ct0_e - D2ct1_e).transpose();
    H.block(9, 6, 1, 3) = D2ct1_e.transpose();

    H.block(10, 0, 1, 3) = -D2ct0_f.transpose();
    H.block(10, 3, 1, 3) = (D2ct0_f - D2ct1_f).transpose();
    H.block(10, 6, 1, 3) = D2ct1_f.transpose();
}

// Constructor

BendEnergy::BendEnergy(const std::vector<BendTwistSpring>& springs,
                       const RobotState& initial_state,
                       std::function<Eigen::MatrixXd(const Eigen::MatrixXd&)> get_strain)
    : ElasticEnergy(
        // Convert stiff_EI and voronoi_len to K matrix (Nx2)
        [&]() {
            const int N = static_cast<int>(springs.size());
            Eigen::MatrixXd K(N, 2);
            for (int i = 0; i < N; ++i) {
                K(i, 0) = springs[i].stiff_EI[0] / springs[i].voronoi_len;
                K(i, 1) = springs[i].stiff_EI[1] / springs[i].voronoi_len;
            }
            return K;
        }(),
        // nodes_ind
        [&]() {
            std::vector<std::array<int, 3>> out;
            for (const auto& s : springs) out.push_back(s.nodes_ind);
            return out;
        }(),
        // ind
        [&]() {
            std::vector<std::vector<int>> out;
            for (const auto& s : springs) out.push_back(s.ind);
            return out;
        }(),
        initial_state,
        get_strain
    )
{

    // Extract _edges_ind and _sgn
    N = springs.size();
    std::cout << "N: " << N << std::endl;
    _edges_ind.resize(N, 2);
    _sgn.resize(N, 2);

    for (int i = 0; i < N; ++i) {
        _edges_ind(i, 0) = springs[i].edges_ind[0];
        _edges_ind(i, 1) = springs[i].edges_ind[1];
        _sgn(i, 0) = springs[i].sgn[0];
        _sgn(i, 1) = springs[i].sgn[1];
    }

    // Initialize _sign_grad (N x 11)
    _sign_grad = Eigen::MatrixXd::Ones(N, 11);
    _sign_grad.col(9) = _sgn.col(0).cast<double>();
    _sign_grad.col(10) = _sgn.col(1).cast<double>();

    // Initialize _sign_hess: outer product of sign_grad with itself (N x 11 x 11)
    _sign_hess.reserve(N);
    for (int i = 0; i < N; ++i) {
        Eigen::MatrixXd outer = _sign_grad.row(i).transpose() * _sign_grad.row(i);
        _sign_hess.push_back(outer);
    }
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
 BendEnergy::_get_adjusted_material_directors(const std::vector<Eigen::Vector3d>& m1,
                                                  const std::vector<Eigen::Vector3d>& m2) const
{
    Eigen::MatrixXd m1e(N, 3), m2e(N, 3), m1f(N, 3), m2f(N, 3);

    for (int i = 0; i < N; ++i) {
        int e0 = _edges_ind(i, 0);
        int e1 = _edges_ind(i, 1);
        int s0 = _sgn(i, 0);
        int s1 = _sgn(i, 1);

        m1e.row(i) = m1[e0].transpose(); // (1,3)
        m2e.row(i) = (m2[e0] * static_cast<double>(s0)).transpose(); // (1,3)
        m1f.row(i) = m1[e1].transpose(); // (1,3)
        m2f.row(i) = (m2[e1] * static_cast<double>(s1)).transpose(); // (1,3)
    }
    return std::make_tuple(m1e, m2e, m1f, m2f);
}

Eigen::MatrixXd BendEnergy::get_strain(const RobotState& state) const{
    // 1. Get node positions (n0p, n1p, n2p)
    auto node_pos = this->get_node_positions(state.q0); // shape: M x N x 3
    int M = node_pos.size();
// n0p, n1p, n2p = self._get_node_pos(state.q)
    Eigen::MatrixXd n0p(N, 3), n1p(N, 3), n2p(N, 3);
    for (int i = 0; i < N; ++i) {
        int n0 = _node_dof_ind[3*i + 0];
        int n1 = _node_dof_ind[3*i + 1];
        int n2 = _node_dof_ind[3*i + 2];
        n0p.row(i) = node_pos[0][n0];
        n1p.row(i) = node_pos[0][n1];
        n2p.row(i) = node_pos[0][n2];
    }
    // 2. Get adjusted material directors
    auto [m1e, m2e, m1f, m2f] = _get_adjusted_material_directors(state.m1, state.m2);

    // 3. Precompute edge tangents
    Eigen::MatrixXd ee = n1p - n0p;
    Eigen::MatrixXd ef = n2p - n1p;

    Eigen::VectorXd norm_e = ee.rowwise().norm();
    Eigen::VectorXd norm_f = ef.rowwise().norm();
    std::cout << "norm_e: " << norm_e.transpose() << std::endl;
    std::cout << "norm_f: " << norm_f.transpose() << std::endl;

    Eigen::MatrixXd te = ee.array().colwise() / norm_e.array();
    Eigen::MatrixXd tf = ef.array().colwise() / norm_f.array();

    std::cout << "te row 0: " << te.row(0) << std::endl;
    std::cout << "tf row 0: " << tf.row(0) << std::endl;
    std::cout << "dot: " << te.row(0).dot(tf.row(0)) << std::endl;

    // 4. chi and kb
    Eigen::VectorXd chi = (te.array() * tf.array()).rowwise().sum() + 1.0;
    Eigen::VectorXd chi_inv = chi.cwiseInverse();
    std::cout << "chi: " << chi.transpose() << std::endl;

    Eigen::MatrixXd cross_te_tf(N, 3);
    for (int i = 0; i < N; ++i) {
        cross_te_tf.row(i) = Eigen::Vector3d(te.row(i)).cross(Eigen::Vector3d(tf.row(i)));    
    }

    Eigen::MatrixXd kb = (2.0 * cross_te_tf).array() * chi_inv.replicate(1, 3).array(); // Broadcast explicity with .replicate() and .array() to match shapes

    std::cout << "kb.shape: " << kb.rows() << "x" << kb.cols() << "\n";
    std::cout << "m2e.shape: " << m2e.rows() << "x" << m2e.cols() << "\n";
    std::cout << "m2f.shape: " << m2f.rows() << "x" << m2f.cols() << "\n";
    // 5. Compute curvatures
    Eigen::VectorXd kappa1 = 0.5 * ((kb.array() * (m2e + m2f).array()).rowwise().sum());
    Eigen::VectorXd kappa2 = -0.5 * ((kb.array() * (m1e + m1f).array()).rowwise().sum());

    // 6. Stack into N x 2 strain matrix
    Eigen::MatrixXd strain(N, 2);
    strain.col(0) = kappa1;
    strain.col(1) = kappa2;

    return strain;
}

void batch_assign_blocks(std::vector<Eigen::MatrixXd>& H_out,
                         const std::vector<Eigen::Matrix3d>& D2De2,
                         const std::vector<Eigen::Matrix3d>& D2DeDf,
                         const std::vector<Eigen::Matrix3d>& D2DfDe,
                         const std::vector<Eigen::Matrix3d>& D2Df2,
                         const Eigen::VectorXd& D2t1,
                         const Eigen::VectorXd& D2t2,
                         const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& D2ct) {

    const int N = static_cast<int>(H_out.size());
    for (int i = 0; i < N; ++i) {
        Eigen::MatrixXd& H = H_out[i];

        // Positional block 3x3 sections
        assign_block(H, 0, 0, D2De2[i]);
        assign_block(H, 0, 3, -D2De2[i] + D2DeDf[i]);
        assign_block(H, 0, 6, -D2DeDf[i]);

        assign_block(H, 3, 0, -D2De2[i] + D2DfDe[i]);
        assign_block(H, 3, 3, D2De2[i] - D2DeDf[i] - D2DfDe[i] + D2Df2[i]);
        assign_block(H, 3, 6, D2DeDf[i] - D2Df2[i]);

        assign_block(H, 6, 0, -D2DfDe[i]);
        assign_block(H, 6, 3, D2DfDe[i] - D2Df2[i]);
        assign_block(H, 6, 6, D2Df2[i]);

        // Twist diagonal terms
        H(9, 9) = D2t1(i);
        H(10, 10) = D2t2(i);

        // Coupled terms
        assign_twist_blocks(H, 0, 9, -D2ct[i].first);
        assign_twist_blocks(H, 3, 9, D2ct[i].first);
        assign_twist_blocks(H, 6, 9, Eigen::Vector3d::Zero());

        assign_twist_blocks(H, 0, 10, -D2ct[i].second);
        assign_twist_blocks(H, 3, 10, D2ct[i].second);
        assign_twist_blocks(H, 6, 10, Eigen::Vector3d::Zero());
    }
}

// Compute second derivative tensors for each spring
std::tuple<
    std::vector<Eigen::Matrix3d>, std::vector<Eigen::Matrix3d>,
    std::vector<Eigen::Matrix3d>, std::vector<Eigen::Matrix3d>,
    Eigen::VectorXd, Eigen::VectorXd,
    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>,
    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>
> compute_second_derivatives(
    const Eigen::MatrixXd& te,
    const Eigen::MatrixXd& tf,
    const Eigen::MatrixXd& tilde_t,
    const Eigen::MatrixXd& tilde_d1,
    const Eigen::MatrixXd& tilde_d2,
    const Eigen::MatrixXd& kb,
    const Eigen::MatrixXd& m1e,
    const Eigen::MatrixXd& m1f,
    const Eigen::MatrixXd& m2e,
    const Eigen::MatrixXd& m2f,
    const Eigen::VectorXd& kappa,
    const Eigen::VectorXd& chi,
    const Eigen::VectorXd& norm_e,
    const Eigen::VectorXd& norm_f,
    bool is_kappa1
) {
    int N = te.rows();
    Eigen::VectorXd norm2_e = norm_e.array().square();
    Eigen::VectorXd norm2_f = norm_f.array().square();
    
    auto Id3 = broadcast_identity(N);
    auto tt_o_tt = batch_outer(tilde_t, tilde_t);

    Eigen::MatrixXd tf_c_dt(N, 3), te_c_dt(N, 3);
    std::vector<Eigen::Matrix3d> tf_c_dt_o_tt, tt_o_tf_c_dt, te_c_dt_o_tt, tt_o_te_c_dt;
    std::vector<Eigen::Matrix3d> kb_o_de, kb_o_df;
    std::vector<Eigen::Matrix3d> D2De2(N), D2Df2(N), D2DeDf(N);
    Eigen::VectorXd D2t1(N), D2t2(N);
    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> D2ct;

    for (int i = 0; i < N; ++i) {
        Eigen::Vector3d dt_e = is_kappa1 ? tilde_d2.row(i) : tilde_d1.row(i);
        Eigen::Vector3d dt_f = dt_e;

        // Cross products
        tf_c_dt.row(i) = Eigen::Vector3d(tf.row(i)).cross(dt_e).transpose();
        te_c_dt.row(i) = Eigen::Vector3d(te.row(i)).cross(dt_f).transpose();

        // Outer products
        Eigen::Vector3d tfcdt = tf_c_dt.row(i);
        Eigen::Vector3d ttilde = tilde_t.row(i);
        tf_c_dt_o_tt.push_back(tfcdt * ttilde.transpose());
        tt_o_tf_c_dt.push_back(ttilde * tfcdt.transpose());
        Eigen::Vector3d tecdt = te_c_dt.row(i);
        te_c_dt_o_tt.push_back(tecdt * ttilde.transpose());
        tt_o_te_c_dt.push_back(ttilde * tecdt.transpose());

        Eigen::Vector3d de = is_kappa1 ? m2e.row(i) : m1e.row(i);
        Eigen::Vector3d df = is_kappa1 ? m2f.row(i) : m1f.row(i);

        // Outer products for kb
        Eigen::Vector3d kb_row = kb.row(i);
        kb_o_de.push_back(kb_row * de.transpose());
        kb_o_df.push_back(kb_row * df.transpose());

        double k = kappa(i), c = chi(i), ne = norm_e(i), nf = norm_f(i);

        D2De2[i] = (2 * k * tt_o_tt[i] - tf_c_dt_o_tt[i] - tt_o_tf_c_dt[i]) / norm2_e(i)
                   - (k / (c * norm2_e(i))) * (Id3[i] - te.row(i).transpose() * te.row(i))
                   + (0.5 / norm2_e(i)) * kb_o_de[i];

        D2Df2[i] = (2 * k * tt_o_tt[i] + te_c_dt_o_tt[i] + tt_o_te_c_dt[i]) / norm2_f(i)
                   - (k / (c * norm2_f(i))) * (Id3[i] - tf.row(i).transpose() * tf.row(i))
                   + (0.5 / norm2_f(i)) * kb_o_df[i];

        Eigen::Matrix3d D2DeDf_i;
        if (is_kappa1) {
            Eigen::Matrix3d a = -tf_c_dt_o_tt[i];
            Eigen::Matrix3d b = tt_o_te_c_dt[i];
            Eigen::Matrix3d c = cross_mat_batch(tilde_d2)[i];
            D2DeDf_i = a + b - c;
        } else {
            Eigen::Matrix3d a = tf_c_dt_o_tt[i];
            Eigen::Matrix3d b = -tt_o_te_c_dt[i];
            Eigen::Matrix3d c = cross_mat_batch(tilde_d1)[i];
            D2DeDf_i = a + b + c;
        }
        D2DeDf[i] = (-k / (c * norm_e(i) * norm_f(i))) * (Id3[i] + Eigen::Vector3d(te.row(i)) * Eigen::Vector3d(tf.row(i)).transpose())
                    + (1.0 / (norm_e(i) * norm_f(i))) * (2 * k * tt_o_tt[i] + D2DeDf_i);

        D2t1(i) = is_kappa1 ? -0.5 * kb.row(i).dot(de) : 0.5 * kb.row(i).dot(de);
        D2t2(i) = is_kappa1 ? -0.5 * kb.row(i).dot(df) : 0.5 * kb.row(i).dot(df);

        Eigen::Vector3d t = tilde_t.row(i);
        Eigen::Vector3d d0 = is_kappa1 ? m1e.row(i) : m2e.row(i);
        Eigen::Vector3d d1 = is_kappa1 ? m1f.row(i) : m2f.row(i);

        // Fix dot and cross products
        double kb_dot_d0 = kb.row(i).dot(d0);
        double kb_dot_d1 = kb.row(i).dot(d1);
        double kb_dot_df = kb.row(i).dot(df);
        double kb_dot_de = kb.row(i).dot(de);
        Eigen::Vector3d tf_vec = Eigen::Vector3d(tf.row(i));
        Eigen::Vector3d te_vec = Eigen::Vector3d(te.row(i));
        // c0, c1, c2, c3 with correct dot and cross
        Eigen::Vector3d c0 = (0.5 * kb_dot_d0 * t) - (tf_vec.cross(d0)) / c;
        Eigen::Vector3d c1 = (0.5 * kb_dot_d1 * t) - (tf_vec.cross(d1)) / c;
        Eigen::Vector3d c2 = (0.5 * kb_dot_d0 * t) + (te_vec.cross(d0)) / c;
        Eigen::Vector3d c3 = (0.5 * kb_dot_d1 * t) + (te_vec.cross(d1)) / c;
        if (is_kappa1) {
            D2ct.emplace_back(c0, c1);
        } else {
            D2ct.emplace_back(c2, c3);
        }
    }

    return {D2De2, D2DeDf, D2DeDf, D2Df2, D2t1, D2t2, D2ct, D2ct};
}



std::pair<Eigen::MatrixXd, std::vector<Eigen::MatrixXd>>
BendEnergy::grad_hess_strain(const RobotState& state) const {
    const int N = static_cast<int>(_edges_ind.size());

    // Node positions
    auto node_pos = this->get_node_positions(state.q0);
    int M = node_pos.size();
    Eigen::MatrixXd n0p(N, 3), n1p(N, 3), n2p(N, 3);
    for (int i = 0; i < N; ++i) {
        int n0 = _node_dof_ind[3*i + 0];
        int n1 = _node_dof_ind[3*i + 1];
        int n2 = _node_dof_ind[3*i + 2];
        n0p.row(i) = node_pos[0][n0];
        n1p.row(i) = node_pos[0][n1];
        n2p.row(i) = node_pos[0][n2];
    }

    // Directors
    auto [m1e, m2e, m1f, m2f] = _get_adjusted_material_directors(state.m1, state.m2);

    // Geometry
    Eigen::MatrixXd ee = n1p - n0p, ef = n2p - n1p;
    Eigen::VectorXd norm_e = ee.rowwise().norm(), norm_f = ef.rowwise().norm();
    Eigen::MatrixXd te = ee.array().colwise() / norm_e.array();
    Eigen::MatrixXd tf = ef.array().colwise() / norm_f.array();
    Eigen::VectorXd chi = (te.array() * tf.array()).rowwise().sum() + 1.0;
    Eigen::VectorXd chi_inv = chi.cwiseInverse();

    // kb
    Eigen::MatrixXd kb(N, 3);
    for (int i = 0; i < N; ++i)
        kb.row(i) = Eigen::Vector3d(te.row(i)).cross(Eigen::Vector3d(tf.row(i))) * 2.0 * chi_inv(i);

    Eigen::MatrixXd tilde_t = (te + tf).array().colwise() * chi_inv.array();
    Eigen::MatrixXd tilde_d1 = (m1e + m1f).array().colwise() * chi_inv.array();
    Eigen::MatrixXd tilde_d2 = (m2e + m2f).array().colwise() * chi_inv.array();

    Eigen::VectorXd kappa1 = 0.5 * ((kb.array() * (m2e + m2f).array()).rowwise().sum());
    Eigen::VectorXd kappa2 = -0.5 * ((kb.array() * (m1e + m1f).array()).rowwise().sum());

    // First derivative of kappa1
    Eigen::MatrixXd Dkappa1De = (-kappa1).asDiagonal() * tilde_t;
    Eigen::MatrixXd Dkappa1Df = (-kappa1).asDiagonal() * tilde_t;
    Eigen::MatrixXd Dkappa2De = (-kappa2).asDiagonal() * tilde_t;
    Eigen::MatrixXd Dkappa2Df = (-kappa2).asDiagonal() * tilde_t;
    // Temporary storage for cross products
    Eigen::MatrixXd cross_te_tilde_d2(N, 3), cross_tf_tilde_d2(N, 3);
    Eigen::MatrixXd cross_te_tilde_d1(N, 3), cross_tf_tilde_d1(N, 3);
    for (int i = 0; i < N; ++i) {
        cross_te_tilde_d2.row(i) = Eigen::Vector3d(te.row(i)).cross(Eigen::Vector3d(tilde_d2.row(i))).transpose();
        cross_tf_tilde_d2.row(i) = Eigen::Vector3d(tf.row(i)).cross(Eigen::Vector3d(tilde_d2.row(i))).transpose();
        cross_te_tilde_d1.row(i) = Eigen::Vector3d(te.row(i)).cross(Eigen::Vector3d(tilde_d1.row(i))).transpose();
        cross_tf_tilde_d1.row(i) = Eigen::Vector3d(tf.row(i)).cross(Eigen::Vector3d(tilde_d1.row(i))).transpose();
    }
    Dkappa1De += cross_te_tilde_d2;
    Dkappa1De = Dkappa1De.array().colwise() / norm_e.array();
    Dkappa1Df -= cross_tf_tilde_d2;
    Dkappa1Df = Dkappa1Df.array().colwise() / norm_f.array();

    Dkappa2De -= cross_te_tilde_d1;
    Dkappa2De = Dkappa2De.array().colwise() / norm_e.array();
    Dkappa2Df += cross_tf_tilde_d1;
    Dkappa2Df = Dkappa2Df.array().colwise() / norm_f.array();

    Eigen::MatrixXd gradKappa = Eigen::MatrixXd::Zero(N, 22);
    for (int i = 0; i < N; ++i) {
        gradKappa.block(i, 0, 1, 3) = -Dkappa1De.row(i);
        gradKappa.block(i, 3, 1, 3) = Dkappa1De.row(i) - Dkappa1Df.row(i);
        gradKappa.block(i, 6, 1, 3) = Dkappa1Df.row(i);
        gradKappa(i, 9) = -0.5 * kb.row(i).dot(state.m1[_edges_ind(i, 0)]);
        gradKappa(i, 10) = -0.5 * kb.row(i).dot(state.m1[_edges_ind(i, 1)]);

        gradKappa.block(i, 11 + 0, 1, 3) = -Dkappa2De.row(i);
        gradKappa.block(i, 11 + 3, 1, 3) = Dkappa2De.row(i) - Dkappa2Df.row(i);
        gradKappa.block(i, 11 + 6, 1, 3) = Dkappa2Df.row(i);
        gradKappa(i, 11 + 9) = -0.5 * kb.row(i).dot(state.m2[_edges_ind(i, 0)]);
        gradKappa(i, 11 + 10) = -0.5 * kb.row(i).dot(state.m2[_edges_ind(i, 1)]);
    }

    std::vector<Eigen::MatrixXd> hessKappa(N, Eigen::MatrixXd::Zero(11, 22));

    // Compute second derivatives for kappa1 and kappa2
    auto [D2De2_1, D2DeDf_1, D2DfDe_1, D2Df2_1, D2t1_1, D2t2_1, D2ct1, _] =
        compute_second_derivatives(te, tf, tilde_t, tilde_d1, tilde_d2,
                                kb, m1e, m1f, m2e, m2f, kappa1, chi, norm_e, norm_f,
                                /*is_kappa1=*/true);

    auto [D2De2_2, D2DeDf_2, D2DfDe_2, D2Df2_2, D2t1_2, D2t2_2, D2ct2, __] =
        compute_second_derivatives(te, tf, tilde_t, tilde_d1, tilde_d2,
                                kb, m1e, m1f, m2e, m2f, kappa2, chi, norm_e, norm_f,
                                /*is_kappa1=*/false);

    // Fill in the Hessian matrices
    batch_assign_blocks(hessKappa, D2De2_1, D2DeDf_1, D2DfDe_1, D2Df2_1, D2t1_1, D2t2_1, D2ct1);
    for (int i = 0; i < N; ++i) {
        Eigen::MatrixXd H2 = Eigen::MatrixXd::Zero(11, 11);
        assign_second_derivatives(
            D2De2_2[i], D2DeDf_2[i], D2DfDe_2[i], D2Df2_2[i],
            D2t1_2(i), D2t2_2(i),
            D2ct2[i].first, Eigen::Vector3d::Zero(),
            D2ct2[i].second, Eigen::Vector3d::Zero(),
            H2);
        hessKappa[i].block(0, 11, 11, 11) = H2;
    }


    return {gradKappa, hessKappa};
}

