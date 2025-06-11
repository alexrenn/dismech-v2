#include "frame_util.h"
#include <iostream>


std::vector<Eigen::Vector3d> parallel_transport(
    const std::vector<Eigen::Vector3d>& u,
    const std::vector<Eigen::Vector3d>& t_start,
    const std::vector<Eigen::Vector3d>& t_end) {

    size_t N = u.size();
    std::vector<Eigen::Vector3d> result(N);

    for (size_t i = 0; i < N; ++i) {
        const Eigen::Vector3d& ti = t_start[i];
        const Eigen::Vector3d& tj = t_end[i];
        const Eigen::Vector3d& ui = u[i];

        // b = ti x tj
        Eigen::Vector3d b = ti.cross(tj);
        double b_norm = b.norm();
        bool is_degenerate = (b_norm < 1e-10);

        if (is_degenerate) {
            result[i] = ui;  // No transport needed
            continue;
        }

        Eigen::Vector3d b_normalized = b / b_norm;

        // Orthogonalize b_normalized against ti
        double dot_bt = b_normalized.dot(ti);
        Eigen::Vector3d b_ortho = b_normalized - dot_bt * ti;

        double b_ortho_norm = b_ortho.norm();
        if (b_ortho_norm < 1e-10) b_ortho_norm = 1.0; // safe normalization
        Eigen::Vector3d b_ortho_normalized = b_ortho / b_ortho_norm;

        // Compute transport basis vectors
        Eigen::Vector3d n1 = ti.cross(b_ortho_normalized);
        Eigen::Vector3d n2 = tj.cross(b_ortho_normalized);

        // Project ui onto (ti, n1, b_ortho) and transport
        double comp_ti   = ui.dot(ti);
        double comp_n1   = ui.dot(n1);
        double comp_bort = ui.dot(b_ortho_normalized);

        result[i] = comp_ti * tj + comp_n1 * n2 + comp_bort * b_ortho_normalized;
    }

    return result;
}

Eigen::Vector3d parallel_transport(const Eigen::Vector3d& u,
                                   const Eigen::Vector3d& t_start,
                                   const Eigen::Vector3d& t_end) {
    Eigen::Vector3d b = t_start.cross(t_end);
    double b_norm = b.norm();

    if (b_norm < 1e-10) {
        return u;
    }

    Eigen::Vector3d b_normalized = b / b_norm;
    Eigen::Vector3d b_ortho = b_normalized - b_normalized.dot(t_start) * t_start;

    double b_ortho_norm = b_ortho.norm();
    if (b_ortho_norm < 1e-10) {
        return u;
    }

    Eigen::Vector3d b_ortho_normalized = b_ortho / b_ortho_norm;

    Eigen::Vector3d n1 = t_start.cross(b_ortho_normalized);
    Eigen::Vector3d n2 = t_end.cross(b_ortho_normalized);

    double comp_ti   = u.dot(t_start);
    double comp_n1   = u.dot(n1);
    double comp_bort = u.dot(b_ortho_normalized);

    return comp_ti * t_end + comp_n1 * n2 + comp_bort * b_ortho_normalized;
}



void compute_tfc_midedge(
    const std::vector<Eigen::Matrix3d>& p_s,                // [N faces] each 3x3
    const std::vector<Eigen::Matrix3d>& tau0_s,             // [N faces] each 3x3
    const std::vector<std::array<int, 3>>& s_s,          // [N faces] 3 signs
    std::vector<Eigen::Matrix3d>& ts,                       // output: tangents [N, 3x3]
    std::vector<std::array<double, 3>>& fs,                 // output: force projections
    std::vector<std::array<double, 3>>& cs                  // output: scalar coefficients
) {
    size_t N = p_s.size();
    ts.resize(N);
    fs.resize(N);
    cs.resize(N);

    for (size_t i = 0; i < N; ++i) {
        const auto& p = p_s[i];
        const auto& tau0 = tau0_s[i];
        const auto& s = s_s[i];

        // Sign-adjusted tau0
        Eigen::Matrix3d tau;
        for (int j = 0; j < 3; ++j) {
            tau.col(j) = static_cast<double>(s[j]) * tau0.col(j);  // Convert int to double for multiplication
        }

        Eigen::Vector3d vi = p.col(2) - p.col(1);
        Eigen::Vector3d vj = p.col(0) - p.col(2);
        Eigen::Vector3d vk = p.col(1) - p.col(0);

        double li = vi.norm();
        double lj = vj.norm();
        double lk = vk.norm();

        Eigen::Vector3d normal = vk.cross(vi);
        double norm_normal = normal.norm();
        double A = 0.5 * norm_normal;
        Eigen::Vector3d unit_norm = normal / norm_normal;

        Eigen::Vector3d t_i = vi.cross(unit_norm);
        Eigen::Vector3d t_j = vj.cross(unit_norm);
        Eigen::Vector3d t_k = vk.cross(unit_norm);

        Eigen::Vector3d t_i_n = t_i.normalized();
        Eigen::Vector3d t_j_n = t_j.normalized();
        Eigen::Vector3d t_k_n = t_k.normalized();

        double dot_i = t_i_n.dot(static_cast<double>(s[0]) * tau0.col(0));
        double dot_j = t_j_n.dot(static_cast<double>(s[1]) * tau0.col(1));
        double dot_k = t_k_n.dot(static_cast<double>(s[2]) * tau0.col(2));

        double c_i = 1.0 / (A * li * dot_i);
        double c_j = 1.0 / (A * lj * dot_j);
        double c_k = 1.0 / (A * lk * dot_k);

        double f_i = unit_norm.dot(static_cast<double>(s[0]) * tau0.col(0));
        double f_j = unit_norm.dot(static_cast<double>(s[1]) * tau0.col(1));
        double f_k = unit_norm.dot(static_cast<double>(s[2]) * tau0.col(2));

        // Store results
        Eigen::Matrix3d t_out;
        t_out.col(0) = t_i;
        t_out.col(1) = t_j;
        t_out.col(2) = t_k;
        ts[i] = t_out;

        fs[i] = {f_i, f_j, f_k};
        cs[i] = {c_i, c_j, c_k};
    }
}

    Eigen::VectorXi compute_reference_twist_util(
    const std::vector<std::array<int, 2>>& edges,
    const std::vector<std::array<int, 2>>& sgn,
    const std::vector<Eigen::Vector3d>& a1,
    const std::vector<Eigen::Vector3d>& tangent,
    const Eigen::VectorXi& ref_twist)
{
    size_t n = ref_twist.size();
    Eigen::VectorXi updated_twist(n);

    for (size_t i = 0; i < n; ++i) {
        int e0 = edges[i][0];
        int e1 = edges[i][1];

        Eigen::Vector3d t0 = tangent[e0] * static_cast<double>(sgn[i][0]);
        Eigen::Vector3d t1 = tangent[e1] * static_cast<double>(sgn[i][1]);

        Eigen::Vector3d u0 = a1[e0];
        Eigen::Vector3d u1 = a1[e1];

        // Parallel transport u0 along t0 → t1
        Eigen::Vector3d ut = parallel_transport(u0, t0, t1);

        // Rotate transported u0 by ref_twist angle around t1
        double twist_angle = static_cast<double>(ref_twist(i));
        ut = rotate_axis_angle(ut, t1, twist_angle);

        // Compute signed angle between rotated u0 and u1
        double angle = signed_angle(ut, u1, t1);

        // Add angle to ref_twist, and round to integer (if needed)
        updated_twist(i) = static_cast<int>(std::round(twist_angle + angle));
    }

    return updated_twist;
    }

double signed_angle(const Eigen::Vector3d& u,
                const Eigen::Vector3d& v,
                const Eigen::Vector3d& n)
{
    Eigen::Vector3d w = u.cross(v);
    double norm_w = w.norm();

    double dot_uv = u.dot(v);
    double safe_denominator = (std::abs(dot_uv) < 1e-10) ? 1.0 : dot_uv;

    double angle = std::atan2(norm_w, safe_denominator);
    double sign = std::copysign(1.0, n.dot(w));  // gives correct sign

    return angle * sign;
}

    Eigen::Vector3d rotate_axis_angle(const Eigen::Vector3d& v,
                                  const Eigen::Vector3d& axis,
                                  double theta)
    {
        double cos_theta = std::cos(theta);
        double sin_theta = std::sin(theta);
        double dot_product = axis.dot(v);

        return cos_theta * v +
            sin_theta * axis.cross(v) +
            (1.0 - cos_theta) * dot_product * axis;
    }

    