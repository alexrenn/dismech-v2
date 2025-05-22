#include "frame_util.h"


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
