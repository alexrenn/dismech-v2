#include "stretch_energy.h"
#include <iostream>

StretchEnergy::StretchEnergy(const std::vector<StretchSpring>& springs, const RobotState& initial_state)
    : ElasticEnergy(
        // Construct K matrix
        [&]() {
            const int N = static_cast<int>(springs.size());
            Eigen::VectorXd K(N);
            for (int i = 0; i < N; ++i) {
                K[i] = springs[i].ref_len * springs[i].EA;
            }
            return K;
        }(),
        // nodes_ind
        [&]() {
            std::vector<std::vector<int>> out;
            for (const auto& s : springs) {
                out.push_back({s.nodes_ind[0], s.nodes_ind[1]});
            }
            return out;
        }(),
        // ind
        [&]() {
            std::vector<std::vector<int>> out;
            for (const auto& s : springs) out.push_back(s.ind);
            return out;
        }(),
        initial_state),
      springs_ref_(springs)
{
    inv_l_k_.resize(springs.size());
    for (size_t i = 0; i < springs.size(); ++i) {
        inv_l_k_[i] = 1.0 / springs[i].ref_len;
    }
}

Eigen::MatrixXd StretchEnergy::get_strain(const RobotState& state) const {
    Eigen::MatrixXd strain(springs_ref_.size(), 1);
    auto node_pos = get_node_positions(state.q0);
    for (size_t i = 0; i < springs_ref_.size(); ++i) {
        Eigen::Vector3d edge = node_pos[1][i] - node_pos[0][i];
        double edge_len = edge.norm();
        strain(i, 0) = edge_len * inv_l_k_[i] - 1.0;
    }
    return strain;
}
std::pair<Eigen::MatrixXd, std::vector<Eigen::MatrixXd>> StretchEnergy::grad_hess_strain(const RobotState& state) const {
    Eigen::MatrixXd grad_eps(springs_ref_.size(), 6);
    std::vector<Eigen::MatrixXd> hess_eps(springs_ref_.size(), Eigen::MatrixXd::Zero(6, 6));

    auto node_pos = get_node_positions(state.q0);  // only call once

    for (size_t i = 0; i < springs_ref_.size(); ++i) {
        const Eigen::Vector3d& x0 = node_pos[0][i];  // first node in spring i
        const Eigen::Vector3d& x1 = node_pos[1][i];  // second node in spring i

        Eigen::Vector3d edge = x1 - x0;
        double edge_len = edge.norm();
        Eigen::Vector3d tangent = edge / edge_len;
        double eps = edge_len * inv_l_k_[i] - 1.0;

        Eigen::Vector3d eps_unit = tangent * inv_l_k_[i];
        grad_eps.block<1, 3>(i, 0) = -eps_unit.transpose();
        grad_eps.block<1, 3>(i, 3) = eps_unit.transpose();

        Eigen::Matrix3d edge_outer = edge * edge.transpose();
        double len3 = std::pow(edge_len, 3);
        Eigen::Matrix3d M = ((inv_l_k_[i] - 1.0 / edge_len) * Eigen::Matrix3d::Identity() + edge_outer / len3) * (2.0 * inv_l_k_[i]);
        Eigen::Matrix3d M2 = M - 2.0 * (eps_unit * eps_unit.transpose());

        if (std::abs(eps) > 1e-10) {
            M2 /= eps;
        } else {
            M2.setZero();
        }
        M2 *= 0.5;

        hess_eps[i].block<3, 3>(0, 0) = M2;
        hess_eps[i].block<3, 3>(3, 3) = M2;
        hess_eps[i].block<3, 3>(0, 3) = -M2;
        hess_eps[i].block<3, 3>(3, 0) = -M2;
    }

    return {grad_eps, hess_eps};
}
