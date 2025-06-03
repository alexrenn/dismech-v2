#include "elastic_energy.h"
#include "softRobot.h"
#include "../eigenIncludes.h"
#include <iostream>

ElasticEnergy::ElasticEnergy(const Eigen::MatrixXd& K,
                             const Eigen::MatrixXi& nodes_ind,
                             const Eigen::MatrixXi& ind,
                             const RobotState& initial_state,
                             std::function<Eigen::MatrixXd(const Eigen::MatrixXd&)> get_strain)
    : K_(K),
      n_K_(K.cols() > 1 ? K.cols() : 1),
      initial_state_(initial_state),
      ind_(ind) {
    
    node_dof_ind_ = SoftRobot::map_node_to_dof(nodes_ind);
    n_nodes_ = nodes_ind.cols();

    int stencil_n_dof = ind_.cols();
    for (int i = 0; i < ind_.rows(); ++i) {
        for (int j = 0; j < stencil_n_dof; ++j) {
            for (int k = 0; k < stencil_n_dof; ++k) {
                rows_.push_back(ind_(i, j));
                cols_.push_back(ind_(i, k));
            }
        }
    }

    if (get_strain) {
        Eigen::MatrixXd initial_pos = get_node_pos(initial_state.q);
        nat_strain_ = get_strain(initial_pos);
    }
}

Eigen::MatrixXd ElasticEnergy::get_node_pos(const Eigen::VectorXd& q) const {
    Eigen::MatrixXd positions(n_nodes_, 3);
    for (int i = 0; i < n_nodes_; ++i) {
        for (int d = 0; d < 3; ++d) {
            positions(i, d) = q[node_dof_ind_(3 * i + d)];
        }
    }
    return positions;
}

double ElasticEnergy::get_energy_linear_elastic(const RobotState& state, bool output_scalar) const {
    Eigen::MatrixXd strain = get_strain(state);
    Eigen::MatrixXd del_strain = (strain - nat_strain_).reshaped(strain.rows(), n_K_);

    if (output_scalar) {
        return 0.5 * (K_.reshaped(strain.rows(), n_K_).array() * del_strain.array().square()).sum();
    } else {
        // Per-element energy would be returned here
        return 0.0; // Adjust as needed
    }
}

std::pair<Eigen::VectorXd, Eigen::SparseMatrix<double>> ElasticEnergy::grad_hess_energy_linear_elastic_sparse(const RobotState& state) const {
    Eigen::MatrixXd strain = get_strain(state);
    auto [grad_strain, hess_strain] = grad_hess_strain(state);

    Eigen::MatrixXd del_strain = (strain - nat_strain_).reshaped(strain.rows(), n_K_);
    Eigen::MatrixXd gradE_strain = K_.reshaped(del_strain.rows(), n_K_).array() * del_strain.array();

    // TODO: Complete gradient and Hessian computation based on analytical expressions
    Eigen::VectorXd grad(state.q.size());
    grad.setZero();

    for (int i = 0; i < ind_.rows(); ++i) {
        for (int j = 0; j < ind_.cols(); ++j) {
            grad[ind_(i, j)] -= gradE_strain(i, j % n_K_);
        }
    }

    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < rows_.size(); ++i) {
        triplets.emplace_back(rows_[i], cols_[i], 0);  // Placeholder: Fill with computed -hess_energy
    }

    Eigen::SparseMatrix<double> hess(state.q.size(), state.q.size());
    hess.setFromTriplets(triplets.begin(), triplets.end());
    return {grad, hess};
}

void ElasticEnergy::fdm_check_grad_hess_strain(const RobotState& state,
                                               Eigen::MatrixXd& grad_FDM,
                                               Eigen::Tensor<double, 3>& hess_FDM,
                                               bool& grad_match,
                                               bool& hess_match) const {
    constexpr double change = 1e-8;
    Eigen::MatrixXd strain = get_strain(state);
    auto [grad_strain, hess_strain] = grad_hess_strain(state);

    grad_FDM = Eigen::MatrixXd::Zero(grad_strain.rows(), grad_strain.cols());
    hess_FDM = Eigen::Tensor<double, 3>(hess_strain.dimension(0), hess_strain.dimension(1), hess_strain.dimension(2));
    hess_FDM.setZero();

    for (int i = 0; i < grad_strain.cols(); ++i) {
        Eigen::VectorXd q_perturbed = state.q;
        for (int j = 0; j < ind_.rows(); ++j) {
            q_perturbed[ind_(j, i)] += change;
        }

        RobotState perturbed_state = state;
        perturbed_state.q = q_perturbed;

        Eigen::MatrixXd perturbed_strain = get_strain(perturbed_state);
        auto [perturbed_grad_strain, _] = grad_hess_strain(perturbed_state);

        grad_FDM.col(i) = (perturbed_strain - strain) / change;
        hess_FDM.chip(i, 2) = (perturbed_grad_strain - grad_strain) / change;
    }

    grad_match = (grad_FDM - grad_strain).cwiseAbs().mean() < 1e-3 * grad_strain.cwiseAbs().mean();
    hess_match = (hess_FDM - hess_strain).abs().mean() < 1e-3 * hess_strain.abs().mean();
}
