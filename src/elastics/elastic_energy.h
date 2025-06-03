#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <memory>
#include <tuple>
#include "../mechanics/RobotState.h"

class ElasticEnergy {
public:
    ElasticEnergy(const Eigen::MatrixXd& K,
                  const Eigen::MatrixXi& nodes_ind,
                  const Eigen::MatrixXi& ind,
                  const RobotState& initial_state,
                  std::function<Eigen::MatrixXd(const Eigen::MatrixXd&)> get_strain = nullptr);

    virtual ~ElasticEnergy() = default;

    virtual Eigen::MatrixXd get_strain(const RobotState& state) const = 0;
    virtual std::pair<Eigen::MatrixXd, Eigen::Tensor<double, 3>> grad_hess_strain(const RobotState& state) const = 0;

    double get_energy_linear_elastic(const RobotState& state, bool output_scalar = true) const;

    std::pair<Eigen::VectorXd, Eigen::SparseMatrix<double>> grad_hess_energy_linear_elastic_sparse(const RobotState& state) const;

    void fdm_check_grad_hess_strain(const RobotState& state,
                                    Eigen::MatrixXd& grad_FDM,
                                    Eigen::Tensor<double, 3>& hess_FDM,
                                    bool& grad_match, bool& hess_match) const;

protected: //not accessible to the user or derived classes
    Eigen::MatrixXd K_;
    int n_K_;
    Eigen::VectorXi node_dof_ind_;
    int n_nodes_;
    RobotState initial_state_;
    Eigen::MatrixXi ind_;

    std::vector<int> rows_, cols_;
    Eigen::MatrixXd nat_strain_;

    Eigen::MatrixXd get_node_pos(const Eigen::VectorXd& q) const;
};
