#ifndef ELASTIC_ENERGY_H
#define ELASTIC_ENERGY_H

#include <vector>
#include <array>
#include <memory>

#include "../mechanics/robotState.h"
#include "../mechanics/softRobot.h" 
#include "../mechanics/environment.h"
#include "../mechanics/geometry.h"
#include "../eigenIncludes.h"

class ElasticEnergy {
public:
    ElasticEnergy(const Eigen::MatrixXd& K,
                  const std::vector<std::vector<int>>& nodes_ind,
                  const std::vector<std::vector<int>>& ind,
                  const RobotState& initial_state,
                  std::function<Eigen::MatrixXd(const Eigen::MatrixXd&)> get_strain = nullptr);

    virtual ~ElasticEnergy() = default;

    // Abstract strain-related methods
    virtual Eigen::MatrixXd get_strain(const RobotState& state) const = 0;
    virtual std::pair<Eigen::MatrixXd, std::vector<Eigen::MatrixXd>> grad_hess_strain(const RobotState& state) const = 0;

    void set_nat_strain(const Eigen::MatrixXd& strain);

    // one function in python
    double get_energy_linear_elastic_scalar(const RobotState& state) const;
    Eigen::MatrixXd get_energy_linear_elastic_matrix(const RobotState& state) const;

    std::pair<Eigen::VectorXd, Eigen::MatrixXd> grad_hess_energy_linear_elastic_dense(const RobotState& state) const;
    std::pair<Eigen::VectorXd, Eigen::SparseMatrix<double>> grad_hess_energy_linear_elastic_sparse(const RobotState& state) const;


protected:
    Eigen::MatrixXd _K;
    int _n_K;
    int _n_nodes;
    int _n_elems;

    std::vector<std::vector<int>> _nodes_ind;
    std::vector<int> _node_dof_ind;
    std::vector<std::vector<int>> _ind;

    std::vector<int> _rows;
    std::vector<int> _cols;

    RobotState _initial_state;
    Eigen::MatrixXd _nat_strain;

    void compute_post_init_strain(std::function<Eigen::MatrixXd(const Eigen::MatrixXd&)> get_strain);

    std::vector<std::vector<Eigen::Vector3d>> get_node_positions(const Eigen::VectorXd& q) const;
};

#endif