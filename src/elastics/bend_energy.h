#ifndef BEND_ENERGY_H
#define BEND_ENERGY_H

#include "elastic_energy.h"
#include "../springs/bend_twist_spring.h" // You need to define this struct/class
#include "../mechanics/robotState.h"
#include <vector>
#include <array>
#include <Eigen/Dense>

class BendEnergy : public ElasticEnergy {
public:
    BendEnergy(const std::vector<BendTwistSpring>& springs, const RobotState& initial_state, 
               std::function<Eigen::MatrixXd(const Eigen::MatrixXd&)> get_strain = nullptr);

    // Override virtuals from ElasticEnergy
    Eigen::MatrixXd get_strain(const RobotState& state) const override;
    std::pair<Eigen::MatrixXd, std::vector<Eigen::MatrixXd>> grad_hess_strain(const RobotState& state) const override;

protected:
    Eigen::MatrixXi _sgn;         // N x 2
    Eigen::MatrixXi _edges_ind;   // N x 2
    Eigen::MatrixXd _sign_grad;   // N x 11
    int N;

    std::vector<Eigen::MatrixXd> _sign_hess;

    // Helper
    std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> _get_adjusted_material_directors(
            const std::vector<Eigen::Vector3d>& m1, const std::vector<Eigen::Vector3d>& m2) const;
};

#endif // BEND_ENERGY_H