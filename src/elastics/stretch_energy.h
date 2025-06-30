#ifndef STRETCH_ENERGY_H
#define STRETCH_ENERGY_H

#include "../eigenIncludes.h"
#include "elastic_energy.h"
#include "../springs/stretch_spring.h"
#include "../mechanics/robotState.h"
#include <vector>
#include <array>

class StretchEnergy : public ElasticEnergy {
public:
    StretchEnergy(const std::vector<StretchSpring>& springs, const RobotState& initial_state);

    Eigen::MatrixXd get_strain(const RobotState& state) const override;

    std::pair<Eigen::MatrixXd, std::vector<Eigen::MatrixXd>> grad_hess_strain(const RobotState& state) const override;

    Eigen::VectorXd get_K() const;

private:
    Eigen::VectorXd inv_l_k_; 
    const std::vector<StretchSpring>& springs_ref_;

    // std::vector<std::vector<Eigen::Vector3d>>  get_node_pos(const Eigen::VectorXd& q) const;
};

#endif // STRETCH_ENERGY_H
