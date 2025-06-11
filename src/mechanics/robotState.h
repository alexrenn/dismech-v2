#ifndef ROBOT_STATE_H
#define ROBOT_STATE_H

#include "../eigenIncludes.h"
#include <vector>

struct RobotState {
    Eigen::VectorXd q0, u;
    std::vector<Eigen::Vector3d> a, a1, a2, m1, m2, tau;
    Eigen::VectorXi ref_twist, free_dof;

    // RobotState(); // default
    RobotState(const Eigen::VectorXd& q0,
               const Eigen::VectorXd& u,
               const std::vector<Eigen::Vector3d>& a,
               const std::vector<Eigen::Vector3d>& a1,
               const std::vector<Eigen::Vector3d>& a2,
               const std::vector<Eigen::Vector3d>& m1,
               const std::vector<Eigen::Vector3d>& m2,
               const Eigen::VectorXi& ref_twist,
               const std::vector<Eigen::Vector3d>& tau,
               const Eigen::VectorXi& free_dof);
};

#endif // ROBOT_STATE_H