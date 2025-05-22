#ifndef ROBOT_STATE_H
#define ROBOT_STATE_H

#include "../eigenIncludes.h"
#include <vector>

struct RobotState {
    Eigen::VectorXd q0, u, tau;
    Eigen::VectorXi ref_twist, free_dof;
    std::vector<Eigen::Vector3d> a, a1, a2, m1, m2;

    // RobotState(); // default
    RobotState(const Eigen::VectorXd& q0,
               const Eigen::VectorXd& u,
               const std::vector<Eigen::Vector3d>& a,
               const std::vector<Eigen::Vector3d>& a1,
               const std::vector<Eigen::Vector3d>& a2,
               const std::vector<Eigen::Vector3d>& m1,
               const std::vector<Eigen::Vector3d>& m2,
               const Eigen::VectorXi& ref_twist,
               const Eigen::VectorXd& tau,
               const Eigen::VectorXi& free_dof);

    // static RobotState init(const Eigen::VectorXd& q0,
    //                    const Eigen::VectorXd& u,
    //                    const std::vector<Eigen::VectorXd>& a,
    //                    const std::vector<Eigen::VectorXd>& a1,
    //                    const std::vector<Eigen::VectorXd>& a2,
    //                    const std::vector<Eigen::VectorXd>& m1,
    //                    const std::vector<Eigen::VectorXd>& m2,
    //                    const Eigen::VectorXi& ref_twist,
    //                    const Eigen::VectorXd& free_dof);
};

#endif // ROBOT_STATE_H