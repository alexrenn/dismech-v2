#include "robotState.h"

RobotState::RobotState(const Eigen::VectorXd& q0,
                       const Eigen::VectorXd& u_in,
                       const std::vector<Eigen::Vector3d>& a_in,
                       const std::vector<Eigen::Vector3d>& a1_in,
                       const std::vector<Eigen::Vector3d>& a2_in,
                       const std::vector<Eigen::Vector3d>& m1_in,
                       const std::vector<Eigen::Vector3d>& m2_in,
                       const Eigen::VectorXi& ref_twist_in,
                       const Eigen::VectorXd& tau_in,
                       const Eigen::VectorXi& free_dof_in)
{}
