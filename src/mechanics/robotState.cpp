#include "robotState.h"

RobotState::RobotState(const Eigen::VectorXd& q0_in,
                       const Eigen::VectorXd& u_in,
                       const std::vector<Eigen::Vector3d>& a_in,
                       const std::vector<Eigen::Vector3d>& a1_in,
                       const std::vector<Eigen::Vector3d>& a2_in,
                       const std::vector<Eigen::Vector3d>& m1_in,
                       const std::vector<Eigen::Vector3d>& m2_in,
                       const Eigen::VectorXi& ref_twist_in,
                       const std::vector<Eigen::Vector3d>& tau_in,
                       const Eigen::VectorXi& free_dof_in)
    : q0(q0_in), u(u_in), a(a_in), a1(a1_in), a2(a2_in), 
      m1(m1_in), m2(m2_in), ref_twist(ref_twist_in), 
      tau(tau_in), free_dof(free_dof_in) {}
