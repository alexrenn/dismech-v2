#include <Eigen/Dense>
#include <vector>

struct RobotState {
    Eigen::VectorXd q, u, a, a1, a2, m1, m2, ref_twist, tau;
    std::vector<int> free_dof;

    RobotState(); // default
    RobotState(const Eigen::VectorXd& q0,
               const Eigen::VectorXd& a1,
               const Eigen::VectorXd& a2,
               const Eigen::VectorXd& m1,
               const Eigen::VectorXd& m2,
               const Eigen::VectorXd& ref_twist,
               const Eigen::VectorXd& tau);

    static RobotState init(const Eigen::VectorXd& q0,
                           const Eigen::VectorXd& a1,
                           const Eigen::VectorXd& a2,
                           const Eigen::VectorXd& m1,
                           const Eigen::VectorXd& m2,
                           const Eigen::VectorXd& ref_twist,
                           const Eigen::VectorXd& tau);
};
