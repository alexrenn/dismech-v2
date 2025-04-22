#include <Eigen/Dense>
#include <array>
#include <vector>

struct RobotState {
    Eigen::VectorXd q;
    Eigen::VectorXd u;
    Eigen::VectorXd a;
    Eigen::VectorXd a1;
    Eigen::VectorXd a2;
    Eigen::VectorXd m1;
    Eigen::VectorXd m2;
    Eigen::VectorXd ref_twist;
    Eigen::VectorXd tau;
    std::vector<int> free_dof;

    // Constructor to initialize RobotState
    RobotState(const Eigen::VectorXd& q0,
               const Eigen::VectorXd& a1,
               const Eigen::VectorXd& a2,
               const Eigen::VectorXd& m1,
               const Eigen::VectorXd& m2,
               const Eigen::VectorXd& ref_twist,
               const Eigen::VectorXd& tau)
        : q(q0),
          u(Eigen::VectorXd::Zero(q0.size())),
          a(Eigen::VectorXd::Zero(q0.size())),
          a1(a1),
          a2(a2),
          m1(m1),
          m2(m2),
          ref_twist(ref_twist),
          tau(tau),
          free_dof(std::vector<int>(q0.size())) {
        // Fill the free_dof with indices from 0 to q0.size()-1
        for (int i = 0; i < q0.size(); ++i) {
            free_dof[i] = i;
        }
    }

    // Static initialization method
    static RobotState init(const Eigen::VectorXd& q0,
                           const Eigen::VectorXd& a1,
                           const Eigen::VectorXd& a2,
                           const Eigen::VectorXd& m1,
                           const Eigen::VectorXd& m2,
                           const Eigen::VectorXd& ref_twist,
                           const Eigen::VectorXd& tau) {
        return RobotState(q0, a1, a2, m1, m2, ref_twist, tau);
    }
};
