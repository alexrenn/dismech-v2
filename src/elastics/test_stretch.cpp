#include "stretch_energy.h"
#include <iostream>
#include "../eigenIncludes.h"
#include "../mechanics/robotState.h"
#include "../springs/stretch_spring.h"

int main() {
    // === Define springs ===
    std::vector<std::array<int, 2>> node_pairs = { {0, 1}, {2, 3} };
    std::vector<double> ref_lens = {1.0, 1.0};
    std::vector<double> EAs = {10.0, 20.0};

    auto map_node_to_dof = [](int node) {
        return std::array<int, 3>{3 * node, 3 * node + 1, 3 * node + 2};
    };

    std::vector<StretchSpring> springs;
    for (size_t i = 0; i < node_pairs.size(); ++i) {
        springs.emplace_back(node_pairs[i], ref_lens[i], EAs[i], map_node_to_dof);
    }

    // === Initial robot state ===
    Eigen::VectorXd q0(12);
    q0 << 0,0,0, 1,0,0, 0,1,0, 1,1,0; // 4 nodes in XY plane
    Eigen::VectorXd u = q0;
    std::vector<Eigen::Vector3d> a(4, Eigen::Vector3d::Zero());
    std::vector<Eigen::Vector3d> a1(4, Eigen::Vector3d::Zero());
    std::vector<Eigen::Vector3d> a2(4, Eigen::Vector3d::Zero());
    std::vector<Eigen::Vector3d> m1(4, Eigen::Vector3d::UnitY());
    std::vector<Eigen::Vector3d> m2(4, Eigen::Vector3d::UnitZ());
    Eigen::VectorXi ref_twist = Eigen::VectorXi::Zero(4);
    std::vector<Eigen::Vector3d> tau(4, Eigen::Vector3d::Zero());
    Eigen::VectorXi free_dof = Eigen::VectorXi::LinSpaced(12, 0, 11);

    RobotState state(q0, u, a, a1, a2, m1, m2, ref_twist, tau, free_dof);

    // === Create StretchEnergy object ===
    StretchEnergy stretch_energy(springs, state);  // using first constructor (vector version)

    // === Run get_strain ===
    Eigen::MatrixXd strain = stretch_energy.get_strain(state);
    std::cout << "Strain:\n" << strain << "\n\n";

    // === Run grad_hess_strain ===
    auto [grad, hess] = stretch_energy.grad_hess_strain(state);
    for (int i = 0; i < hess.size(); ++i) {
        std::cout << "Gradient of strain (row " << i << "):\n" << grad.row(i) << "\n\n";
        std::cout << "Hessian block (" << i << "):\n" << hess[i] << "\n";
    }

    std::cout << "Gradient shape: " << grad.rows() << " x " << grad.cols() << std::endl;
    std::cout << "Hessian shape: " << hess.size() << " x " << hess[0].rows() << " x " << hess[0].cols() << std::endl;

    return 0;
}
