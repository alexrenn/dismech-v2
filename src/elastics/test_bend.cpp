/* File used for testing bend energy class*/

#include "bend_energy.h"
#include <iostream>
#include "../eigenIncludes.h"
#include "../springs/bend_twist_spring.h"
#include "../mechanics/robotState.h"
#include "elastic_energy.h"
#include <functional>
#include <vector>
#include <array>

int main() {
    // === Define springs ===
    // Two springs: first uses nodes 0,1,2; second uses nodes 3,4,5
    std::array<int, 5> nodes_edges_index1 = {0, 0, 1, 1, 2}; // [n0, e0, n1, e1, n2]
    std::array<int, 5> nodes_edges_index2 = {3, 2, 4, 3, 5}; // [n0, e0, n1, e1, n2]
    std::array<int, 2> signs = {1, 1};
    std::vector<double> ref_len_arr = {1.0, 1.0, 1.0, 1.0}; // edge lengths (enough for both springs)
    std::vector<double> EI = {1.0, 1.0};          // stiffness
    double GJ = 0.0;

    // Dummy mapping functions
    auto map_node_to_dof = [](const std::vector<int>& v) {
        std::vector<std::array<int, 3>> out;
        for (size_t i = 0; i + 2 < v.size(); i += 3) {
            out.push_back({v[i], v[i+1], v[i+2]});
        }
        return out;
    };
    auto map_edge_to_dof = [](int e) { return e; };

    BendTwistSpring spring1(
        nodes_edges_index1,
        signs,
        ref_len_arr,
        EI,
        GJ,
        map_node_to_dof,
        map_edge_to_dof
    );
    BendTwistSpring spring2(
        nodes_edges_index2,
        signs,
        ref_len_arr,
        EI,
        GJ,
        map_node_to_dof,
        map_edge_to_dof
    );
    std::vector<BendTwistSpring> springs = {spring1, spring2};

    // === Initial robot state ===
    // 6 nodes, 3 DOFs each = 18 DOFs
    Eigen::VectorXd q0(18);
    q0 << 1,2,3, 4,5,6, 7,8,9, 10, 11, 12, 13, 14, 15, 16, 17, 18; // 3 nodes, each with 3 coordinates
    Eigen::VectorXd u = q0;
    std::vector<Eigen::Vector3d> a(6, Eigen::Vector3d::Zero());
    std::vector<Eigen::Vector3d> a1(6, Eigen::Vector3d::Zero());
    std::vector<Eigen::Vector3d> a2(6, Eigen::Vector3d::Zero());
    std::vector<Eigen::Vector3d> m1(6, Eigen::Vector3d::UnitY());
    std::vector<Eigen::Vector3d> m2(6, Eigen::Vector3d::UnitZ());
    Eigen::VectorXi ref_twist = Eigen::VectorXi::Zero(6);
    std::vector<Eigen::Vector3d> tau(6, Eigen::Vector3d::Zero());
    Eigen::VectorXi free_dof = Eigen::VectorXi::LinSpaced(18, 0, 17);

    RobotState state(q0, u, a, a1, a2, m1, m2, ref_twist, tau, free_dof);

    // === Create BendEnergy object ===
    BendEnergy bend_energy(springs, state,
        [](const Eigen::MatrixXd&) {
            return Eigen::MatrixXd(); // fallback to internal strain computation
        });

    // Eigen::VectorXd q(18);
    // q << 1,2,3, 4,5,6, 7,8,9, 10, 11, 12, 13, 14, 15, 16, 17, 18; // 3 nodes, each with 3 coordinates
    // int N_nodes = 6;
    // // M = 3
    // // N = 2
    // auto node_pos = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>>(q.data(), N_nodes, 3);
    // std::cout << "node_pos:\n" << node_pos << "\n\n";

    // === Run get_strain ===
    Eigen::MatrixXd strain = bend_energy.get_strain(state);
    std::cout << "Strain:\n" << strain << "\n\n";

    // // === Run grad_hess_strain ===
    auto [grad, hess] = bend_energy.grad_hess_strain(state);
    for (int i = 0; i < springs.size(); ++i) {
        std::cout << "Gradient of strain (row " << i << "):\n" << grad.row(i) << "\n\n";
        std::cout << "Hessian block (" << i << "):\n" << hess[i] << "\n";
    }

    std::cout << "Gradient shape: " << grad.rows() << " x " << grad.cols() << std::endl;
    std::cout << "Hessian shape: " << hess.size() << " x " << hess[0].rows() << " x " << hess[0].cols() << std::endl;

    return 0;
}