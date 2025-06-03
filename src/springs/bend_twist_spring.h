#ifndef BEND_TWIST_SPRING_H
#define BEND_TWIST_SPRING_H

#include <vector>
#include <array>
#include <functional>
#include <numeric>

class BendTwistSpring {
public:
    std::vector<double> stiff_EI;
    double stiff_GJ;

    std::array<int, 3> nodes_ind;   // [n0, n1, n2]
    std::array<int, 2> edges_ind;   // [e0, e1]
    std::array<int, 2> sgn;         // [sign0, sign1]
    std::vector<int> ind;           // DOF indices
    double voronoi_len;

    BendTwistSpring(
        const std::array<int, 5>& nodes_edges_index,            // [n0, e0, n1, e1, n2]
        const std::array<int, 2>& signs,                         // [s0, s1]
        const std::vector<double>& ref_len_arr,                 // edge lengths
        const std::vector<double>& EI,                          // stiffness matrix (vectorized)
        double GJ,                                              // torsional stiffness
        const std::function<std::vector<std::array<int, 3>>(const std::vector<int>&)>& map_node_to_dof,
        const std::function<int(int)>& map_edge_to_dof)
        : stiff_EI(EI), stiff_GJ(GJ), sgn(signs)
    {
        // Extract node and edge indices
        nodes_ind = { nodes_edges_index[0], nodes_edges_index[2], nodes_edges_index[4] };
        edges_ind = { nodes_edges_index[1], nodes_edges_index[3] };

        // Map node DOFs: returns vector of arrays of 3 ints
        std::vector<int> node_ids = { nodes_ind[0], nodes_ind[1], nodes_ind[2] };
        std::vector<std::array<int, 3>> node_dofs = map_node_to_dof(node_ids);

        // Flatten to ind
        ind.reserve(9 + 2); // 3 nodes × 3 DOFs + 2 edge DOFs
        for (const auto& dof_triplet : node_dofs) {
            ind.insert(ind.end(), dof_triplet.begin(), dof_triplet.end());
        }

        // Add edge DOFs
        ind.push_back(map_edge_to_dof(edges_ind[0]));
        ind.push_back(map_edge_to_dof(edges_ind[1]));

        // Compute Voronoi length
        voronoi_len = 0.5 * (ref_len_arr[edges_ind[0]] + ref_len_arr[edges_ind[1]]);
    }
};

#endif