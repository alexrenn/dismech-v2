#ifndef STRETCH_SPRING_H
#define STRETCH_SPRING_H
#include <array>
#include <vector>
#include <functional>

class StretchSpring {
public:
    double EA;                     // Axial stiffness
    double ref_len;                // Rest length
    std::array<int, 2> nodes_ind;  // 2 node indices
    std::vector<int> ind;          // DOF indices

    StretchSpring(
        const std::array<int, 2>& nodes_ind,                    // [n0, n1]
        double ref_len,
        double EA,
        const std::function<std::array<int, 3>(int)>& map_node_to_dof)
        : EA(EA), ref_len(ref_len), nodes_ind(nodes_ind)
    {
        // Each node has 3 DOFs, so 2 nodes → 6 total
        std::array<int, 3> dof0 = map_node_to_dof(nodes_ind[0]);
        std::array<int, 3> dof1 = map_node_to_dof(nodes_ind[1]);

        ind.reserve(6);
        ind.insert(ind.end(), dof0.begin(), dof0.end());
        ind.insert(ind.end(), dof1.begin(), dof1.end());
    }
};
#endif