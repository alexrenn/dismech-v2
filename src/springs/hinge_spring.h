#ifndef HINGE_SPRING_H
#define HINGE_SPRING_H
#pragma once
#include <vector>
#include <array>
#include <functional>

class HingeSpring {
public:
    double kb;
    std::array<int, 4> nodes_ind;
    std::vector<int> ind;

    HingeSpring(const std::array<int, 4>& nodes_ind, double kb,
                const std::function<std::array<int, 3>(int)>& map_node_to_dof)
        : kb(kb), nodes_ind(nodes_ind)
    {
        std::array<int, 3> dof0 = map_node_to_dof(nodes_ind[0]);
        std::array<int, 3> dof1 = map_node_to_dof(nodes_ind[1]);
        std::array<int, 3> dof2 = map_node_to_dof(nodes_ind[2]);

        ind.reserve(9);
        ind.insert(ind.end(), dof0.begin(), dof0.end());
        ind.insert(ind.end(), dof1.begin(), dof1.end());
        ind.insert(ind.end(), dof2.begin(), dof2.end());
    }
};


#endif