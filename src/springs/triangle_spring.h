#ifndef TRIANGLE_SPRING_H
#define TRIANGLE_SPRING_H

#include <array>
#include <vector>
#include <functional>
#include "../eigenIncludes.h"

class TriangleSpring {
public:
    double kb;
    double nu;

    std::array<int, 3> nodes_ind;
    std::array<int, 3> edges_ind;
    std::array<int, 3> face_edges;
    std::array<int, 3> sgn;

    std::array<double, 3> ref_len;

    std::vector<double> A;
    Eigen::Matrix3d init_ts;
    std::array<double, 3> init_fs;
    std::array<double, 3> init_cs;
    std::array<double, 3> init_xis;

    std::vector<int> ind;

    TriangleSpring(
        const std::array<int, 3>& nodes_ind,
        const std::array<int, 3>& edges_ind,
        const std::array<int, 3>& face_edges,
        const std::array<int, 3>& signs,
        const std::vector<double>& ref_len_arr,
        const std::vector<double>& A, // RADHA why in a vector
        const Eigen::Matrix3d& init_ts,
        const std::array<double, 3> & init_fs,
        const std::array<double, 3> & init_cs,
        const std::array<double, 3> & init_xis,
        double kb,
        double nu,
        const std::function<std::array<int, 3>(int)>& map_node_to_dof,
        const std::function<int(int)>& map_face_edge_to_dof)
        : kb(kb), nu(nu),
          nodes_ind(nodes_ind),
          edges_ind(edges_ind),
          face_edges(face_edges),
          sgn(signs),
          A(A), init_ts(init_ts), init_fs(init_fs),
          init_cs(init_cs), init_xis(init_xis)
    {
        // Select 3 reference lengths using face_edges as indices
        for (int i = 0; i < 3; ++i) {
            ref_len[i] = ref_len_arr[face_edges[i]];
        }

        // DOF mapping
        std::array<int, 3> dof0 = map_node_to_dof(nodes_ind[0]);
        std::array<int, 3> dof1 = map_node_to_dof(nodes_ind[1]);
        std::array<int, 3> dof2 = map_node_to_dof(nodes_ind[2]);

        int edge_dof0 = map_face_edge_to_dof(edges_ind[0]);
        int edge_dof1 = map_face_edge_to_dof(edges_ind[1]);
        int edge_dof2 = map_face_edge_to_dof(edges_ind[2]);

        ind.reserve(9 + 3); // 3 nodes × 3 DOFs + 3 edge DOFs
        ind.insert(ind.end(), dof0.begin(), dof0.end());
        ind.insert(ind.end(), dof1.begin(), dof1.end());
        ind.insert(ind.end(), dof2.begin(), dof2.end());
        ind.push_back(edge_dof0);
        ind.push_back(edge_dof1);
        ind.push_back(edge_dof2);
    }
};


#endif