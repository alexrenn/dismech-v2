#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "../eigenIncludes.h"
#include "params.cpp"
#include <vector>
#include <string>

#define GEOMETRY_INT int

class Geometry {
public:
    Geometry(const std::string& fname);
             
    void from_txt(const std::string& filename);


    std::vector<Eigen::Vector3d> getNodes() const { return nodes_; }
    std::vector<std::array<int, 2> > getEdges() const { return edges_; }
    std::vector<std::array<int, 2> > getRodEdges() const { return rod_edges_; }
    std::vector<std::array<int, 2> > getShellEdges() const { return shell_edges_; }
    std::vector<std::array<int, 2>> getRodShellJointEdges() const { return rod_shell_joint_edges_; }
    std::vector<std::array<int, 2> > getRodShellJointEdgesTotal() const { return rod_shell_joint_edges_total_; }
    std::vector<std::array<int, 3> > getFaceNodes() const { return face_nodes_; }
    std::vector<std::array<int, 3> > getFaceEdges() const { return face_edges_; }
    std::vector<std::array<int, 3> > getFaceShellEdges() const { return face_shell_edges_; }
    std::vector<std::array<int, 2> > getRodStretchSprings() const { return rod_stretch_springs_; }
    std::vector<std::array<int, 2> > getShellStretchSprings() const { return shell_stretch_springs_; }
    std::vector<std::array<int, 5> > getBendTwistSprings() const { return bend_twist_springs_; }
    std::vector<std::array<int, 2> > getBendTwistSigns() const { return bend_twist_signs_; }
    std::vector<std::array<int, 4> > getHinges() const { return hinges_; }
    std::vector<std::array<int, 3> > getSignFaces() const { return sign_faces_; }
    std::vector<Eigen::Vector3d> getFaceUnitNorms() const { return face_unit_norms_; }
    std::vector<double> getTwistAngles() const { return twist_angles_; }

private:
    std::vector<Eigen::Vector3d> nodes_; // double
    std::vector<std::array<int, 2> > edges_, rod_edges_, shell_edges_; 
    std::vector<std::array<int, 2> > rod_shell_joint_edges_, rod_shell_joint_edges_total_; 
    std::vector<std::array<int, 3> >  face_nodes_, face_edges_, face_shell_edges_;
    std::vector<std::array<int, 2> > rod_stretch_springs_, shell_stretch_springs_;
    std::vector<std::array<int, 5> > bend_twist_springs_;
    std::vector<std::array<int, 2> > bend_twist_signs_;
    std::vector<std::array<int, 4> > hinges_;
    std::vector<std::array<int, 3> > sign_faces_;
    std::vector<Eigen::Vector3d> face_unit_norms_; // double
    std::vector<double> twist_angles_; //double
    std::vector<double> As_;
    std::vector<std::array<int, 2> > edge_faces_;
    std::vector<int> third_node_;
    std::vector<std::array<int, 2>> ghost_rod_shell_joint_edges_;

    int n_nodes, n_rod_edges, n_rod_shell_joints, n_faces, n_edges;;

    
    // Track number of shell/hinge edges
    int s_i_ = 0, h_i_ = 0;

    // Methods
    void calculate_shell_and_hinge_edges(int n_faces);
    void trim_unused_values();
    void calculate_ghost_edges(int n_rod_shell_joints, int n_faces);
    void calculate_bend_twist_springs(int n_nodes);
    std::vector<std::array<int, 2> > get_combinations(const std::vector<int>& vec);
    std::vector<std::array<int, 2> >  get_combinations_outof_into(const std::vector<int>& into, const std::vector<int>& outof);
    void sequence_edges(int n_faces, std::vector<std::array<int, 3> >& face_nodes);
    void create_stretch_springs();
    void calculate_face_edges(std::vector<std::array<int, 3> >& face_nodes);
    int find_edge_index(int n1, int n2);
    void calculate_twist_angles();
    
    static void process_temp(int header,
        std::vector<std::vector<double>>& temp,
        std::vector<Eigen::Vector3d>& nodes,
        std::vector<std::array<int, 2> >& edges,
        std::vector<std::array<int, 3> >& triangles);

    static std::pair<std::vector<std::array<int, 2>>, std::vector<std::array<int, 2>>> separate_joint_edges( 
        const std::vector<std::array<int, 3>>& triangles,
        const std::vector<std::array<int, 2>>& edges);
};

#endif
