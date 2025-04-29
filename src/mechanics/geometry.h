#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <Eigen/Dense>
#include "params.cpp"
#include <vector>
#include <string>

#define GEOMETRY_INT int

class Geometry {
public:
    Geometry(const std::vector<Eigen::VectorXd>& nodes,
             const std::vector<std::vector<int>>& edges,
             const std::vector<std::vector<int>>& face_nodes);
             
    Geometry(GeometryData& params); // TODO RADHA

    static Geometry from_txt(const std::string& filename);


    std::vector<Eigen::VectorXd> getNodes() const { return nodes_; }
    std::vector<std::array<int, 2>> getEdges() const { return edges_; }
    std::vector<std::array<int, 2>> getRodEdges() const { return rod_edges_; }
    std::vector<std::array<int, 2>> getShellEdges() const { return shell_edges_; }
    std::vector<std::array<int, 2>> getRodShellJointEdges() const { return rod_shell_joint_edges_; }
    std::vector<std::array<int, 2>> getRodShellJointEdgesTotal() const { return rod_shell_joint_edges_total_; }
    std::vector<std::array<int, 3>> getFaceNodes() const { return face_nodes_; }
    std::vector<std::array<int, 3>> getFaceEdges() const { return face_edges_; }
    std::vector<std::array<int, 3>> getFaceShellEdges() const { return face_shell_edges_; }
    std::vector<std::array<int, 2>> getRodStretchSprings() const { return rod_stretch_springs_; }
    std::vector<std::array<int, 2>> getShellStretchSprings() const { return shell_stretch_springs_; }
    std::vector<std::array<int. 5>> getBendTwistSprings() const { return bend_twist_springs_; }
    std::vector<std::array<int, 2>> getBendTwistSigns() const { return bend_twist_signs_; }
    std::vector<std::array<int, 4>> getHinges() const { return hinges_; }
    std::vector<std::array<int, 3>> getSignFaces() const { return sign_faces_; }
    std::vector<Eigen::VectorXd> getFaceUnitNorms() const { return face_unit_norms_; }
    std::vector<double> getTwistAngles() const { return twist_angles_; }

private:
    std::vector<Eigen::VectorXd> nodes_; // double
    std::vector<std::array<int, 2>> edges_, rod_edges_, shell_edges_; 
    std::vector<std::array<int, 2>> rod_shell_joint_edges_, rod_shell_joint_edges_total_; 
    std::vector<std::array<int, 3>>  face_nodes_, face_edges_, face_shell_edges_;
    std::vector<std::array<int, 2>> rod_stretch_springs_, shell_stretch_springs_;
    std::vector<std::array<int, 5>> bend_twist_springs_;
    std::vector<std::vector<int, 2>> bend_twist_signs_;
    std::vector<std::array<int, 4>> hinges_;
    std::vector<std::array<int, 3>> sign_faces_;
    std::vector<Eigen::VectorXd> face_unit_norms_; // double
    std::vector<double> twist_angles_; //double

    int n_nodes, n_rod_edges, n_rod_shell_joints, n_faces;
    
    //shell
    int n_edges;
    
    // Track number of shell/hinge edges
    int s_i_ = 0, h_i_ = 0;

    // Methods
    void calculate_shell_and_hinge_edges(int n_faces);
    void trim_unused_values();
    void calculate_ghost_edges(int n_rod_shell_joints, int n_faces);
    void calculate_bend_twist_springs(int n_nodes);
    std::vector<std::pair<int, int>> get_combinations(const std::vector<int>& vec);
    std::vector<std::pair<int, int>> get_combinations_outof_into(const std::vector<int>& into, const std::vector<int>& outof);
    void sequence_edges(int n_faces, const Eigen::MatrixXi& face_nodes);
    void create_stretch_springs();
    void calculate_face_edges(const Eigen::MatrixXi& face_nodes);
    int find_edge_index(int n1, int n2);
    void calculate_twist_angles();
    
    static void process_temp_array(int cur_h, const std::vector<std::string>& temp_array, std::vector<Eigen::MatrixXd>& params, const std::vector<int>& h_dtype);
    static std::string trim(const std::string& str);
    static std::vector<std::string> split(const std::string& str, const std::string& delimiter);

    static Eigen::MatrixXd safe_concat(const Eigen::MatrixXd& arr1, const Eigen::MatrixXd& arr2);

    static std::pair<Eigen::MatrixXi, Eigen::MatrixXi> separate_joint_edges(
        const Eigen::MatrixXi& triangles,
        const Eigen::MatrixXi& edges);
    
};

#endif
