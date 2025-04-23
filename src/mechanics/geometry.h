#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <Eigen/Dense>
#include <vector>
#include <string>

#define GEOMETRY_INT int

class Geometry {
public:
    Geometry(const Eigen::MatrixXd& nodes,
             const Eigen::MatrixXi& edges,
             const Eigen::MatrixXi& face_nodes);
    Geometry(const std::vector<Eigen::MatrixXd>& params);

    static Geometry from_txt(const std::string& filename);


    Eigen::MatrixXd getNodes() const { return nodes_; }
    Eigen::MatrixXi getEdges() const { return edges_; }
    Eigen::MatrixXi getRodEdges() const { return rod_edges_; }
    Eigen::MatrixXi getShellEdges() const { return shell_edges_; }
    Eigen::MatrixXi getRodShellJointEdges() const { return rod_shell_joint_edges_; }
    Eigen::MatrixXi getRodShellJointEdgesTotal() const { return rod_shell_joint_edges_total_; }
    Eigen::MatrixXi getFaceNodes() const { return face_nodes_; }
    Eigen::MatrixXi getFaceEdges() const { return face_edges_; }
    Eigen::MatrixXi getFaceShellEdges() const { return face_shell_edges_; }
    Eigen::MatrixXi getRodStretchSprings() const { return rod_stretch_springs_; }
    Eigen::MatrixXi getShellStretchSprings() const { return shell_stretch_springs_; }
    Eigen::MatrixXi getBendTwistSprings() const { return bend_twist_springs_; }
    Eigen::MatrixXi getBendTwistSigns() const { return bend_twist_signs_; }
    Eigen::MatrixXi getHinges() const { return hinges_; }
    Eigen::MatrixXi getSignFaces() const { return sign_faces_; }
    Eigen::MatrixXd getFaceUnitNorms() const { return face_unit_norms_; }
    Eigen::VectorXd getTwistAngles() const { return twist_angles_; }

private:
    Eigen::MatrixXd nodes_;
    Eigen::MatrixXi edges_, rod_edges_, shell_edges_;
    Eigen::MatrixXi rod_shell_joint_edges_, rod_shell_joint_edges_total_;
    Eigen::MatrixXi face_nodes_, face_edges_, face_shell_edges_;
    Eigen::MatrixXi rod_stretch_springs_, shell_stretch_springs_;
    Eigen::MatrixXi bend_twist_springs_, bend_twist_signs_;
    Eigen::MatrixXi hinges_, sign_faces_;
    Eigen::MatrixXd face_unit_norms_;
    Eigen::VectorXd twist_angles_;

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
