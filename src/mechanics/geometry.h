#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <Eigen/Dense>
#include <vector>
#include <string>

class Geometry {
public:
    Geometry(const Eigen::MatrixXd& nodes,
             const Eigen::MatrixXi& edges,
             const Eigen::MatrixXi& face_nodes);

    static Geometry from_txt(const std::string& filename);

    const Eigen::MatrixXd& nodes() const;
    const Eigen::MatrixXi& edges() const;
    const Eigen::MatrixXi& rod_edges() const;
    const Eigen::MatrixXi& shell_edges() const;
    const Eigen::MatrixXi& rod_shell_joint_edges() const;
    const Eigen::MatrixXi& rod_shell_joint_edges_total() const;
    const Eigen::MatrixXi& face_nodes() const;
    const Eigen::MatrixXi& face_edges() const;
    const Eigen::MatrixXi& face_shell_edges() const;
    const Eigen::MatrixXi& rod_stretch_springs() const;
    const Eigen::MatrixXi& shell_stretch_springs() const;
    const Eigen::MatrixXi& bend_twist_springs() const;
    const Eigen::MatrixXi& bend_twist_signs() const;
    const Eigen::MatrixXi& hinges() const;
    const Eigen::MatrixXi& sign_faces() const;
    const Eigen::MatrixXd& face_unit_norms() const;
    const Eigen::VectorXd& twist_angles() const;

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

    static std::pair<Eigen::MatrixXi, Eigen::MatrixXi> separate_joint_edges(
        const Eigen::MatrixXi& triangles,
        const Eigen::MatrixXi& edges);

    static Eigen::MatrixXi safe_concat(
        const Eigen::MatrixXi& arr1,
        const Eigen::MatrixXi& arr2);
};

#endif
