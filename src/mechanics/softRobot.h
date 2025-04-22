#ifndef SOFT_ROBOT_H
#define SOFT_ROBOT_H

#include "eigenIncludes.h"

class SoftRobot {
    public:
        SoftRobot(const GeomParams& geom, const Material& material, const Geometry&geo,
        const SimParams& sim_params, const Environment& env);

    private:
        // Parameters
        SimParams sim_params_;
        Environment env_;

        // Geometry
        int n_nodes, n_edges_rod_only, n_edges_shell_only;
        int n_edges_, n_edges_dof_, n_faces_;
        MatrixXd nodes_;
        MatrixXi edges;
        MatrixXi face_nodes_shell_;
        VectorXd twist_angles_;
        int n_dof_;
        VectorXd q0_;

        // Mid-edge data
        VectorXd tau0_;
        MatrixXd init_ts_, init_fs_, init_cs_, init_xis_;

        // Helper functions
        void _init_geometry(const Geometry& geo);
        void _init_stiffness(const GeomParams& geom, const Material& material);
        void _init_state(const Geometry& geo);
        void _init_springs(const Geometry& geo);
        Eigen::VectorXd _get_mass_matrix(const GeomParams& geom, const Material& material);

        Eigen::VectorXd _get_ref_len() const;
        Eigen::VectorXd _get_voronoi_ref_len() const;
        Eigen::VectorXd _get_voronoi_area() const;
        Eigen::VectorXd _get_face_area() const;

        // TODO: Finish adding more helper functions
        Eigen::VectorXd ref_len;
        Eigen::VectorXd voronoi_ref_len;
        Eigen::VectorXd voronoi_area;
        Eigen::VectorXd face_area;
        Eigen::VectorXd mass_matrix;

        int n_nodes;
}

#endif