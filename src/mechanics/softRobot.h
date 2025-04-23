#ifndef SOFT_ROBOT_H
#define SOFT_ROBOT_H

#include "../eigenIncludes.h"
#include "../mechanics/robotState.h"
#include "../mechanics/params.cpp"
#include "../mechanics/environment.cpp"
#include "../mechanics/geometry.h"


class SoftRobot {
    public:
        SoftRobot(const GeomParams& geom, const Material& material, const Geometry&geo,
        const SimParams& sim_params, const Environment& env, const RobotState& state);

        void scale_mass_matrix(const Eigen::VectorXi& nodes, double scale);
        void _init_curvature_midedge(const Geometry& geo);

        std::shared_ptr<SoftRobot> update(
            std::optional<Eigen::VectorXd> q = std::nullopt,
            std::optional<Eigen::VectorXd> u = std::nullopt,
            std::optional<Eigen::VectorXd> a = std::nullopt,
            std::optional<Eigen::VectorXd> a1 = std::nullopt,
            std::optional<Eigen::VectorXd> a2 = std::nullopt,
            std::optional<Eigen::VectorXd> m1 = std::nullopt,
            std::optional<Eigen::VectorXd> m2 = std::nullopt,
            std::optional<Eigen::VectorXd> ref_twist = std::nullopt,
            std::optional<Eigen::VectorXi> free_dof = std::nullopt
        ) const;

        Eigen::MatrixXd _compute_tangent(const Eigen::VectorXd& q);
        
        std::shared_ptr<SoftRobot> _fix_dof(const Eigen::VectorXi& new_fixed_dof) const;
        Eigen::VectorXi _get_intermediate_edge_dof(const Eigen::VectorXi& nodes) const;
        
        // Utility Helpers
        Eigen::VectorXi setdiff1d(const Eigen::VectorXi& a, const Eigen::VectorXi& b);
        Eigen::VectorXi union1d(const Eigen::VectorXi& a, const Eigen::VectorXi& b);
        std::vector<bool> isin(const Eigen::VectorXi& a, const Eigen::VectorXi& b);

        // Perturb system
        std::shared_ptr<SoftRobot> move_nodes(const std::vector<int>& nodes, const Eigen::VectorXd& perturbation, std::optional<int> axis) const;
        std::shared_ptr<SoftRobot> twist_edges(const std::vector<int>& edges, const Eigen::VectorXd& perturbation) const;
    
        Eigen::VectorXi map_edge_to_dof(const std::vector<int>& edge_nums) const;
        Eigen::VectorXi map_face_edge_to_dof(const std::vector<int>& edge_nums) const;

        const Eigen::VectorXd& ref_len() const { return ref_len; }
        const Eigen::VectorXd& voronoi_area() const { return voronoi_area; }
        const Eigen::VectorXd& face_area() const { return face_area; }
        const Eigen::MatrixXd& mass_matrix() const { return mass_matrix; }

        // Getter for node_dof_indices (Node DOF indices matrix)
        Eigen::MatrixXi get_node_dof_indices() const {
            Eigen::MatrixXi indices(n_nodes_, 3);
            for (int i = 0; i < n_nodes_; ++i) {
                indices(i, 0) = 3 * i;
                indices(i, 1) = 3 * i + 1;
                indices(i, 2) = 3 * i + 2;
            }
            return indices;
        }

        // Getter for end_node_dof_index (First edge DOF index after node DOFs)
        int get_end_node_dof_index() const {
            return 3 * n_nodes_;
        }

        // Getter for q0 (Initial state vector)
        const Eigen::VectorXd& get_q0() const {
            return q0_;
        }

        // Getter for state (Current state)
        const RobotState& get_state() const {
            return state_;
        }

        // Getter for sim_params
        const SimParams& get_sim_params() const {
            return sim_params_;
        }

        // Getter for env
        const Environment& get_env() const {
            return env_;
        }

        // Getter for n_dof (Total number of degrees of freedom)
        int get_n_dof() const {
            return n_dof_;
        }

        // Getter for bend_twist_springs (List of bend-twist spring elements)
        const std::vector<BendTwistSpring>& get_bend_twist_springs() const {
            return bend_twist_springs_;
        }

        // Getter for stretch_springs (List of stretch spring elements)
        const std::vector<StretchSpring>& get_stretch_springs() const {
            return stretch_springs_;
        }

        // Getter for hinge_springs (List of hinge spring elements)
        const std::vector<HingeSpring>& get_hinge_springs() const {
            return hinge_springs_;
        }

        // Getter for triangle_springs (List of triangle spring elements)
        const std::vector<TriangleSpring>& get_triangle_springs() const {
            return triangle_springs_;
        }

        // Getter for nodes (n_nodes)
        const Eigen::VectorXd& get_nodes() const {
            return nodes_;
        }

        // Getter for edges (n_edges, 2)
        const Eigen::MatrixXi& get_edges() const {
            return edges_;
        }

        // Getter for face_nodes_shell (n_faces, 3)
        const Eigen::MatrixXi& get_face_nodes_shell() const {
            return face_nodes_shell_;
        }

        // Getter for fixed_dof (indices of constrained degrees of freedom)
        Eigen::VectorXi get_fixed_dof() const {
            // np.setdiff1d equivalent in C++ would require a separate implementation
            Eigen::VectorXi fixed_dof = Eigen::VectorXi::Zero(n_dof_ - free_dof_.size());
            int idx = 0;

            for (int i = 0; i < n_dof_; ++i) {
                if (std::find(free_dof_.begin(), free_dof_.end(), i) == free_dof_.end()) {
                    fixed_dof(idx++) = i;
                }
            }

            return fixed_dof.head(idx);  // Resize the result to actual size
        }

         // Fix/free DOF operations
        std::shared_ptr<SoftRobot> free_nodes(const std::vector<int>& nodes, std::optional<int> axis = std::nullopt, bool fix_edges = false) const;
        std::shared_ptr<SoftRobot> fix_nodes(const std::vector<int>& nodes, std::optional<int> axis = std::nullopt, bool fix_edges = false) const;

        std::shared_ptr<SoftRobot> free_edges(const std::vector<int>& edges) const;
        std::shared_ptr<SoftRobot> fix_edges(const std::vector<int>& edges) const;

    private:
        // Parameters
        SimParams sim_params_;
        Environment env_;

        // Geometry
        int n_nodes, n_edges_rod_only, n_edges_shell_only;
        int n_edges_, n_edges_dof_, n_faces_;

        // Eigen::VectorXd nodes_;
        // Eigen::MatrixXi edges;
        // Eigen::MatrixXi face_nodes_shell_;
        // Eigen::VectorXd twist_angles_;

        int n_dof_;
        int n_edges_dof;
        std::vector<int> free_dof_;
        Eigen::VectorXd q0_;
        Eigen::VectorXi fixed_dof;


        // State
        RobotState state_;


        // Mid-edge data
        Eigen::VectorXd tau0_;
        Eigen::MatrixXd init_ts_, init_fs_, init_cs_, init_xis_;

        Eigen::VectorXd ref_len;
        Eigen::VectorXd voronoi_ref_len;
        Eigen::VectorXd voronoi_area;
        Eigen::VectorXd face_area;
        Eigen::VectorXd mass_matrix;

        std::vector<BendTwistSpring> bend_twist_springs_;   // List of bend-twist springs
        std::vector<StretchSpring> stretch_springs_;         // List of stretch springs
        std::vector<HingeSpring> hinge_springs_;             // List of hinge springs
        std::vector<TriangleSpring> triangle_springs_;       // List of triangle springs

        // Helper functions
        void _init_geometry(const Geometry& geo);
        void _init_stiffness(const GeomParams& geom, const Material& material);
        void _init_state(const Geometry& geo);
        void _init_springs(const Geometry& geo);
        void _get_mass_matrix(const GeomParams& geom, const Material& material);

        void _get_ref_len() const;
        void _get_voronoi_ref_len() const;
        void _get_voronoi_area() const;
        void _get_face_area() const;

        std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> init_curvature_midedge(const Geometry& geo);
        Eigen::MatrixXd update_pre_comp_shell(const Eigen::MatrixXd& q);
        std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> compute_space_parallel();
        std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> compute_time_parallel();
        std::pair<Eigen::MatrixXd, Eigen::MatrixXd> compute_material_directors();
        Eigen::VectorXd compute_reference_twist();

}

#endif