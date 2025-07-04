#ifndef SOFT_ROBOT_H
#define SOFT_ROBOT_H

#include "../eigenIncludes.h"
#include "robotState.h"
#include "environment.h"
#include "geometry.h"
#include "../frame_util.h"
#include "../springs/bend_twist_spring.h"
#include "../springs/hinge_spring.h"
#include "../springs/stretch_spring.h"
#include "../springs/triangle_spring.h"
#include "stiffness.h"


class SoftRobot {
    public:
        SoftRobot(const GeomParams& geom, const Material& material, const Geometry&geo,
        const SimParams& sim_params, const Environment& env, const RobotState& state);

        void scale_mass_matrix(const std::vector<Eigen::Vector3d> & nodes, double scale); 
        std::tuple<
            std::vector<Eigen::Matrix3d>,
            std::vector<std::array<double, 3>>,
            std::vector<std::array<double, 3>>,
            std::vector<std::array<double, 3>>
        > _init_curvature_midedge(const Geometry& geo);

        std::vector<Eigen::Vector3d> _compute_tangent(const Eigen::VectorXd& q);
        
        void _fix_dof(const Eigen::VectorXi& new_fixed_dof);
        Eigen::VectorXi _get_intermediate_edge_dof(const std::vector<int>& nodes);
        
        // Utility Helpers
        Eigen::VectorXi setdiff1d(const Eigen::VectorXi& a, const Eigen::VectorXi& b);
        Eigen::VectorXi union1d(const Eigen::VectorXi& a, const Eigen::VectorXi& b);
        std::vector<bool> isin(const Eigen::VectorXi& a, const Eigen::VectorXi& b);

        // Perturb system
        void move_nodes(const std::vector<int>& nodes, const Eigen::VectorXd& perturbation, std::optional<int> axis);
        void twist_edges(const std::vector<int>& edge_indices, const Eigen::VectorXd& perturbation);
    
        Eigen::VectorXi _get_node_dof_mask(const std::vector<int>& nodes, std::optional<int> axis = std::nullopt) const;
        static std::vector<std::array<int, 3>> map_node_to_dof(const std::vector<int>& node_indices);
        Eigen::VectorXi map_edge_to_dof(const std::vector<int>& edge_indices) const;  
        Eigen::VectorXi map_face_edge_to_dof(const std::vector<int>& edge_nums) const;

        //simple getters
        std::vector<double>& get_ref_len()  { return ref_len; }
        std::vector<double>& get_voronoi_ref_len() { return voronoi_ref_len; }
        std::vector<double>& get_voronoi_area() { return voronoi_area; }
        std::vector<double>&  get_face_area() { return face_area; }
        const Eigen::MatrixXd& get_mass_matrix() { return mass_matrix; }
        
        std::vector<std::array<int, 3>> get_node_dof_indices() const;
        Eigen::VectorXi get_fixed_dof() const;

        // Getter for end_node_dof_index (First edge DOF index after node DOFs)
        int get_end_node_dof_index() const { return 3 * n_nodes_;}

        // Getter for q0 (Initial state vector)
        const Eigen::VectorXd& get_q0() const { return q0_; }

        // Getter for state (Current state)
        const RobotState& get_state() const { return state_; }

        // Getter for sim_params
        const SimParams& get_sim_params() const { return sim_params_; }

        // Getter for env
        const Environment& get_env() const { return env_;}

        // Getter for n_dof (Total number of degrees of freedom)
        int get_n_dof() const { return n_dof_; }

       
        // // Getter for bend_twist_springs (List of bend-twist spring elements)
        // const std::vector<BendTwistSpring>& get_bend_twist_springs() const { return bend_twist_springs_; }

        // // Getter for stretch_springs (List of stretch spring elements)
        // const std::vector<StretchSpring>& get_stretch_springs() const { return stretch_springs_; }

        // // Getter for hinge_springs (List of hinge spring elements)
        // const std::vector<HingeSpring>& get_hinge_springs() const { return hinge_springs_; }

        // // Getter for triangle_springs (List of triangle spring elements)
        // const std::vector<TriangleSpring>& get_triangle_springs() const { return triangle_springs_; } 

        // Getter for nodes (n_nodes)
        const  std::vector<Eigen::Vector3d>& get_nodes() { return nodes_; }

        // Getter for edges (n_edges, 2)
        const  std::vector<std::array<int, 2>>& get_edges() { return edges_; }

        // Getter for face_nodes_shell (n_faces, 3)
        const std::vector<std::array<int, 3>>& get_face_nodes_shell() { return face_nodes_shell_; }

        std::vector<Eigen::Vector3d> update_pre_comp_shell(const Eigen::VectorXd& q);

         // Fix/free DOF operations
        void free_nodes(const std::vector<int>& nodes, std::optional<int> axis, bool fix_edges);
        void fix_nodes(const std::vector<int>& nodes, std::optional<int> axis, bool fix_edges);

        void free_edges(const std::vector<int>& edge_indices);
        void fix_edges(const std::vector<int>& edge_indices);

        // Debug function
        void debug();

    private:
        // Parameters
        SimParams sim_params_;
        Environment env_;
        double _TANGENT_THRESHOLD = 1e-10;

        // Geometry
        int n_nodes_, n_edges_rod_only, n_edges_shell_only;
        int n_edges_, n_edges_dof_, n_faces_;

        // types
        std::vector<Eigen::Vector3d> nodes_;
        std::vector<std::array<int, 3> > face_nodes_;
        std::vector<std::array<int, 2>> edges_;
        std::vector<std::array<int, 3>> face_nodes_shell_;
        std::vector<std::array<int, 3> > face_shell_edges_;
        std::vector<std::array<int, 3>> face_edges_;
        std::vector<std::array<int, 3>> sign_faces_;
        std::vector<double>  twist_angles_;
        std::vector<std::array<int, 2> > sign_;
        std::vector<std::array<int, 4> > hinges_;

        int n_dof_;
        Eigen::VectorXi  free_dof_;
        Eigen::VectorXd q0_;
        Eigen::VectorXi fixed_dof_;
        std::vector<Eigen::Vector3d> tangent_;


        // State
        RobotState state_;

        // Mid-edge data
        std::vector<Eigen::Vector3d> tau0_; 
        std::vector<Eigen::Matrix3d> init_ts_;
        std::vector<std::array<double, 3>> init_fs_;
        std::vector<std::array<double, 3>> init_cs_;
        std::vector<std::array<double, 3>> init_xis_;

        std::vector<double> ref_len;
        std::vector<double> voronoi_ref_len;
        std::vector<double>  voronoi_area;
        std::vector<double> face_area;
        Eigen::MatrixXd mass_matrix;

        // Rod stiffness variables
        double EA, EI1, EI2, GJ;

        // Shell stiffness variables
        std::vector<double> ks; // Stiffness for each shell face
        double kb; // Bending stiffness for shell faces
        double nu; // Poisson's ratio for shell faces

        std::vector<std::array<int, 5> > bend_twist_springs_;   // List of bend-twist springs
        std::vector<std::array<int, 2> > rod_stretch_springs_;         // List of stretch springs
        std::vector<std::array<int, 2> > shell_stretch_springs_;
        std::vector<std::array<int, 4> > hinge_springs_;             // List of hinge springs
        std::vector<std::array<int, 2> > triangle_springs_;       // List of triangle springs

        std::vector<std::array<int, 2> >  bend_twist_signs_;

        std::vector<StretchSpring> rod_stretch_springs_object_;
        std::vector<StretchSpring> shell_stretch_springs_object_;
        std::vector<BendTwistSpring> bend_twist_springs_object_;
        std::vector<TriangleSpring> triangle_springs_object_;
        std::vector<HingeSpring> hinge_springs_object_;
        
        std::vector<StretchSpring> stretch_springs_;



        // Helper functions
        void _init_geometry(const Geometry& geo);
        void _init_stiffness(const GeomParams& geom, const Material& material);
        void _init_state(const Geometry& geo);
        void _init_springs(const Geometry& geo);
        void _get_mass_matrix(const GeomParams& geom, const Material& material);

        void _get_ref_len();
        void _get_voronoi_ref_len();
        void _get_voronoi_area();
        void _get_face_area();

        std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> compute_space_parallel(); // TODO check 
        std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> compute_time_parallel(
                                                                                const std::vector<Eigen::Vector3d>& a1_old,
                                                                                const Eigen::VectorXd& q0,
                                                                                const Eigen::VectorXd& q);
        std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> compute_material_directors(
                                                                                const Eigen::VectorXd& q,
                                                                                const std::vector<Eigen::Vector3d>& a1,
                                                                                const std::vector<Eigen::Vector3d>& a2); 
        Eigen::VectorXi compute_reference_twist(const std::vector<BendTwistSpring>& springs,
                                                const Eigen::VectorXd& q,
                                                const std::vector<Eigen::Vector3d>& a1,
                                                const Eigen::VectorXi& ref_twist);

    
        

};

#endif