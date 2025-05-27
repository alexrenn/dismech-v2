#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <memory>

#include "softRobot.h"
#include "robotState.h"

/*
Guide for sparse matricies
https://eigen.tuxfamily.org/dox/group__SparseQuickRefPage.html
*/

class ElasticEnergy {
    public:
        ElasticEnergy(const Eigen::MatrixXd& K,
                      const Eigen::MatrixXi& nodes_ind,
                      const Eigen::MatrixXi& ind,
                      const RobotState& initial_state,
                      std::function<Eigen::VectorXd(const Eigen::MatrixXd&)> get_strain = nullptr)
            : K_(K), ind_(ind), initial_state_(initial_state), get_strain_(get_strain) {
    
            n_K_ = K_.cols();
            node_dof_ind_ = SoftRobot::map_node_to_dof(nodes_ind);
            n_nodes_ = nodes_ind.cols();
            rows_ = create_rows(ind_);
            cols_ = create_cols(ind_);
            
            if (get_strain_) { 
                nat_strain_ = get_strain_(get_node_pos(initial_state_.q));
            }
        }
    
        virtual ~ElasticEnergy() = default;
    
        virtual void post_init() {
            if (!nat_strain_.size()) {
                nat_strain_ = get_strain_(initial_state_.q);
            }
        }

        // TODO: _get_node_pos
        // TODO: set_nat_strain
    
        /*
        After reshaping, K_ will be a matrix with rows representing individual stiffness coefficients
        and columns corresponding to different degrees of freedom (DOF).*/
        virtual Eigen::MatrixXd get_energy_linear_elastic(const RobotState& state, bool output_scalar = true) {
            Eigen::VectorXd strain = get_strain_(state.q);
            Eigen::VectorXd del_strain = (strain - nat_strain_).reshaped(-1, n_K_);
            
            if (output_scalar) {
                return 0.5 * (K_.reshaped(-1, n_K_) * del_strain).squaredNorm(); //squaredNorm is equivalent to np.sum()
            }
            return 0.5 * (K_.reshaped(-1, n_K_) * del_strain).array().square();
        }

         // Abstract method for getting strain (to be implemented by subclasses)
        virtual Eigen::VectorXd getStrain(const RobotState& state) const = 0;

        std::pair<Eigen::VectorXd, Eigen::SparseMatrix<double>> grad_hess_energy_linear_elastic(
            const RobotState& state, bool sparse = false) {
            
            // Compute Gradidents and Hessians of Strain
            Eigen::VectorXd strain = get_strain_(state.q);
            Eigen::MatrixXd grad_strain, hess_strain;
            std::tie(grad_strain, hess_strain) = grad_hess_strain(state);
            
            // Compute difference in strain
            Eigen::VectorXd del_strain = (strain - nat_strain_).reshaped(-1, n_K_);
            // Compute Graident of the energy
            Eigen::VectorXd gradE_strain = (K_.reshaped(-1, n_K_) * del_strain).array();
            
            // Reshaping gradients and Hessians similarly to Python, to handle multiple strain constraints (bending, twisting, etc)
            // Needs to match the expected dimensions of subsequent computations
            gradE_strain = gradE_strain.reshaped(gradE_strain.size(), n_K_);
            grad_strain = grad_strain.reshaped(grad_strain.size(), grad_strain.size(), n_K_);
            hess_strain = hess_strain.reshaped(hess_strain.size(), grad_strain.size(), grad_strain.size(), n_K_);
            
            // Gradient and Hessian terms vectorized
            Eigen::VectorXd grad_energy = gradE_strain.colwise().sum() * grad_strain;
            Eigen::MatrixXd hess_term1 = gradE_strain.transpose()  * hess_strain; 
            Eigen::MatrixXd outer = grad_strain.transpose() * grad_strain;
            Eigen::MatrixXd hess_term2 = K_.reshaped(-1, n_K_) * outer;
            
            // Compute the total Hessian
            Eigen::MatrixXd hess_energy = hess_term1 + hess_term2;
            
            // Compute forces and Hessians
            Eigen::VectorXd Fs = Eigen::VectorXd::Zero(state.q.rows());
            Fs += -grad_energy;
            
            Eigen::SparseMatrix<double> Js(state.q.rows(), state.q.rows());
            if (sparse) {
                Js.setZero();
                // Implement sparse matrix construction
            }
            
            return std::make_pair(Fs, Js);
        }
        
        // TODO : fdm_check_grad_hess_strain
    
    protected:
        Eigen::MatrixXd K_; // Stiffness matrix
        Eigen::MatrixXi ind_; // Indices
        Eigen::VectorXd nat_strain_; // Natural strain
        Eigen::MatrixXi node_dof_ind_; // DOF indices for nodes
        RobotState initial_state_; // Initial robot state
        std::function<Eigen::VectorXd(const Eigen::MatrixXd&)> get_strain_; // Strain function
        Eigen::VectorXi rows_; // Row indices for sparse matrix construction
        Eigen::VectorXi cols_; // Column indices for sparse matrix construction
        int n_K_;
        int n_nodes_;
    
        Eigen::MatrixXd create_rows(const Eigen::MatrixXi& ind) {
            // Implement logic to create row indices for sparse matrix
        }
    
        Eigen::MatrixXd create_cols(const Eigen::MatrixXi& ind) {
            // Implement logic to create column indices for sparse matrix
        }
    
        Eigen::MatrixXd get_node_pos(const Eigen::MatrixXd& q) {
            // Implement the logic to extract node positions
        }
    };
    