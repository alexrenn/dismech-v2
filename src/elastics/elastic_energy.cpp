#include "elastic_energy.h"

#include <numeric>
#include <iostream>

ElasticEnergy::ElasticEnergy(const Eigen::MatrixXd& K,
                             const std::vector<std::vector<int>>& nodes_ind,
                             const std::vector<std::vector<int>>& ind,
                             const RobotState& initial_state,
                             std::function<Eigen::MatrixXd(const Eigen::MatrixXd&)> get_strain)
    : _K(K),
      _ind(ind),
      _initial_state(initial_state)
{
    _n_K = (_K.cols() == 1) ? 1 : _K.cols();
    _nodes_ind = nodes_ind;
    _n_nodes = nodes_ind[0].size(); // in python this is nodes_ind.shape[1]

    std::vector<int> flat_node_indices;
    _n_elems = nodes_ind.size();

    for (int col = 0; col < _n_nodes; ++col) {
        for (int row = 0; row < _n_elems; ++row) {
            flat_node_indices.push_back(nodes_ind[row][col]);
        }
    }
    auto dof_inds = SoftRobot::map_node_to_dof(flat_node_indices);
    for (const auto& arr : dof_inds) {
        _node_dof_ind.insert(_node_dof_ind.end(), arr.begin(), arr.end());
    }

    // Build _rows and _cols for sparse indexing
    int stencil_n_dof = static_cast<int>(_ind[0].size());
    for (const auto& idx_row : _ind) {
        for (int i = 0; i < stencil_n_dof; ++i) {
            for (int j = 0; j < stencil_n_dof; ++j) {
                _rows.push_back(idx_row[i]);
                _cols.push_back(idx_row[j]);
            }
        }
    }

    // Initialize natural strain
    // if (get_strain) {
    //     Eigen::MatrixXd positions = get_node_positions(_initial_state.q0);
    //     _nat_strain = get_strain(positions);
    // } else {
    //     _nat_strain = Eigen::MatrixXd(); // empty, to be handled later
    // }

    // compute_post_init_strain(get_strain);
}

void ElasticEnergy::compute_post_init_strain(std::function<Eigen::MatrixXd(const Eigen::MatrixXd&)> get_strain)
{
    if (_nat_strain.size() == 0) {
        _nat_strain = this->get_strain(_initial_state);
    } else {
        Eigen::MatrixXd fallback = this->get_strain(_initial_state);
        for (int i = 0; i < _nat_strain.rows(); ++i) {
            for (int j = 0; j < _nat_strain.cols(); ++j) {
                if (std::isnan(_nat_strain(i, j))) {
                    _nat_strain(i, j) = fallback(i, j);
                }
            }
        }
    }
}

//    Eigen::MatrixXd ElasticEnergy::get_node_positions(const Eigen::VectorXd& q) const
//    {
//        int N_nodes = _n_nodes;
//     //    return Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>>(q.data(), N_nodes, 3);
//         // return Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>>(q[node_dof_ind], N_nodes, 3);
//     //     Eigen::MatrixXd positions(N_nodes, 3);
//     //     for (int i = 0; i < N_nodes; ++i) {
//     //         positions.row(i) = q.segment<3>(3 * _node_dof_ind[i]);
//     // }
//     // M * N * 3
//     return Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>>(positions.data(), N_nodes, 3);
//    }

double ElasticEnergy::get_energy_linear_elastic_scalar(const RobotState& state) const { // The return type is double because Python's output_scalar=True is the default.
    std::cout << "I AM IN SCALAR VERSION" << std::endl;
    // 1. Get strain
    Eigen::MatrixXd strain = this->get_strain(state);

    // 2. Compute delta strain (difference from natural strain)
    Eigen::MatrixXd del_strain = strain - _nat_strain;

    // 3. Reshape _K and del_strain to make sure their shapes align
    Eigen::MatrixXd K_reshaped = _K;
    if (_K.cols() == 1) {
        K_reshaped = _K.replicate(1, del_strain.cols());
    }

    std::cout << "K_reshaped.rows(): " << K_reshaped.rows() << ", del_strain.rows(): " << del_strain.rows() << std::endl;
    std::cout << "K_reshaped.cols(): " << K_reshaped.cols() << ", del_strain.cols(): " << del_strain.cols() << std::endl;
    assert(K_reshaped.rows() == del_strain.rows() && K_reshaped.cols() == del_strain.cols());

    // 4. Energy calculation
    return 0.5 * (K_reshaped.array() * del_strain.array().square()).sum();
}

Eigen::MatrixXd ElasticEnergy::get_energy_linear_elastic_matrix(const RobotState& state) const {
        // 1. Get strain
    std::cout << "I AM IN VECTOR VERSION" << std::endl;
    Eigen::MatrixXd strain = this->get_strain(state);

    // 2. Compute delta strain (difference from natural strain)
    Eigen::MatrixXd del_strain = strain - _nat_strain;

    // 3. Reshape _K and del_strain to make sure their shapes align
    Eigen::MatrixXd K_reshaped = _K;
    if (_K.cols() == 1) {
        K_reshaped = _K.replicate(1, del_strain.cols());
    }

    std::cout << "K_reshaped.rows(): " << K_reshaped.rows() << ", del_strain.rows(): " << del_strain.rows() << std::endl;
    std::cout << "K_reshaped.cols(): " << K_reshaped.cols() << ", del_strain.cols(): " << del_strain.cols() << std::endl;
    assert(K_reshaped.rows() == del_strain.rows() && K_reshaped.cols() == del_strain.cols());

    Eigen::MatrixXd energy_mat = 0.5 * K_reshaped.array() * del_strain.array().square();
    return energy_mat;  // If needed, change return type to Eigen::MatrixXd
}



std::pair<Eigen::VectorXd, Eigen::MatrixXd> 
ElasticEnergy::grad_hess_energy_linear_elastic_dense(const RobotState& state) const {
    const int n_elements = static_cast<int>(_ind.size());
    const int n_dof = state.q0.size();
    const int stencil_dof = static_cast<int>(_ind[0].size());

    // === 1. Compute strain, gradient, Hessian of strain ===
    Eigen::MatrixXd strain = this->get_strain(state);
    auto [grad_strain, hess_strain] = this->grad_hess_strain(state);

    Eigen::MatrixXd del_strain = strain - _nat_strain;
    Eigen::MatrixXd gradE_strain = _K.array() * del_strain.array();  // (n x k)

    // gradE_strain = gradE_strain.reshape(gradE_strain.rows(), _n_K);

    // === 2. Compute gradient ===
    Eigen::MatrixXd grad_energy = Eigen::MatrixXd::Zero(n_elements, stencil_dof);
    for (int e = 0; e < n_elements; ++e) {
        for (int i = 0; i < stencil_dof; ++i) {
            for (int k = 0; k < _n_K; ++k) {
                grad_energy(e, i) += gradE_strain(e, k) * grad_strain(e * stencil_dof + i, k);
            }
        }
    }

    // === 3. Compute Hessian ===
    Eigen::MatrixXd hess_energy = Eigen::MatrixXd::Zero(n_elements * stencil_dof, stencil_dof);

    for (int e = 0; e < n_elements; ++e) {
        const Eigen::MatrixXd& H = hess_strain[e];  // shape: (stencil_dof^2, n_K)

        for (int i = 0; i < stencil_dof; ++i) {
            for (int j = 0; j < stencil_dof; ++j) {
                double sum = 0.0;
                for (int k = 0; k < _n_K; ++k) {
                    // Term 1: gradE_strain * hess_strain
                    sum += gradE_strain(e, k) * H(i * stencil_dof + j, k);
                    // Term 2: K * outer product of grad_strain
                    sum += _K(e, k) * grad_strain(e * stencil_dof + i, k) * grad_strain(e * stencil_dof + j, k);
                }
                hess_energy(e * stencil_dof + i, j) = sum;
            }
        }
    }

    // === 4. Assemble into global vectors/matrices ===
    Eigen::VectorXd Fs = Eigen::VectorXd::Zero(n_dof);
    Eigen::MatrixXd Js = Eigen::MatrixXd::Zero(n_dof, n_dof);

    for (int e = 0; e < n_elements; ++e) {
        const auto& dof_indices = _ind[e];

        for (int i = 0; i < stencil_dof; ++i) {
            Fs[dof_indices[i]] -= grad_energy(e, i);
            for (int j = 0; j < stencil_dof; ++j) {
                Js(dof_indices[i], dof_indices[j]) -= hess_energy(e * stencil_dof + i, j);
            }
        }
    }

    return {Fs, Js};
}


std::pair<Eigen::VectorXd, Eigen::SparseMatrix<double>>
ElasticEnergy::grad_hess_energy_linear_elastic_sparse(const RobotState& state) const {
    const int n_elements = static_cast<int>(_ind.size());
    const int n_dof = state.q0.size();
    const int stencil_dof = static_cast<int>(_ind[0].size());

    // === 1. Get strain, grad, hess ===
    Eigen::MatrixXd strain = this->get_strain(state);
    auto [grad_strain, hess_strain] = this->grad_hess_strain(state);
    Eigen::MatrixXd del_strain = strain - _nat_strain;

    Eigen::MatrixXd gradE_strain = _K.array() * del_strain.array(); // (n_elements x n_K)

    // === 2. Compute gradient ===
    Eigen::MatrixXd grad_energy = Eigen::MatrixXd::Zero(n_elements, stencil_dof);
    for (int e = 0; e < n_elements; ++e) {
        for (int i = 0; i < stencil_dof; ++i) {
            for (int k = 0; k < _n_K; ++k) {
                grad_energy(e, i) += gradE_strain(e, k) * grad_strain(e * stencil_dof + i, k);
            }
        }
    }

    // === 3. Compute Hessian (flat vector for sparse triplets) ===
    std::vector<Eigen::Triplet<double>> triplets;

    for (int e = 0; e < n_elements; ++e) {
        const auto& dof_indices = _ind[e];
        const Eigen::MatrixXd& H = hess_strain[e];  // shape: (stencil_dof², n_K)

        for (int i = 0; i < stencil_dof; ++i) {
            for (int j = 0; j < stencil_dof; ++j) {
                double value = 0.0;
                for (int k = 0; k < _n_K; ++k) {
                    value += gradE_strain(e, k) * H(i * stencil_dof + j, k);
                    value += _K(e, k) * grad_strain(e * stencil_dof + i, k) * grad_strain(e * stencil_dof + j, k);
                }
                triplets.emplace_back(dof_indices[i], dof_indices[j], -value);
            }
        }
    }

    // === 4. Assemble global gradient and Hessian ===
    Eigen::VectorXd Fs = Eigen::VectorXd::Zero(n_dof);
    for (int e = 0; e < n_elements; ++e) {
        const auto& dof_indices = _ind[e];
        for (int i = 0; i < stencil_dof; ++i) {
            Fs[dof_indices[i]] -= grad_energy(e, i);
        }
    }

    Eigen::SparseMatrix<double> Js(n_dof, n_dof);
    Js.setFromTriplets(triplets.begin(), triplets.end());

    return {Fs, Js};
}

std::vector<std::vector<Eigen::Vector3d>> ElasticEnergy::get_node_positions(const Eigen::VectorXd& q) const {
    int N = _n_elems;
    int M = _n_nodes;

    std::cout << "N: " << N << ", M: " << M << std::endl;

    std::vector<std::vector<Eigen::Vector3d>> positions(M, std::vector<Eigen::Vector3d>(N));
    int node_id = 0;
    for (int i = 0; i < N; ++i) {  // for each spring
        for (int m = 0; m < M; ++m) {
            node_id = _nodes_ind[i][m];     
            std::cout << "i: " << i << ", m: " << m << ", node_id: " << node_id << std::endl;
            positions[m][i] = q.segment<3>(3 * node_id); 
            std::cout << "Spring " << i << ": node_id = " << node_id << " pos = " << positions[m][i].transpose() << std::endl;
        }
    }
    std::cout << "get node positions done" << std::endl;
    return positions;

    // std::vector<std::vector<Eigen::Vector3d>> positions(M, std::vector<Eigen::Vector3d>(N));
    // for (int i = 0; i < N; ++i) {         // Loop over springs/elements
    //     for (int m = 0; m < M; ++m) {     // Loop over nodes per spring
    //         int idx = i * M + m;
    //         int dof_idx = _node_dof_ind[idx];
    //         positions[m][i] = q.segment<3>(3 * dof_idx);  // 3D vector per node
    //     }
    // }
    // return positions;
}