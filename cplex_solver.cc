#include "cplex_solver.h"
#include "util.h"

std::string GetVariableName(const std::string prefix, int dimension,
                            const int* indices) {
  std::string ret = prefix;
  for (int i = 0; i < dimension; ++i) {
    ret += "[" + std::to_string(indices[i]) + "]";
  }
  return ret;
}

MuViNESolver::MuViNESolver(IPGraph* ip, OTNGraph* otn, DWDMGraph* dwdm,
                           IPGraph* vn, std::vector<std::vector<int>>* lc,
                           OverlayMapping<otn_edge_map_t>* ip_otn,
                           OverlayMapping<dwdm_edge_map_t>* otn_dwdm) {
  model_ = IloModel(env_);
  cplex_ = IloCplex(env_);
  constraints_ = IloConstraintArray(env_);
  objective_ = IloIntExpr(env_);
  ip_ = ip;
  otn_ = otn;
  dwdm_ = dwdm;
  vn_ = vn;
  location_constraints_ = lc;
  ip_otn_ = ip_otn;
  otn_dwdm_ = otn_dwdm;
  max_k_ = otn_->module_capacities()->size();
  w_ = dwdm_->num_wavelengths();
  c_ = dwdm_->wavelength_capacity();
  
  // Initialize bandwidth matrix.
  DEBUG("Initializing bw matrix.\n");
  b_uvi_.resize(ip_->node_count());
  cost_uvi_.resize(ip_->node_count());
  for (int u = 0; u < ip_->node_count(); ++u) {
    b_uvi_[u].resize(ip_->node_count());
    cost_uvi_[u].resize(ip_->node_count());
    auto& u_neighbors = ip_->adj_list()->at(u);
    for (auto end_point : u_neighbors) {
      int v = end_point.node_id();
      int order = end_point.order();
      int bandwidth = end_point.bandwidth();
      int cost = end_point.cost();
      if (order >= b_uvi_[u][v].size()) b_uvi_[u][v].resize(order + 1);
      if (order >= cost_uvi_[u][v].size()) cost_uvi_[u][v].resize(order + 1);
      b_uvi_[u][v][order] = bandwidth;
      cost_uvi_[u][v][order] = cost;
    }
    for (int v = 0; v < ip_->node_count(); ++v) {
      if (u == v) continue;
      long cap = std::min(ip_->GetPortCapacity(u), ip_->GetPortCapacity(v));
      b_uvi_[u][v].resize(ip_->GetPortCount(u), cap);
      cost_uvi_[u][v].resize(ip_->GetPortCount(u), 1);
    }
  }

  // m_pq_k_ -> Number of modules of type k installed on link (p, q).
  DEBUG("Initializing m_pq_k_\n");
  m_pq_k_.resize(otn_->node_count());
  for (int p = 0; p < otn_->node_count(); ++p) {
    m_pq_k_[p].resize(otn_->node_count());
    for (int q = 0; q < otn_->node_count(); ++q) {
      m_pq_k_[p][q].resize(max_k_, 0);  
      for (int k = 0; k < max_k_; ++k) {
        m_pq_k_[p][q][k] = otn_->GetNumModulesOnEdge(p, q, k);
      }
    }
  }

  // Initialize the decision variables.
  x_mn_uvi_ = IloIntVar5dArray(env_, vn_->node_count());
  y_mu_ = IloIntVar2dArray(env_, vn_->node_count());
  z_uvi_pqkj_ = IloIntVar7dArray(env_, ip_->node_count());
  gamma_uvi_ = IloIntVar3dArray(env_, ip_->node_count());
  zeta_pq_kj_ = IloIntVar4dArray(env_, otn_->node_count());
  phi_pqkj_l_ = IloIntVar5dArray(env_, otn_->node_count());
  psi_pqkj_abl_ = IloIntVar7dArray(env_, otn_->node_count());

  // Initialize x and y.
  DEBUG("Initializing x and y.\n");
  for (int m = 0; m < vn_->node_count(); ++m) {
    x_mn_uvi_[m] = IloIntVar4dArray(env_, vn_->node_count());
    y_mu_[m] = IloIntVarArray(env_, ip_->node_count(), 0, 1);
    for (int u = 0; u < vn_->node_count(); ++u) {
      int y_indices[] = {m, u};
      auto var_name = GetVariableName("y", 2, y_indices);
      y_mu_[m][u] = IloIntVar(env_, 0, 1, var_name.c_str());
    }
    for (int n = 0; n < vn_->node_count(); ++n) {
      x_mn_uvi_[m][n] = IloIntVar3dArray(env_, ip_->node_count());
      for (int u = 0; u < ip_->node_count(); ++u) {
        x_mn_uvi_[m][n][u] = IloIntVar2dArray(env_, ip_->node_count());
        for (int v = 0; v < ip_->node_count(); ++v) {
          x_mn_uvi_[m][n][u][v] = IloIntVarArray(env_, ip_->GetPortCount(u));
          for (int order = 0; order < ip_->GetPortCount(u); ++order) {
            int x_indices[] = {m, n, u, v, order};
            auto var_name = GetVariableName("x", 5, x_indices);
            x_mn_uvi_[m][n][u][v][order] =
                IloIntVar(env_, 0, 1, var_name.c_str());
          }
        }
      }
    }
  }

  // Initialize z and gamma.
  DEBUG("Initializing z.\n");
  for (int u = 0; u < ip_->node_count(); ++u) {
    z_uvi_pqkj_[u] = IloIntVar6dArray(env_, ip_->node_count());
    gamma_uvi_[u] = IloIntVar2dArray(env_, ip_->node_count());
    for (int v = 0; v < ip_->node_count(); ++v) {
      z_uvi_pqkj_[u][v] = IloIntVar5dArray(env_, ip_->GetPortCount(u));
      gamma_uvi_[u][v] = IloIntVarArray(env_, ip_->GetPortCount(u) + 1, 0, 1);
      for (int order = 0; order < ip_->GetPortCount(u); ++order) {
        z_uvi_pqkj_[u][v][order] = IloIntVar4dArray(env_, otn_->node_count());
        int gamma_indices[] = {u, v, order};
        auto vname = GetVariableName("gamma", 3, gamma_indices);
        gamma_uvi_[u][v][order] = IloIntVar(env_, 0, 1, vname.c_str());
        for (int p = 0; p < otn_->node_count(); ++p) {
          z_uvi_pqkj_[u][v][order][p] =
              IloIntVar3dArray(env_, otn_->node_count());
          for (int q = 0; q < otn_->node_count(); ++q) {
            z_uvi_pqkj_[u][v][order][p][q] = IloIntVar2dArray(env_, max_k_);
            for (int k = 0; k < max_k_; ++k) {
              z_uvi_pqkj_[u][v][order][p][q][k] =
                  IloIntVarArray(env_, m_pq_k_[p][q][k], 0, 1);
              for (int j = 0; j < m_pq_k_[p][q][k]; ++j) {
                int z_indices[] = {u, v, order, p, q, k, j};
                auto var_name = GetVariableName("z", 7, z_indices);
                z_uvi_pqkj_[u][v][order][p][q][k][j] =
                    IloIntVar(env_, 0, 1, var_name.c_str());
              }
            }
          }
        }
      }
    }
  }

  // Initialize zeta, phi, and psi.
  DEBUG("Initializing zeta, phi, and psi.\n");
  for (int p = 0; p < otn_->node_count(); ++p) {
    zeta_pq_kj_[p] = IloIntVar3dArray(env_, otn_->node_count());
    phi_pqkj_l_[p] = IloIntVar4dArray(env_, otn_->node_count());
    psi_pqkj_abl_[p] = IloIntVar6dArray(env_, otn_->node_count());
    for (int q = 0; q < otn_->node_count(); ++q) {
      zeta_pq_kj_[p][q] = IloIntVar2dArray(env_, max_k_);
      phi_pqkj_l_[p][q] = IloIntVar3dArray(env_, max_k_);
      psi_pqkj_abl_[p][q] = IloIntVar5dArray(env_, max_k_);
      for (int k = 0; k < max_k_; ++k) {
        zeta_pq_kj_[p][q][k] = IloIntVarArray(env_, m_pq_k_[p][q][k], 0, 1);
        phi_pqkj_l_[p][q][k] = IloIntVar2dArray(env_, m_pq_k_[p][q][k]);
        psi_pqkj_abl_[p][q][k] = IloIntVar4dArray(env_, m_pq_k_[p][q][k]);
        for (int j = 0; j < m_pq_k_[p][q][k]; ++j) {
          int zeta_indices[] = {p, q, k, j};
          auto var_name = GetVariableName("zeta", 4, zeta_indices);
          zeta_pq_kj_[p][q][k][j] = IloIntVar(env_, 0, 1, var_name.c_str());
          phi_pqkj_l_[p][q][k][j] = IloIntVarArray(env_, w_, 0, 1);
          psi_pqkj_abl_[p][q][k][j] =
              IloIntVar3dArray(env_, dwdm_->node_count());
          for (int l = 0; l < w_; ++l) {
            int phi_indices[] = {p, q, k, j, l};
            var_name = GetVariableName("phi", 5, phi_indices);
            phi_pqkj_l_[p][q][k][j][l] =
                IloIntVar(env_, 0, 1, var_name.c_str());
          }
          for (int a = 0; a < dwdm_->node_count(); ++a) {
            psi_pqkj_abl_[p][q][k][j][a] =
                IloIntVar2dArray(env_, dwdm_->node_count());
            for (int b = 0; b < dwdm_->node_count(); ++b) {
              psi_pqkj_abl_[p][q][k][j][a][b] = IloIntVarArray(env_, w_, 0, 1);
              for (int l = 0; l < w_; ++l) {
                int psi_indices[] = {p, q, k, j, a, b, l};
                var_name = GetVariableName("psi", 7, psi_indices);
                psi_pqkj_abl_[p][q][k][j][a][b][l] = IloIntVar(
                    env_, 0, 1, var_name.c_str());
              }
            }
          }
        }
      }
    }
  }

  // Initialize input variables.
  DEBUG("Initializing input variables.\n");
  l_mu_ = IloInt2dArray(env_, vn_->node_count());
  tau_up_ = IloInt2dArray(env_, ip_->node_count());
  omega_pq_kj_ = IloInt4dArray(env_, otn_->node_count());
  ip_link_uvi_ = IloInt3dArray(env_, ip_->node_count());

  // Initialize location constraint.
  DEBUG("Initializing location constraint.\n");
  for (int m = 0; m < vn_->node_count(); ++m) {
    l_mu_[m] = IloIntArray(env_, ip_->node_count(), 0, 1);
    for (int u = 0; u < ip_->node_count(); ++u) {
      l_mu_[m][u] = 0;
    }
    for (int j = 0; j < location_constraints_->at(m).size(); ++j) {
      int u = location_constraints_->at(m)[j];
      l_mu_[m][u] = 1;
    }
  }

  // Initialize IP link input variable.
  DEBUG("Initializing IP link input variable.\n");
  for (int u = 0; u < ip_->node_count(); ++u) {
    ip_link_uvi_[u] = IloInt2dArray(env_, ip_->node_count());
    for (int v = 0; v < ip_->node_count(); ++v) {
      ip_link_uvi_[u][v] = IloIntArray(env_, ip_->node_count(), 0, 1);
      for (int order = 0; order < ip_->GetPortCount(u); ++order) {
        ip_link_uvi_[u][v][order] = 0;
      }
    }
    auto& u_neighbors = ip_->adj_list()->at(u);
    for (auto end_point : u_neighbors) {
      int v = end_point.node_id();
      int order = end_point.order();
      ip_link_uvi_[u][v][order] = 1;
    }
  }

  // Initialize IP to OTN attachment variable.
  DEBUG("Initializing IP to OTN attachment variable.\n");
  for (int u = 0; u < ip_->node_count(); ++u) {
    tau_up_[u] = IloIntArray(env_, otn_->node_count(), 0, 1);
    for (int p = 0; p < otn_->node_count(); ++p) {
      tau_up_[u][p] = 0;
    }
    tau_up_[u][ip_otn_->node_map[u]] = 1;
  }

  // Initialize omega variable.
  DEBUG("Initializing omega variable.\n");
  for (int p = 0; p < otn_->node_count(); ++p) {
    omega_pq_kj_[p] = IloInt3dArray(env_, otn_->node_count());
    for (int q = 0; q < otn_->node_count(); ++q) {
      omega_pq_kj_[p][q] = IloInt2dArray(env_, max_k_);
      for (int k = 0; k < max_k_; ++k) {
        int n_mods = otn_->GetNumModulesOnEdge(p, q, k);
        omega_pq_kj_[p][q][k] = IloIntArray(env_, n_mods);
        for (int j = 0; j < n_mods; ++j) {
          if (otn_->module_capacities()->at(k) > 
              otn_->GetModuleResidualCapacity(p, q, k, j))
            omega_pq_kj_[p][q][k][j] = 1;
          else
            omega_pq_kj_[p][q][k][j] = 0;
        }
      }
    }
  }
}

void MuViNESolver::BuildModel() {
  // Constraint (1), (2)
  for (int m = 0; m < vn_->node_count(); ++m) {
    auto& m_neighbors = vn_->adj_list()->at(m);
    IloIntExpr sum(env_);
    for (int u = 0; u < ip_->node_count(); ++u) {
      sum += y_mu_[m][u];
      constraints_.add(y_mu_[m][u] <= l_mu_[m][u]);
    }
    constraints_.add(sum == 1);
  }

  // Constraint (3)
  for (int u = 0; u < ip_->node_count(); ++u) {
    IloIntExpr sum(env_);
    for (int m = 0; m < vn_->node_count(); ++m) {
      sum += y_mu_[m][u];
    }
    constraints_.add(sum <= 1);
  }

  // Constraint (4), (5), (6)
  for (int m = 0; m < vn_->node_count(); ++m) {
    auto& m_neighbors = vn_->adj_list()->at(m);
    for (auto vend_point : m_neighbors) {
      int n = vend_point.node_id();
      IloIntExpr sum_c4(env_);
      for (int u = 0; u < ip_->node_count(); ++u) {
        for (int v = 0; v < ip_->node_count(); ++v) {
          if (u == v) continue;
          int p_uv = std::min(ip_->GetPortCount(u), ip_->GetPortCount(v));
          for (int order = 0; order < p_uv; ++order) {
            constraints_.add(x_mn_uvi_[m][n][u][v][order] <=
                             gamma_uvi_[u][v][order] + gamma_uvi_[v][u][order] +
                                 ip_link_uvi_[u][v][order]);
          }
          IloIntExpr sum_c5(env_);
          for (int order = 0; order < ip_->GetPortCount(u); ++order) {
            sum_c4 += x_mn_uvi_[m][n][u][v][order];
            sum_c5 += x_mn_uvi_[m][n][u][v][order];
          }
          constraints_.add(sum_c4 >= 1);
          constraints_.add(sum_c5 <= 1);
        }
      }
    }
  }

  // Constraint (7)
  for (int u = 0; u < ip_->node_count(); ++u) {
    for (int v = 0; v < ip_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_->GetPortCount(u); ++order) {
        IloIntExpr sum(env_);
        for (int m = 0; m < vn_->node_count(); ++m) {
          auto& m_neighbors = vn_->adj_list()->at(m);
          for (auto vend_point : m_neighbors) {
            int n = vend_point.node_id();
            long b_mn = vend_point.bandwidth();
            sum += (x_mn_uvi_[m][n][u][v][order] * b_mn);
          }
        }
        constraints_.add(sum <= b_uvi_[u][v][order]);
      }
    }
  }

  // Constraint (8)
  for (int m = 0; m < vn_->node_count(); ++m) {
    auto& m_neighbors = vn_->adj_list()->at(m);
    for (auto vend_point : m_neighbors) {
      int n = vend_point.node_id();
      for (int u = 0; u < ip_->node_count(); ++u) {
        IloIntExpr sum(env_);
        for (int v = 0; v < ip_->node_count(); ++v) {
          if (u == v) continue;
          int p_uv = std::min(ip_->GetPortCount(u), ip_->GetPortCount(v));
          for (int order = 0; order < p_uv; ++order) {
            sum +=
                (x_mn_uvi_[m][n][u][v][order] - x_mn_uvi_[m][n][v][u][order]);
          }
        }
        constraints_.add(sum <= y_mu_[m][u] - y_mu_[n][u]);
      }
    }
  }

  // Constraint (9)
  for (int u = 0; u < ip_->node_count(); ++u) {
    IloIntExpr sum(env_);
    for (int v = 0; v < ip_->node_count(); ++v) {
      if (u == v) continue;
      int p_uv = std::min(ip_->GetPortCount(u), ip_->GetPortCount(v));
      for (int order = 0; order < p_uv; ++order) {
        sum += (gamma_uvi_[u][v][order] + gamma_uvi_[v][u][order] +
                ip_link_uvi_[u][v][order]);
      }
    }
    constraints_.add(sum <= ip_->GetPortCount(u));
  }

  // Constraint (10)
  for (int u = 0; u < ip_->node_count(); ++u) {
    for (int v = 0; v < ip_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_->GetPortCount(u); ++order) {
        constraints_.add(gamma_uvi_[u][v][order] + ip_link_uvi_[u][v][order] <=
                         1);
      }
    }
  }

  // Constraint (11), (12)
  for (int u = 0; u < ip_->node_count(); ++u) {
    for (int v = 0; v < ip_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_->GetPortCount(u); ++order) {
        for (int p = 0; p < otn_->node_count(); ++p) {
          auto& p_neighbors = otn_->adj_list()->at(p);
          for (auto end_point : p_neighbors) {
            int q = end_point.node_id();
            for (int k = 0; k < max_k_; ++k) {
              for (int j = 0; j < m_pq_k_[p][q][k]; ++j) {
                // (11)
                constraints_.add(z_uvi_pqkj_[u][v][order][p][q][k][j] <=
                                 gamma_uvi_[u][v][order]);
                // (12)
                constraints_.add(z_uvi_pqkj_[u][v][order][p][q][k][j] <=
                                 zeta_pq_kj_[p][q][k][j] + 
                                 omega_pq_kj_[p][q][k][j]);
              }
            }
          }
        }
      }
    }
  }

  // Constraint (13), (14), (15)
  for (int p = 0; p < otn_->node_count(); ++p) {
    auto& p_neighbors = otn_->adj_list()->at(p);
    for (auto end_point : p_neighbors) {
      int q = end_point.node_id();
      for (int k = 0; k < max_k_; ++k) {
        for (int j = 0; j < m_pq_k_[p][q][k]; ++j) {
          // (15)
          constraints_.add(zeta_pq_kj_[p][q][k][j] + omega_pq_kj_[p][q][k][j] <=
                           1);
          IloIntExpr sum(env_);
          for (int u = 0; u < ip_->node_count(); ++u) {
            for (int v = 0; v < ip_->node_count(); ++v) {
              if (u == v) continue;
              for (int order = 0; order < ip_->GetPortCount(u); ++order) {
                sum +=
                    z_uvi_pqkj_[u][v][order][p][q][k][j] * b_uvi_[u][v][order];
                // (14)
                constraints_.add(zeta_pq_kj_[p][q][k][j] <=
                                z_uvi_pqkj_[u][v][order][p][q][k][j]);
              }
            }
          }
          // (13)
          long module_res_cap = otn_->GetModuleResidualCapacity(p, q, k, j);
          constraints_.add(sum <= module_res_cap);
        }
      }
    }
  }

  // Constraint (16)
  for (int u = 0; u < ip_->node_count(); ++u) {
    for (int v = 0; v < ip_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_->GetPortCount(u); ++order) {
        for (int p = 0; p < otn_->node_count(); ++p) {
          auto& p_neighbors = otn_->adj_list()->at(p);
          IloIntExpr sum(env_);
          for (auto end_point : p_neighbors) {
            int q = end_point.node_id();
            for (int k = 0; k < max_k_; ++k) {
              for (int j = 0; j < m_pq_k_[p][q][k]; ++j) {
                sum += (z_uvi_pqkj_[u][v][order][p][q][k][j] -
                        z_uvi_pqkj_[u][v][order][q][p][k][j]);
              }
            }
          }
          if (tau_up_[u][p] == 1) {
            constraints_.add(sum == gamma_uvi_[u][v][order]);
          } else if (tau_up_[v][p] == 1) {
            constraints_.add(sum == -gamma_uvi_[u][v][order]);
          } else {
            constraints_.add(sum == 0);
          }
        }
      }
    }
  }

  // Constraint (17), (18)
  for (int p = 0; p < otn_->node_count(); ++p) {
    auto& p_neighbors = otn_->adj_list()->at(p);
    for (auto end_point : p_neighbors) {
      int q = end_point.node_id();
      for (int k = 0; k < max_k_; ++k) {
        for (int j = 0; j < m_pq_k_[p][q][k]; ++j) {
          IloIntExpr sum(env_);
          for (int l = 0; l < w_; ++l) {
            sum += phi_pqkj_l_[p][q][k][j][l];
          }
          // (17)
          constraints_.add(sum == zeta_pq_kj_[p][q][k][j]);
          for (int a = 0; a < dwdm_->node_count(); ++a) {
            auto& a_neighbors = dwdm_->adj_list()->at(a);
            for (auto aend_point : a_neighbors) {
              int b = aend_point.node_id();
              for (int l = 0; l < aend_point.is_lambda_free().size(); ++l) {
                // (18)
                if (!aend_point.is_lambda_free()[l]) continue;
                constraints_.add(psi_pqkj_abl_[p][q][k][j][a][b][l] <=
                                 phi_pqkj_l_[p][q][k][j][l]);
              }
            }
          }
        }
      }
    }
  }

  // Constraint (19)
  for (int p = 0; p < otn_->node_count(); ++p) {
    auto& p_neighbors = otn_->adj_list()->at(p);
    for (auto end_point : p_neighbors) {
      int q = end_point.node_id();
      for (int a = 0; a < dwdm_->node_count(); ++a) {
        auto& a_neighbors = dwdm_->adj_list()->at(a);
        for (auto aend_point : a_neighbors) {
          int b = aend_point.node_id();
          for (int l = 0; l < aend_point.is_lambda_free().size(); ++l) {
            if (!aend_point.is_lambda_free()[l]) continue;
            IloIntExpr sum(env_);
            for (int k = 0; k < max_k_; ++k) {
              for (int j = 0; j < m_pq_k_[p][q][k]; ++j) {
                sum += psi_pqkj_abl_[p][q][k][j][a][b][l] * 
                  otn_->module_capacities()->at(k);
              }
            }
            constraints_.add(sum <= c_);
          }
        }
      }
    }
  }

  // Constraint (20)
  for (int p = 0; p < otn_->node_count(); ++p) {
    auto& p_neighbors = otn_->adj_list()->at(p);
    for (auto end_point : p_neighbors) {
      int q = end_point.node_id();
      for (int k = 0; k < max_k_; ++k) {
        for (int j = 0; j < m_pq_k_[p][q][k]; ++j) {
          auto& dwdm_path =
              otn_dwdm_->edge_map[otn_edge_t(p, q, k, j)];
          for (auto link : dwdm_path) {
            int a = link.first, b = link.second;
            auto wavelength_mask = dwdm_->GetWavelengthMask(a, b);
            for (int l = 0; l < wavelength_mask.size(); ++l) {
              if (!wavelength_mask[l]) continue;
              constraints_.add(psi_pqkj_abl_[p][q][k][j][a][b][l] ==
                              phi_pqkj_l_[p][q][k][j][l]);
            }
          }
        }
      }
    }
  }

  // Objective function.
  // Component 1: Cost of embedding VLinks.
  for (int m = 0; m < vn_->node_count(); ++m) {
    auto& m_neighbors = vn_->adj_list()->at(m);
    for (auto vend_point : m_neighbors) {
      int n = vend_point.node_id();
      long b_mn = vend_point.bandwidth();
      for (int u = 0; u < ip_->node_count(); ++u) {
        for (int v = 0; v < ip_->node_count(); ++v) {
          if (u == v) continue;
          for (int order = 0; order < ip_->GetPortCount(u); ++order) {
            objective_ += (x_mn_uvi_[m][n][u][v][order] * b_mn * cost_uvi_[u][v][order]);
          }
        }
      }
    }
  }

  // Component 2: Cost of creating new IP links.
  for (int u = 0; u < ip_->node_count(); ++u) {
    for (int v = 0; v < ip_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_->GetPortCount(u); ++order) {
        for (int p = 0; p < otn_->node_count(); ++p) {
          auto& p_neighbors = otn_->adj_list()->at(p);
          for (auto pend_point : p_neighbors) {
            int q = pend_point.node_id();
            for (int k = 0; k < max_k_; ++k) {
              for (int j = 0; j < m_pq_k_[p][q][k]; ++j) {
                objective_ += (z_uvi_pqkj_[u][v][order][p][q][k][j] * b_uvi_[u][v][order] *  otn_->module_cost()->at(k));
              }
            }
          }
        }
      }
    }
  }

  // Component 3: Cost of activating new modules, i.e., routing new wavelengths.
  for (int p = 0; p < otn_->node_count(); ++p) {
    auto& p_neighbors = otn_->adj_list()->at(p);
    for (auto pend_point : p_neighbors) {
      int q = pend_point.node_id();
      for (int k = 0; k < max_k_; ++k) {
        for (int j = 0; j < m_pq_k_[p][q][k]; ++j) {
          for (int a = 0; a < dwdm_->node_count(); ++a) {
            auto& a_neighbors = dwdm_->adj_list()->at(a);
            for (auto aend_point : a_neighbors) {
              int b = aend_point.node_id();
              for (int l = 0; l < aend_point.is_lambda_free().size(); ++l) {
                if (!aend_point.is_lambda_free()[l]) continue;
                objective_ += (psi_pqkj_abl_[p][q][k][j][a][b][l] * aend_point.cost());
              }
            }
          }
        }        
      }
    }
  }
  constraints_.add(objective_ > 0);
  model_.add(constraints_);
  model_.add(IloMinimize(env_, objective_));
}


bool MuViNESolver::Solve() {
  int n_threads = sysconf(_SC_NPROCESSORS_ONLN);
  if (n_threads < 64) n_threads = 64;
  cplex_.setParam(IloCplex::Threads, n_threads);
  cplex_.exportModel("muvine.lp");
  bool success = cplex_.solve();
  if (cplex_.getStatus() == IloAlgorithm::Infeasible) {
    IloConstraintArray infeasible(env_);
    IloNumArray preferences(env_);
    infeasible.add(constraints_);
    for (int i = 0; i < infeasible.getSize(); ++i) preferences.add(1.0);
    if (cplex_.refineConflict(infeasible, preferences)) {
      IloCplex::ConflictStatusArray conflict = cplex_.getConflict(infeasible);
      env_.getImpl()->useDetailedDisplay(IloTrue);
      std::cout << "Conflict: " << std::endl;
      for (IloInt j = 0; j < infeasible.getSize(); ++j) {
        if (conflict[j] == IloCplex::ConflictMember) {
          std::cout << "Proved: " << infeasible[j] << std::endl;
        } else if (conflict[j] == IloCplex::ConflictPossibleMember) {
          std::cout << "Possible: " << infeasible[j] << std::endl;
        }
      }
    }
  }
  return success;
}
