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

MuViNESolver::MuViNESolver(
    IPGraph* ip, OTNGraph* otn, DWDMGraph* dwdm,
    IPGraph* vn, std::vector<std::vector<int>>* lc,
    OverlayMapping<otn_edge_map_t>* ip_otn,
    OverlayMapping<dwdm_edge_map_t>* otn_dwdm) {
  model_ = IloModel(env_);
  cplex_ = IloCplex(model_);
  constraints_ = IloConstraintArray(env_);
  objective_ = IloIntExpr(env_);
  ip_ = ip;
  otn_ = otn;
  dwdm_ = dwdm;
  vn_ = vn;
  location_constraints_ = lc;
  ip_otn_ = ip_otn;
  otn_dwdm_ = otn_dwdm;
  max_k_ = otn_->interface_info()->size();
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

  // m_p_k_ -> Number of modules of type k installed on node p.
  DEBUG("Initializing m_p_k_\n");
  m_p_k_.resize(otn_->node_count());
  for (int p = 0; p < otn_->node_count(); ++p) {
    m_p_k_[p].resize(max_k_, 0);
    for (int k = 0; k < max_k_; ++k) {
      m_p_k_[p][k] = otn_->GetInterfaceCount(p, k);
    }
  }

  // q_pq_kjl_ -> bandwidth matrix for OTN layer
  q_pq_kjl_.resize(otn_->node_count());
  for (int p = 0; p < otn_->node_count(); ++p) {
    q_pq_kjl_[p].resize(otn_->node_count());
    for (int q = 0; q < otn_->node_count(); ++q) {
      q_pq_kjl_[p][q].resize(max_k_);
      for (int k = 0; k < max_k_; ++k) {
        q_pq_kjl_[p][q][k].resize(m_p_k_[p][k]);
        for (int j = 0; j < m_p_k_[p][k]; ++j) {
          q_pq_kjl_[p][q][k][j].resize(
              m_p_k_[q][k], otn_->interface_info()->at(k).capacity);
        }
      }
    }   
    auto& p_neighbors = otn_->adj_list()->at(p);
    for (auto end_point : p_neighbors) {
      int q = end_point.node_id();
      int k = end_point.intf_type();
      int j = end_point.src_intf();
      int l = end_point.dst_intf();
      long rbw = end_point.residual_bandwidth();
      q_pq_kjl_[p][q][k][j][l] = rbw;
    }
  }
  

  // Initialize the decision variables.
  x_mn_uvi_ = IloIntVar5dArray(env_, vn_->node_count());
  y_mu_ = IloIntVar2dArray(env_, vn_->node_count());
  z_uvi_pqkjl_ = IloIntVar8dArray(env_, ip_->node_count());
  gamma_uvi_ = IloIntVar3dArray(env_, ip_->node_count());
  zeta_pq_kjl_ = IloIntVar5dArray(env_, otn_->node_count());
  phi_pqkjl_l_ = IloIntVar6dArray(env_, otn_->node_count());
  psi_pqkjl_abl_ = IloIntVar8dArray(env_, otn_->node_count());

  // Initialize x and y.
  DEBUG("Initializing x and y.\n");
  for (int m = 0; m < vn_->node_count(); ++m) {
    x_mn_uvi_[m] = IloIntVar4dArray(env_, vn_->node_count());
    y_mu_[m] = IloIntVarArray(env_, ip_->node_count(), 0, 1);
    for (int u = 0; u < ip_->node_count(); ++u) {
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
    z_uvi_pqkjl_[u] = IloIntVar7dArray(env_, ip_->node_count());
    gamma_uvi_[u] = IloIntVar2dArray(env_, ip_->node_count());
    for (int v = 0; v < ip_->node_count(); ++v) {
      z_uvi_pqkjl_[u][v] = IloIntVar6dArray(env_, ip_->GetPortCount(u));
      gamma_uvi_[u][v] = IloIntVarArray(env_, ip_->GetPortCount(u) + 1, 0, 1);
      for (int order = 0; order < ip_->GetPortCount(u); ++order) {
        z_uvi_pqkjl_[u][v][order] = IloIntVar5dArray(env_, otn_->node_count());
        int gamma_indices[] = {u, v, order};
        auto vname = GetVariableName("gamma", 3, gamma_indices);
        gamma_uvi_[u][v][order] = IloIntVar(env_, 0, 1, vname.c_str());
        for (int p = 0; p < otn_->node_count(); ++p) {
          z_uvi_pqkjl_[u][v][order][p] =
              IloIntVar4dArray(env_, otn_->node_count());
          for (int q = 0; q < otn_->node_count(); ++q) {
            z_uvi_pqkjl_[u][v][order][p][q] = IloIntVar3dArray(env_, max_k_);
            for (int k = 0; k < max_k_; ++k) {
              z_uvi_pqkjl_[u][v][order][p][q][k] =
                  IloIntVar2dArray(env_, m_p_k_[p][k]);
              for (int j = 0; j < m_p_k_[p][k]; ++j) {
                z_uvi_pqkjl_[u][v][order][p][q][k][j] =
                  IloIntVarArray(env_, m_p_k_[q][k], 0, 1);
                for (int l = 0; l < m_p_k_[q][k]; ++l) {
                  int z_indices[] = {u, v, order, p, q, k, j, l};
                  auto var_name = GetVariableName("z", 8, z_indices);
                  z_uvi_pqkjl_[u][v][order][p][q][k][j][l] =
                      IloIntVar(env_, 0, 1, var_name.c_str());
                }
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
    zeta_pq_kjl_[p] = IloIntVar4dArray(env_, otn_->node_count());
    phi_pqkjl_l_[p] = IloIntVar5dArray(env_, otn_->node_count());
    psi_pqkjl_abl_[p] = IloIntVar7dArray(env_, otn_->node_count());
    for (int q = 0; q < otn_->node_count(); ++q) {
      zeta_pq_kjl_[p][q] = IloIntVar3dArray(env_, max_k_);
      phi_pqkjl_l_[p][q] = IloIntVar4dArray(env_, max_k_);
      psi_pqkjl_abl_[p][q] = IloIntVar6dArray(env_, max_k_);
      for (int k = 0; k < max_k_; ++k) {
        zeta_pq_kjl_[p][q][k] = IloIntVar2dArray(env_, m_p_k_[p][k]);
        phi_pqkjl_l_[p][q][k] = IloIntVar3dArray(env_, m_p_k_[p][k]);
        psi_pqkjl_abl_[p][q][k] = IloIntVar5dArray(env_, m_p_k_[p][k]);
        for (int j = 0; j < m_p_k_[p][k]; ++j) {
          zeta_pq_kjl_[p][q][k][j] = IloIntVarArray(env_, m_p_k_[q][k], 0, 1);
          phi_pqkjl_l_[p][q][k][j] = IloIntVar2dArray(env_, m_p_k_[q][k]);
          psi_pqkjl_abl_[p][q][k][j] = IloIntVar4dArray(env_, m_p_k_[q][k]);
          for (int l = 0; l < m_p_k_[q][k]; ++l) {
            int zeta_indices[] = {p, q, k, j, l};
            auto var_name = GetVariableName("zeta", 5, zeta_indices);
            zeta_pq_kjl_[p][q][k][j][l] = IloIntVar(env_, 0, 1, var_name.c_str());
            phi_pqkjl_l_[p][q][k][j][l] = IloIntVarArray(env_, w_, 0, 1);
            psi_pqkjl_abl_[p][q][k][j][l] =
                IloIntVar3dArray(env_, dwdm_->node_count());
            for (int ll = 0; ll < w_; ++ll) {
              int phi_indices[] = {p, q, k, j, l, ll};
              var_name = GetVariableName("phi", 6, phi_indices);
              phi_pqkjl_l_[p][q][k][j][l][ll] =
                  IloIntVar(env_, 0, 1, var_name.c_str());
            }
            for (int a = 0; a < dwdm_->node_count(); ++a) {
              psi_pqkjl_abl_[p][q][k][j][l][a] =
                  IloIntVar2dArray(env_, dwdm_->node_count());
              for (int b = 0; b < dwdm_->node_count(); ++b) {
                psi_pqkjl_abl_[p][q][k][j][l][a][b] = IloIntVarArray(env_, w_, 0, 1);
                for (int ll = 0; ll < w_; ++ll) {
                  int psi_indices[] = {p, q, k, j, l, a, b, ll};
                  var_name = GetVariableName("psi", 8, psi_indices);
                  psi_pqkjl_abl_[p][q][k][j][l][a][b][ll] =
                      IloIntVar(env_, 0, 1, var_name.c_str());
                }
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
  nu_pa_ = IloInt2dArray(env_, otn_->node_count());
  omega_pq_kjl_ = IloInt5dArray(env_, otn_->node_count());
  ip_link_uvi_ = IloInt3dArray(env_, ip_->node_count());

  // Initialize location constraint.
  DEBUG("Initializing location constraint.\n");
  for (int m = 0; m < vn_->node_count(); ++m) {
    DEBUG("IP node count = %d\n", ip_->node_count());
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

  // Initialize OTN to DWDM attachment variable.
  DEBUG("Initializing OTN to DWDM attachment variable.\n");
  for (int p = 0; p < otn_->node_count(); ++p) {
    nu_pa_[p] = IloIntArray(env_, dwdm_->node_count(), 0, 1);
    for (int a = 0; a < dwdm_->node_count(); ++a) {
      nu_pa_[p][a] = 0;
    }
    nu_pa_[p][otn_dwdm_->node_map[p]] = 1;
  }

  // Initialize omega variable.
  DEBUG("Initializing omega variable.\n");
  for (int p = 0; p < otn_->node_count(); ++p) {
    omega_pq_kjl_[p] = IloInt4dArray(env_, otn_->node_count());
    for (int q = 0; q < otn_->node_count(); ++q) {
      omega_pq_kjl_[p][q] = IloInt3dArray(env_, max_k_);
      for (int k = 0; k < max_k_; ++k) {
        omega_pq_kjl_[p][q][k] = IloInt2dArray(env_, m_p_k_[p][k]);
        for (int j = 0; j < m_p_k_[p][k]; ++j) {
          omega_pq_kjl_[p][q][k][j] = IloIntArray(env_, m_p_k_[q][k], 0, 1);
          for (int l = 0; l < m_p_k_[q][k]; ++l) {
            omega_pq_kjl_[p][q][k][j][l] = 0;
          }
        }
      }
    }
  }

  DEBUG("Populating omega.\n");
  for (auto it = ip_otn_->edge_map.begin(); it != ip_otn_->edge_map.end();
       ++it) {
    auto& otn_path = it->second;
    for (auto otn_edge : otn_path) {
      int p = otn_edge.first;
      int q = otn_edge.second;
      int k = otn_edge.interface_type;
      int j = otn_edge.src_idx;
      int l = otn_edge.dst_idx;
      omega_pq_kjl_[p][q][k][j][l] = 1;
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

      // Constraint (2).
      constraints_.add(y_mu_[m][u] <= l_mu_[m][u]);
    }

    // Constraint (1).
    constraints_.add(sum == 1);
  }
  printf("Constraints (1), (2) added\n");

  // Constraint (3)
  for (int u = 0; u < ip_->node_count(); ++u) {
    IloIntExpr sum(env_);
    for (int m = 0; m < vn_->node_count(); ++m) {
      sum += y_mu_[m][u];
    }
    constraints_.add(sum <= 1);
  }
  printf("Constraints (3) added\n");

  // Constraint (4), (5), (6)
  for (int m = 0; m < vn_->node_count(); ++m) {
    auto& m_neighbors = vn_->adj_list()->at(m);
    for (auto vend_point : m_neighbors) {
      int n = vend_point.node_id();
      if (m < n) continue;
      IloIntExpr sum_c5(env_);
      for (int u = 0; u < ip_->node_count(); ++u) {
        for (int v = 0; v < ip_->node_count(); ++v) {
          if (u == v) continue;
          int p_uv = std::min(ip_->GetPortCount(u), ip_->GetPortCount(v));
          for (int order = 0; order < p_uv; ++order) {
            // (4)
            constraints_.add(x_mn_uvi_[m][n][u][v][order] <=
                             gamma_uvi_[u][v][order] + gamma_uvi_[v][u][order] +
                                 ip_link_uvi_[u][v][order]);
          }
          IloIntExpr sum_c6(env_);
          for (int order = 0; order < ip_->GetPortCount(u); ++order) {
            sum_c5 += x_mn_uvi_[m][n][u][v][order];
            sum_c6 += x_mn_uvi_[m][n][u][v][order];
            for (int other_order = 0; other_order < ip_->GetPortCount(v);
                 ++other_order) {
              // Do not map (m, n) on both (u, v) and (v, u).
              constraints_.add(
                  IloIfThen(env_, x_mn_uvi_[m][n][u][v][order] == 1,
                            x_mn_uvi_[m][n][v][u][other_order] == 0));
              constraints_.add(
                  IloIfThen(env_, x_mn_uvi_[m][n][v][u][other_order] == 1,
                            x_mn_uvi_[m][n][u][v][order] == 0));
              // If (u, v, order) is activated then do not activate any of (v,
              // u, other_order).
              constraints_.add(
                  gamma_uvi_[u][v][order] + gamma_uvi_[v][u][other_order] <= 1);
            }
          }
          // (6)
          constraints_.add(sum_c6 <= 1);
        }
      }
      // (5)
      constraints_.add(sum_c5 >= 1);
    }
  }  
  printf("Constraints (4), (5), (6) added\n");

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
            if (m < n) continue;
            long b_mn = vend_point.bandwidth();
            sum += (x_mn_uvi_[m][n][u][v][order] * b_mn);
          }
        }
        constraints_.add(sum <= b_uvi_[u][v][order]);
      }
    }
  }
  printf("Constraints (7) added\n");

  // Constraint (8)
  for (int m = 0; m < vn_->node_count(); ++m) {
    auto& m_neighbors = vn_->adj_list()->at(m);
    for (auto vend_point : m_neighbors) {
      int n = vend_point.node_id();
      if (m < n) continue;
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
        constraints_.add((sum == (y_mu_[m][u] - y_mu_[n][u])));
      }
    }
  }
  printf("Constraints (8) added\n");

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
  printf("Constraints (9) added\n");

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
  printf("Constraints (10) added\n");

  // Constraint (11), (12)
  for (int u = 0; u < ip_->node_count(); ++u) {
    for (int v = 0; v < ip_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_->GetPortCount(u); ++order) {
        for (int p = 0; p < otn_->node_count(); ++p) {
          for (int q = 0; q < otn_->node_count(); ++q) {
            if (p == q) continue;
            for (int k = 0; k < max_k_; ++k) {
              for (int j = 0; j < m_p_k_[p][k]; ++j) {
                for (int l = 0; l < m_p_k_[q][k]; ++l) {
                  // (11)
                  constraints_.add(z_uvi_pqkjl_[u][v][order][p][q][k][j][l] <=
                                   gamma_uvi_[u][v][order]);
                  // (12)
                  constraints_.add(z_uvi_pqkjl_[u][v][order][p][q][k][j][l] <=
                                   zeta_pq_kjl_[p][q][k][j][l] +
                                       omega_pq_kjl_[p][q][k][j][l]);
                }
              }
            }
          }
        }
      }
    }
  }
  printf("Constraints (11), (12) added\n");

  // Constraint (13), (15)
  for (int p = 0; p < otn_->node_count(); ++p) {
    for (int q = 0; q < otn_->node_count(); ++q) {
      if (p == q) continue;
      for (int k = 0; k < max_k_; ++k) {
        for (int j = 0; j < m_p_k_[p][k]; ++j) {
          for (int l = 0; l < m_p_k_[q][k]; ++l) {
            // (15)
            constraints_.add(zeta_pq_kjl_[p][q][k][j][l] + 
                omega_pq_kjl_[p][q][k][j][l] <= 1);
            IloIntExpr sum(env_);
            for (int u = 0; u < ip_->node_count(); ++u) {
              for (int v = 0; v < ip_->node_count(); ++v) {
                if (u == v) continue;
                for (int order = 0; order < ip_->GetPortCount(u); ++order) {
                  sum += z_uvi_pqkjl_[u][v][order][p][q][k][j][l] * 
                      b_uvi_[u][v][order];
                }
              }
            }
            // (13)
            constraints_.add(sum <= q_pq_kjl_[p][q][k][j][l]);
          }
        }
      }
    }
  }
  printf("Constraints (13), (15) added\n");

  // Constraint (14).
  for (int p = 0; p < otn_->node_count(); ++p) {
    for (int k = 0; k < max_k_; ++k) {
      IloIntExpr sum(env_);
      for (int q = 0; q < otn_->node_count(); ++q) {
        if (p == q) continue;
        for (int j = 0; j < m_p_k_[p][k]; ++j) {
          for (int l = 0; l < m_p_k_[q][k]; ++l) {
            sum += (zeta_pq_kjl_[p][q][k][j][l] + 
                zeta_pq_kjl_[q][p][k][l][j] + omega_pq_kjl_[p][q][k][j][l]);
          }
        }
      }
      // (14)
      constraints_.add(sum <= m_p_k_[p][k]);
    }
  }
  printf("Constraints (14) added\n");

  // Constraint (16)
  for (int u = 0; u < ip_->node_count(); ++u) {
    for (int v = 0; v < ip_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_->GetPortCount(u); ++order) {
        for (int p = 0; p < otn_->node_count(); ++p) {
          IloIntExpr sum(env_);
          for (int q = 0; q < otn_->node_count(); ++q) {
            if (p == q) continue;
            for (int k = 0; k < max_k_; ++k) {
              for (int j = 0; j < m_p_k_[p][k]; ++j) {
                for (int l = 0; l < m_p_k_[q][k]; ++l) {
                  sum += (z_uvi_pqkjl_[u][v][order][p][q][k][j][l] -
                          z_uvi_pqkjl_[u][v][order][q][p][k][j][l]);
                  for (int other_k = 0; other_k < max_k_; ++other_k) {
                    for (int other_j = 0; other_j < m_p_k_[p][other_k];
                         ++other_j) {
                      for (int other_l = 0; other_l < m_p_k_[q][other_k]; ++other_l) {
                        constraints_.add(IloIfThen(
                            env_, z_uvi_pqkjl_[u][v][order][p][q][k][j][l] == 1,
                            z_uvi_pqkjl_[u][v][order][q][p][other_k][other_l][other_j] == 0));
                        constraints_.add(IloIfThen(
                            env_,
                            z_uvi_pqkjl_[u][v][order][q][p][other_k][other_l][other_j] == 1,
                            z_uvi_pqkjl_[u][v][order][p][q][k][j][l] == 0));
                        constraints_.add(zeta_pq_kjl_[p][q][k][j][l] +
                                             zeta_pq_kjl_[q][p][other_k][other_l][other_j] <=
                                         1);
                      }
                    }
                  }
                }
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
  printf("Constraints (16) added\n");

  // Constraint (17), (18)
  for (int p = 0; p < otn_->node_count(); ++p) {
    for (int q = 0; q < otn_->node_count(); ++q) {
      if (p == q) continue;
      for (int k = 0; k < max_k_; ++k) {
        for (int j = 0; j < m_p_k_[p][k]; ++j) {
          for (int l = 0; l < m_p_k_[q][k]; ++l) {
            IloIntExpr sum(env_);
            for (int ll = 0; ll < w_; ++ll) {
              sum += phi_pqkjl_l_[p][q][k][j][l][ll];
            }

            // (17)
            constraints_.add(sum == zeta_pq_kjl_[p][q][k][j][l]);

            for (int a = 0; a < dwdm_->node_count(); ++a) {
            auto& a_neighbors = dwdm_->adj_list()->at(a);
              for (auto aend_point : a_neighbors) {
                int b = aend_point.node_id();
                for (int ll = 0; ll < w_; ++ll) {
                  // (18)
                  constraints_.add(psi_pqkjl_abl_[p][q][k][j][l][a][b][ll] <=
                                 phi_pqkjl_l_[p][q][k][j][l][ll]);
                }
              }
            }
          }
        }
      }
    }
  }
  printf("Constraints (17), (18) added\n");

  // Constraint (20)
  for (int a = 0; a < dwdm_->node_count(); ++a) {
    auto& a_neighbors = dwdm_->adj_list()->at(a);
    for (auto aend_point : a_neighbors) {
      int b = aend_point.node_id();
      for (int ll = 0; ll < w_; ++ll) {
        IloIntExpr sum(env_);
        for (int p = 0; p < otn_->node_count(); ++p) {
          for (int q = 0; q < otn_->node_count(); ++q) {
            if (p == q) continue;
            for (int k = 0; k < max_k_; ++k) {
              for (int j = 0; j < m_p_k_[p][k]; ++j) {
                for (int l = 0; l < m_p_k_[q][k]; ++l) {
                  sum += psi_pqkjl_abl_[p][q][k][j][l][a][b][ll] *
                          otn_->interface_info()->at(k).capacity;
                }
              }
            }
          }
          constraints_.add(sum <= aend_point.lambda_residual_bandwidth()[ll]);
        }
      }
    }
  }
  printf("Constraints (20) added\n");

  // Constraint (21)
  for (int p = 0; p < otn_->node_count(); ++p) {
    for (int q = 0; q < otn_->node_count(); ++q) {
      if (p == q) continue;
      for (int k = 0; k < max_k_; ++k) {
        for (int j = 0; j < m_p_k_[p][k]; ++j) {
          for (int l = 0; l < m_p_k_[q][k]; ++l) {
            for (int a = 0; a < dwdm_->node_count(); ++a) {
              IloIntExpr sum(env_);
              auto a_neighbors = dwdm_->adj_list()->at(a);
              for (auto aend_point : a_neighbors) {
                int b = aend_point.node_id();
                for (int ll = 0; ll < w_; ++ll) {
                  sum += (psi_pqkjl_abl_[p][q][k][j][l][a][b][ll] - 
                      psi_pqkjl_abl_[p][q][k][j][l][b][a][ll]);
                }
              }
              if (nu_pa_[p][a] == 1) {
                constraints_.add(sum == zeta_pq_kjl_[p][q][k][j][l]);
              } else if (nu_pa_[q][a] == 1) {
                constraints_.add(sum == -zeta_pq_kjl_[p][q][k][j][l]);
              } else constraints_.add(sum == 0);
            }
          }
        }
      }
    }
  }
  printf("Constraints (21) added\n");

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
            objective_ +=
                (x_mn_uvi_[m][n][u][v][order] * b_mn * cost_uvi_[u][v][order]);
          }
        }
      }
    }
  }
  printf("Objective component 1 added\n");

  // Component 2: Cost of creating new IP links.
  for (int u = 0; u < ip_->node_count(); ++u) {
    for (int v = 0; v < ip_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_->GetPortCount(u); ++order) {
        for (int p = 0; p < otn_->node_count(); ++p) {
          for (int q = 0; q < otn_->node_count(); ++q) {
            if (p == q) continue;
            for (int k = 0; k < max_k_; ++k) {
              for (int j = 0; j < m_p_k_[p][k]; ++j) {
                for (int l = 0; l < m_p_k_[q][k]; ++l) {
                  objective_ +=
                      (z_uvi_pqkjl_[u][v][order][p][q][k][j][l] *
                       b_uvi_[u][v][order] * otn_->interface_info()->at(k).cost);
                }
              }
            }
          }
        }
      }
    }
  }
  printf("Objective component 2 added\n");

  // Component 3: Cost of activating new modules, i.e., routing new wavelengths.
  for (int p = 0; p < otn_->node_count(); ++p) {
    for (int q = 0; q < otn_->node_count(); ++q) {
      for (int k = 0; k < max_k_; ++k) {
        for (int j = 0; j < m_p_k_[p][k]; ++j) {
          for (int l = 0; l < m_p_k_[q][k]; ++l) {
            for (int a = 0; a < dwdm_->node_count(); ++a) {
              auto& a_neighbors = dwdm_->adj_list()->at(a);
              for (auto aend_point : a_neighbors) {
                int b = aend_point.node_id();
                for (int ll = 0; ll < w_; ++ll) {
                  objective_ +=
                      (psi_pqkjl_abl_[p][q][k][j][l][a][b][ll] * aend_point.cost());
                }
              }
            }
          }
        }
      }
    }
  }
  printf("Objective component 3 added\n");
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
