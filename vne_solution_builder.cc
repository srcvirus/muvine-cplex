#include "vne_solution_builder.h"
#include "util.h"

#include <math.h>
#include <algorithm>

void VNESolutionBuilder::PrintVLinkMapping(const char *filename) {
  FILE *outfile = nullptr;
  if (filename) outfile = fopen(filename, "w");
  const IloCplex &cplex = vne_solver_ptr_->cplex();
  const IloIntVar5dArray &x_mn_uvi = vne_solver_ptr_->x_mn_uvi();
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    const auto &m_neighbors = vn_topology_->adj_list()->at(m);
    for (const auto vend_point : m_neighbors) {
      int n = vend_point.node_id();
      if (m < n) continue;
      for (int u = 0; u < ip_topology_->node_count(); ++u) {
        for (int v = 0; v < ip_topology_->node_count(); ++v) {
          if (u == v) continue;
          for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
            if (fabs(cplex.getValue(x_mn_uvi[m][n][u][v][order]) - 1) < EPS) {
              printf("Virtual link (%d, %d) --> IP link (%d, %d, %d)\n", m, n,
                     u, v, order);
              if (outfile) {
                fprintf(outfile, "%d,%d,%d,%d,%d\n", m, n, u, v, order);
              }
            }
          }
        }
      }
    }
  }
  if (outfile) fclose(outfile);
}

void VNESolutionBuilder::PrintVNodeMapping(const char *filename) {
  FILE *outfile = nullptr;
  if (filename) outfile = fopen(filename, "w");
  const IloCplex &cplex = vne_solver_ptr_->cplex();
  const IloIntVar2dArray &y_m_u = vne_solver_ptr_->y_m_u();
  for (int m = 0; m < vn_topology_->node_count(); ++m) {
    for (int u = 0; u < ip_topology_->node_count(); ++u) {
      if (fabs(cplex.getValue(y_m_u[m][u]) - 1) < EPS) {
        printf("Virtual node %d --> IP node %d\n", m, u);
        if (outfile) {
          fprintf(outfile, "%d,%d\n", m, u);
        }
      }
    }
  }
  if (outfile) fclose(outfile);
}

void VNESolutionBuilder::PrintNewIPLinks(const char *filename) {
  FILE *outfile = nullptr;
  if (filename) outfile = fopen(filename, "w");
  const IloCplex &cplex = vne_solver_ptr_->cplex();
  const IloIntVar8dArray &z_uvi_pqkjl = vne_solver_ptr_->z_uvi_pqkjl();
  const IloIntVar3dArray &gamma_uvi = vne_solver_ptr_->gamma_uvi();
  for (int u = 0; u < ip_topology_->node_count(); ++u) {
    for (int v = 0; v < ip_topology_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
        if (fabs(cplex.getValue(gamma_uvi[u][v][order]) - 0) < EPS) continue;
        printf("gamma_uvi[%d][%d][%d] = 1\n", u, v, order);
        for (int p = 0; p < otn_topology_->node_count(); ++p) {
          for (int q = 0; q < otn_topology_->node_count(); ++q) {
            if (p == q) continue;
            const int kNumModuleTypes = otn_topology_->interface_info()->size();
            for (int k = 0; k < kNumModuleTypes; ++k) {
              for (int j = 0; j < otn_topology_->interfaces_installed()->at(p)[k]; ++j) {
                for (int l = 0; l < otn_topology_->interfaces_installed()->at(q)[k]; ++l) {
                  if (fabs(cplex.getValue(z_uvi_pqkjl[u][v][order][p][q][k][j][l]) -
                           1) < EPS) {
                    printf(
                        "New IP link (%d, %d, %d) --> OTN link (%d, %d, %d, "
                        "%d, %d)\n",
                        u, v, order, p, q, k, j, l);
                    if (outfile) {
                      fprintf(outfile, "%d,%d,%d,%d,%d,%d,%d,%d\n", u, v, order, p,
                              q, k, j, l);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (outfile) fclose(outfile);
}

void VNESolutionBuilder::PrintNewOTNLinks(const char *filename) {
  FILE *outfile = nullptr;
  if (filename) outfile = fopen(filename, "r");
  const IloCplex &cplex = vne_solver_ptr_->cplex();
  const IloIntVar5dArray &zeta_pq_kjl = vne_solver_ptr_->zeta_pq_kjl();
  const IloIntVar6dArray &phi_pqkjl_l = vne_solver_ptr_->phi_pqkjl_l();
  const IloIntVar8dArray &psi_pqkjl_abl = vne_solver_ptr_->psi_pqkjl_abl();
  const int kNumModuleTypes = otn_topology_->interface_info()->size();
  const int kNumLambdas = dwdm_topology_->num_wavelengths();
  for (int p = 0; p < otn_topology_->node_count(); ++p) {
    for (int q = 0; q < otn_topology_->node_count(); ++q) {
      if (p == q) continue;
      for (int k = 0; k < kNumModuleTypes; ++k) {
        for (int j = 0; j < otn_topology_->interfaces_installed()->at(p)[k]; ++j) {
          for (int l = 0; l < otn_topology_->interfaces_installed()->at(q)[k]; ++l) {
            if (fabs(cplex.getValue(zeta_pq_kjl[p][q][k][j][l]) - 0) < EPS) continue;
            DEBUG("zeta_pq_kj[%d][%d][%d][%d] = 1\n", p, q, k, j);
            for (int ll = 0; ll < kNumLambdas; ++ll) {
              if (fabs(cplex.getValue(phi_pqkjl_l[p][q][k][j][l][ll]) - 0) < EPS)
                continue;
              DEBUG("phi_pqkj_l[%d][%d][%d][%d][%d] = 1\n", p, q, k, j, l);
              for (int a = 0; a < dwdm_topology_->node_count(); ++a) {
                auto& a_neighbors = dwdm_topology_->adj_list()->at(a);
                for (auto& aend_point : a_neighbors) {
                  int b = aend_point.node_id();                  
                  if (fabs(cplex.getValue(psi_pqkjl_abl[p][q][k][j][l][a][b][ll]) - 1) <
                    EPS) {
                    printf(
                        "New OTN Module %d of type %d on (%d, %d) --> DWDM Link "
                        "(%d, %d) with Lambda = %d\n",
                        j, k, p, q, a, b, l);
                    if (outfile) {
                      printf("%d,%d,%d,%d,%d,%d,%d\n", p, q, k, j, a, b, l);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (outfile) fclose(outfile);
}

void VNESolutionBuilder::PrintSolutionStatus(const char *filename) {
  const IloCplex &cplex = vne_solver_ptr_->cplex();
  std::cout << "Solution status = " << cplex.getStatus() << std::endl;
  if (filename) {
    std::ofstream ofs(filename);
    ofs << cplex.getStatus() << std::endl;
    ofs.close();
  }
}

void VNESolutionBuilder::PrintCost(const char *filename) {
  FILE *outfile = nullptr;
  if (filename) outfile = fopen(filename, "w");
  const IloCplex &cplex = vne_solver_ptr_->cplex();
  printf("Cost = %lf\n", cplex.getObjValue());
  if (outfile) {
    fprintf(outfile, "%lf\n", cplex.getObjValue());
    fclose(outfile);
  }
}
