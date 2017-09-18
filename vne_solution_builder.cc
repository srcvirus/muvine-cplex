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
                fprintf(outfile,"%d,%d,%d,%d,%d\n", m, n, u, v, order);
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
  const IloIntVar7dArray &z_uvi_pqkj = vne_solver_ptr_->z_uvi_pqkj();
  const IloIntVar3dArray &gamma_uvi = vne_solver_ptr_->gamma_uvi();
  for (int u = 0; u < ip_topology_->node_count(); ++u) {
    for (int v = 0; v < ip_topology_->node_count(); ++v) {
      if (u == v) continue;
      for (int order = 0; order < ip_topology_->GetPortCount(u); ++order) {
        if (fabs(cplex.getValue(gamma_uvi[u][v][order]) - 0) < EPS)
          continue;
        printf("gamma_uvi[%d][%d][%d] = 1\n", u, v, order);
        for (int p = 0; p < otn_topology_->node_count(); ++p) {
          const auto &p_neighbors = otn_topology_->adj_list()->at(p);
          for (const auto end_point : p_neighbors) {
            int q = end_point.node_id();
            const int kNumModuleTypes = otn_topology_->module_capacities()->size();
            for (int k = 0; k < kNumModuleTypes; ++k) {
              const int kNumModulesOnLink = otn_topology_->GetNumModulesOnEdge(p, q, k);
              for (int j = 0; j < kNumModulesOnLink; ++j) {
                if (fabs(cplex.getValue(z_uvi_pqkj[u][v][order][p][q][k][j]) - 1) < EPS) {
                  printf("New IP link (%d, %d, %d) --> OTN link (%d, %d, %d, %d)\n", u, v,
                         order, p, q, k, j);
                  if (outfile) {
                    fprintf(outfile, "%d,%d,%d,%d,%d,%d,%d\n", u, v, order, p, q, k, j);
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

void VNESolutionBuilder::PrintNewOTNModules(const char *filename) {
  FILE *outfile = nullptr;
  if (filename) outfile = fopen(filename, "r");
  const IloCplex &cplex = vne_solver_ptr_->cplex();
  const IloIntVar4dArray &zeta_pq_kj = vne_solver_ptr_->zeta_pq_kj();
  const IloIntVar5dArray &phi_pqkj_l = vne_solver_ptr_->phi_pqkj_l();
  const IloIntVar7dArray &psi_pqkj_abl = vne_solver_ptr_->psi_pqkj_abl();
  const int kNumModuleTypes = otn_topology_->module_capacities()->size();
  const int kNumLambdas = dwdm_topology_->num_wavelengths();
  for (int p = 0; p < otn_topology_->node_count(); ++p) {
    for (int q = 0; q < otn_topology_->node_count(); ++q) {
      if (p == q) continue;
      for (int k = 0; k < kNumModuleTypes; ++k) {
        const int kNumModulesOnLink = otn_topology_->GetNumModulesOnEdge(p, q, k);
        for (int j = 0; j < kNumModulesOnLink; ++j) {
          if (fabs(cplex.getValue(zeta_pq_kj[p][q][k][j]) - 0) < EPS)
            continue;
          for (int l = 0; l < kNumLambdas; ++l) {
            if (fabs(cplex.getValue(phi_pqkj_l[p][q][k][j][l]) - 0) < EPS)
              continue;
            for (int a = 0; a < dwdm_topology_->node_count(); ++a) {
              for (int b = 0; b < dwdm_topology_->node_count(); ++b) {
                if (a == b) continue;
                if (fabs(cplex.getValue(psi_pqkj_abl[p][q][k][j][a][b][l]) - 1) < EPS) {
                  printf("New OTN Module %d of type %d on (%d, %d) --> DWDM Link (%d, %d) with Lambda = %d\n", j, k, p, q, a, b, l);
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
