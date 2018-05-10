#ifndef CPLEX_SOLVER_H_
#define CPLEX_SOLVER_H_

#include <ilcplex/ilocplex.h>
#include "datastructure.h"

typedef IloArray<IloIntVarArray> IloIntVar2dArray;
typedef IloArray<IloIntVar2dArray> IloIntVar3dArray;
typedef IloArray<IloIntVar3dArray> IloIntVar4dArray;
typedef IloArray<IloIntVar4dArray> IloIntVar5dArray;
typedef IloArray<IloIntVar5dArray> IloIntVar6dArray;
typedef IloArray<IloIntVar6dArray> IloIntVar7dArray;
typedef IloArray<IloIntVar7dArray> IloIntVar8dArray;

typedef IloArray<IloIntArray> IloInt2dArray;
typedef IloArray<IloInt2dArray> IloInt3dArray;
typedef IloArray<IloInt3dArray> IloInt4dArray;
typedef IloArray<IloInt4dArray> IloInt5dArray;

using std::string;

class MuViNESolver {
 public:
  MuViNESolver() {}
  MuViNESolver(IPGraph* ip, OTNGraph* otn, DWDMGraph* dwdm, IPGraph* vn,
               std::vector<std::vector<int>>* lc,
               OverlayMapping<otn_edge_map_t>* ip_otn,
               OverlayMapping<dwdm_edge_map_t>* otn_dwdm);
  IloEnv& env() { return env_; }
  IloModel& model() { return model_; }
  IloCplex& cplex() { return cplex_; }
  IloConstraintArray& constraints() { return constraints_; }
  IloIntVar5dArray& x_mn_uvi() { return x_mn_uvi_; }
  IloIntVar2dArray& y_m_u() { return y_mu_; }
  IloIntVar3dArray& gamma_uvi() { return gamma_uvi_; }
  IloIntVar8dArray& z_uvi_pqkjl() { return z_uvi_pqkjl_; }
  IloIntVar5dArray& zeta_pq_kjl() { return zeta_pq_kjl_; }
  IloIntVar6dArray& phi_pqkjl_l() { return phi_pqkjl_l_; }
  IloIntVar8dArray& psi_pqkjl_abl() { return psi_pqkjl_abl_; }
  void BuildModel();
  bool Solve();

 private:
  // CPLEX related variables.
  IloEnv env_;
  IloModel model_;
  IloCplex cplex_;
  IloConstraintArray constraints_;
  IloNumArray preferences_;

  // Input data structures.
  IPGraph* ip_;
  OTNGraph* otn_;
  DWDMGraph* dwdm_;
  IPGraph* vn_;
  std::vector<std::vector<int>>* location_constraints_;
  OverlayMapping<otn_edge_map_t>* ip_otn_;
  OverlayMapping<dwdm_edge_map_t>* otn_dwdm_;

  // Input constants.
  int max_k_;
  int w_;
  int c_;

  // Derived data structures.
  std::vector<std::vector<std::vector<long>>> b_uvi_;
  std::vector<std::vector<std::vector<int>>> cost_uvi_;
  std::vector<std::vector<int>> m_p_k_;
  std::vector<std::vector<std::vector<std::vector<std::vector<long>>>>> q_pq_kjl_;

  // Input variables.
  // Location Constraint.
  IloInt2dArray l_mu_;

  // IP node to OTN node attachment.
  IloInt2dArray tau_up_;

  // OTN node to DWDM switch attachment.
  IloInt2dArray nu_pa_;

  // Indicates if OTN link (p, q, k, j, l) already exists or not.
  IloInt5dArray omega_pq_kjl_;

  // Indicates if IP Link (u, v, i) already exists or not.
  IloInt3dArray ip_link_uvi_;

  // Decision variables.
  // x: VLink to IP Link mapping.
  // y: VNode to IP Node mapping.
  // z: IP Link to OTN Link mapping.
  // zeta: Establishment of new OTN link.
  // phi: Allocation of wavelength to an OTN link.
  // psi: Allocation of wavelength on fiber.
  IloIntVar5dArray x_mn_uvi_;
  IloIntVar2dArray y_mu_;
  IloIntVar3dArray gamma_uvi_;
  IloIntVar8dArray z_uvi_pqkjl_;
  IloIntVar5dArray zeta_pq_kjl_;
  IloIntVar6dArray phi_pqkjl_l_;
  IloIntVar8dArray psi_pqkjl_abl_;

  // Objective function.
  IloIntExpr objective_;
};

#endif  // CPLEX_SOLVER_H_
