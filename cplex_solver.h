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

typedef IloArray<IloIntArray> IloInt2dArray;
typedef IloArray<IloInt2dArray> IloInt3dArray;
typedef IloArray<IloInt3dArray> IloInt4dArray;

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
  IloIntVar7dArray& z_uvi_pqkj() { return z_uvi_pqkj_; }
  IloIntVar4dArray& zeta_pq_kj() { return zeta_pq_kj_; }
  IloIntVar5dArray& phi_pqkj_l() { return phi_pqkj_l_; }
  IloIntVar7dArray& psi_pqkj_abl() { return psi_pqkj_abl_; }
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
  std::vector<std::vector<std::vector<int>>> m_pq_k_;

  // Input variables.
  IloInt2dArray l_mu_;
  IloInt2dArray tau_up_;
  IloInt4dArray omega_pq_kj_;
  IloInt3dArray ip_link_uvi_;

  // Decision variables.
  // x: VLink to IP Link mapping.
  // y: VNode to IP Node mapping.
  // z: IP Link to OTN Link mapping.
  // zeta: Activation of ODU module.
  // phi: Allocation of wavelength to an ODU module.
  // psi: Allocation of wavelength on fiber.
  IloIntVar5dArray x_mn_uvi_;
  IloIntVar2dArray y_mu_;
  IloIntVar3dArray gamma_uvi_;
  IloIntVar7dArray z_uvi_pqkj_;
  IloIntVar4dArray zeta_pq_kj_;
  IloIntVar5dArray phi_pqkj_l_;
  IloIntVar7dArray psi_pqkj_abl_;

  // Objective function.
  IloIntExpr objective_;
};

#endif  // CPLEX_SOLVER_H_
