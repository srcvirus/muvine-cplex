#ifndef VNE_SOLUTION_BUILDER_H_
#define VNE_SOLUTION_BUILDER_H_

#include "cplex_solver.h"

class VNESolutionBuilder {
 public:
  VNESolutionBuilder(MuViNESolver *vne_solver_ptr, IPGraph *ip_topology,
                     OTNGraph *otn_topology, DWDMGraph *dwdm_topology,
                     IPGraph *vn_topology,
                     OverlayMapping<ip_edge_map_t> *otn_link_mapping)
      : vne_solver_ptr_(vne_solver_ptr),
        ip_topology_(ip_topology),
        otn_topology_(otn_topology),
        dwdm_topology_(dwdm_topology),
        vn_topology_(vn_topology),
        otn_link_mapping_(otn_link_mapping) {}

  // Prints virtual node to IP node mapping on stdout. If filename is not NULL
  // then the same output is written to the corresponding file as well. Each
  // line in the output has the following format:
  // Virtual node <vnode_id> --> IP node <ip_node_id>
  void PrintVNodeMapping(const char *filename);

  // Prints virtual link to IP link mapping on stdout. If filename is not NULL
  // then the same output is written to the corresponding file as well. Each
  // line in the output has the following format:
  // Virtual link (<src>, <dst>) --> IP link (<src>, <dst>, <order>).
  void PrintVLinkMapping(const char *filename);

  // Prints status of running the solver, i.e., Optimal, Infeasible, etc.
  void PrintSolutionStatus(const char *filename);

  // Prints the value of objective function obtained by the solver.
  void PrintCost(const char *filename);

  // Prints new IP links and their mapping on OTN on stdout. If filename is not
  // NULL then output is written to the corresponding file as well. Each line in
  // the output has the following format:
  // New IP Link (<src>, <dst>, <order>) --> OTN link (<otn_src>, <otn_dst>,
  // <module_type>, <module_instance>).
  // Order is used to break ties for parallel links.
  void PrintNewIPLinks(const char *filename);

  // Prints the routing information for newly activated OTN modules on stdout.
  // If filename is not NULL then output is written to the corresponding file as
  // well. Each line in the output has the following format:
  // (<otn_src>, <otn_dst>, <module_type>, <module_instance>) -> (<dwdm_src>,
  // <dwdm_dst>, <lambda>)
  void PrintNewOTNModules(const char *filename);

 private:
  MuViNESolver *vne_solver_ptr_;
  IPGraph *ip_topology_;
  OTNGraph *otn_topology_;
  DWDMGraph *dwdm_topology_;
  IPGraph *vn_topology_;
  OverlayMapping<ip_edge_map_t> *otn_link_mapping_;
};

#endif  // VNE_SOLUTION_BUILDER_H_
