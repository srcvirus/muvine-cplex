#ifndef DATASTRUCTURE_H_
#define DATASTRUCTURE_H_

#include <math.h>
#include <stdlib.h>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#define INF 99999999
#define MAXN 1000
#define MAX_PARALLEL_LINKS 10
#define NIL -1

using std::map;
using std::string;
using std::unique_ptr;
using std::vector;

// Configuration parameters for the optimizer.
typedef struct configuration {
  string otn_topology_file;
  string ip_topology_file;
  string dwdm_topology_file;
  string vn_topology_file;
  string vn_location_file;
  string ip_link_mapping_file;
  string ip_port_info_file;
  string ip_node_mapping_file;
  string otn_link_mapping_file;
  string otn_node_mapping_file;
} configuration;

// Wrapper for a matrix datastructure.
template <typename T>
struct matrix_t {
  std::vector<std::vector<T>> matrix;
  matrix_t() {}
  matrix_t(int rows, int columns, T fill_value = T())
      : matrix(rows, std::vector<T>(columns, fill_value)) {}
  const int num_rows() const { return this->matrix.size(); }
  const int num_columns() const { 
    this->matrix.size() > 0 ? this->matrix[0].size() : 0; 
  }
};

// A generic edge endpoint. At a bare minimum an edge endpoint consists of the
// node_id_ of the endpoint.
class EdgeEndpoint {
 public:
  EdgeEndpoint() : node_id_(NIL) {}
  EdgeEndpoint(int node_id) : node_id_(node_id) {}
  const int node_id() const { return node_id_; }

 protected:
  int node_id_;
};

// IP network specific edge endpoint that extends EdgeEndpoint Class.
class IPEdgeEndpoint : public EdgeEndpoint {
 public:
  IPEdgeEndpoint()
      : EdgeEndpoint(),
        order_(NIL),
        bandwidth_(0),
        residual_bandwidth_(0),
        cost_(INF) {}
  IPEdgeEndpoint(int n, int o, long bw, int c)
      : EdgeEndpoint(n),
        order_(o),
        bandwidth_(bw),
        residual_bandwidth_(bw),
        cost_(c) {}
  int& order() { return order_; }
  long& bandwidth() { return bandwidth_; }
  long& residual_bandwidth() { return residual_bandwidth_; }
  int& cost() { return cost_; }
  string GetDebugString() const { return ""; }

 private:
  // In case of parallel links order_ differentiates multiple links between
  // the same pair of nodes.
  int order_;

  // Bandwidth of the link.
  long bandwidth_;

  // Remaining bandwidth on the link.
  long residual_bandwidth_;

  // Cost of allocating unit bandwidth from the link.
  int cost_;
};

// OTN specific edge endpoint that extends EdgeEndpoint class. It contains
// necessary metadata to indentify the multiplexing units.
class OTNEdgeEndpoint : public EdgeEndpoint {
 public:
  OTNEdgeEndpoint() : EdgeEndpoint() {}
  OTNEdgeEndpoint(int n, int k, int source_idx, int dst_idx, long bw, int c)
      : EdgeEndpoint(n), 
      intf_type_(k), 
      src_intf_(source_idx), 
      dst_intf_(dst_idx), 
      bandwidth_(bw), 
      residual_bandwidth_(bw),
      cost_(c) {}

  int& intf_type() { return intf_type_; }
  int& src_intf() { return src_intf_; }
  int& dst_intf() { return dst_intf_; }
  long& bandwidth() { return bandwidth_; }
  long& residual_bandwidth() { return residual_bandwidth_; }
  int& cost() { return cost_; }

 private:
  // Type of the OTN interface used to establish the OTN link.
  int intf_type_;

  // Index of the source interface.
  int src_intf_;

  // Index of the destination interface.
  int dst_intf_;

  // Max bandwidth capacity of the link.
  long bandwidth_;

  // Residual bandwidth of this OTN link.
  long residual_bandwidth_;

  // Cost of allocating unit bandwidth.
  int cost_;
};

// DWDM network specific edge endpoint that extends EdgeEndpoint class. Contains
// necessary metadata to identify the allocated wavelengths from a fiber.
class DWDMEdgeEndpoint : public EdgeEndpoint {
 public:
  DWDMEdgeEndpoint() : EdgeEndpoint() {}
  DWDMEdgeEndpoint(int node_id, const std::vector<long>& lambda_res_bw, int c)
      : EdgeEndpoint(node_id), 
      lambda_residual_bandwidth_(lambda_res_bw), 
      cost_(c) {}
  std::vector<long>& lambda_residual_bandwidth() { return lambda_residual_bandwidth_; }
  int& cost() { return cost_; }

 private:
  // lambda_residual_bandwidth[l] is the residual bandwidth for wavelength l
  // that can be used to multiplex more OTN connections.
  std::vector<long> lambda_residual_bandwidth_;

  // Cost of using a wavelength from this DWDM link.
  int cost_;
};

// Represents an IP network where the only resources are bandwidth. In the
// optimizer this class is used to instantiate the IP network and also the VN
// request.
class IPGraph {
 public:
  IPGraph() {
    adj_list_ = unique_ptr<std::vector<std::vector<IPEdgeEndpoint>>>(
        new std::vector<std::vector<IPEdgeEndpoint>>);
    port_counts_ = unique_ptr<std::vector<int>>(new std::vector<int>());
    port_capacities_ = unique_ptr<std::vector<int>>(new std::vector<int>());
    node_count_ = edge_count_ = 0;
  }

  // Accessor methods.
  int node_count() { return node_count_; }
  int edge_count() { return edge_count_; }
  const std::vector<std::vector<IPEdgeEndpoint>>* adj_list() const {
    return static_cast<const std::vector<std::vector<IPEdgeEndpoint>>*>(
        adj_list_.get());
  }

  void SetNodeCount(int node_count) {
    this->node_count_ = node_count;
    adj_list_->resize(node_count);
  }

  // u and v are 0-based identifiers of an edge endpoint. An edge is
  // bi-directional, i.e., calling Graph::AddEdge with u = 1, v = 3 will add
  // both (1, 3) and (3, 1) in the graph.
  void AddEdge(int u, int v, long bw, int cost) {
    if (adj_list_->size() < u + 1) adj_list_->resize(u + 1);
    if (adj_list_->size() < v + 1) adj_list_->resize(v + 1);
    int order = 0;
    auto& neighbors = adj_list_->at(u);
    for (auto& end_point : neighbors) {
      if (end_point.node_id() == v) ++order;
    }
    adj_list_->at(u).push_back(IPEdgeEndpoint(v, order, bw, cost));
    adj_list_->at(v).push_back(IPEdgeEndpoint(u, order, bw, cost));
    ++edge_count_;
    node_count_ = adj_list_->size();
  }

  int GetEdgeCost(int u, int v, int order = 0) const {
    auto& neighbors = adj_list_->at(u);
    for (auto end_point : neighbors) {
      if (end_point.node_id() == v && end_point.order() == order)
        return end_point.cost();
    }
    return -1;
  }

  long GetEdgeBandwidth(int u, int v, int order = 0) const {
    auto& neighbors = adj_list_->at(u);
    for (auto end_point : neighbors) {
      if (end_point.node_id() == v && end_point.order() == order)
        return end_point.bandwidth();
    }
    return -1;
  }

  void SetEdgeBandwidth(int u, int v, long bw, int order = 0) {
    auto& neighbors = adj_list_->at(u);
    for (auto end_point_it = neighbors.begin(); end_point_it != neighbors.end();
         ++end_point_it) {
      if (end_point_it->node_id() == v && end_point_it->order() == order) {
        end_point_it->bandwidth() = bw;
        break;
      }
    }
  }

  long GetEdgeResidualBandwidth(int u, int v, int order = 0) const {
    auto& neighbors = adj_list_->at(u);
    for (auto end_point : neighbors) {
      if (end_point.node_id() == v && end_point.order() == order)
        return end_point.residual_bandwidth();
    }
    return -1;
  }

  void SetEdgeResidualBandwidth(int u, int v, long rbw, int order = 0) {
    auto& neighbors = adj_list_->at(u);
    for (auto end_point_it = neighbors.begin(); end_point_it != neighbors.end();
         ++end_point_it) {
      if (end_point_it->node_id() == v && end_point_it->order() == order) {
        end_point_it->residual_bandwidth() = rbw;
        break;
      }
    }
  }

  // Returns the cumulative bandwidth of all links originating at node u.
  long GetTotalNodeBandwidth(int u) const {
    auto& neighbors = adj_list_->at(u);
    long total_bw = 0;
    for (auto end_point : neighbors) {
      total_bw += end_point.bandwidth();
    }
    return total_bw;
  }

  inline int GetNodeDegree(int u) const { return adj_list_->at(u).size(); }
  inline int GetPortCapacity(int u) const { return port_capacities_->at(u); }
  inline int GetPortCount(int u) const { return port_counts_->at(u); }
  inline int GetResidualPortCount(int u) const {
    return port_counts_->at(u) - GetNodeDegree(u);
  }

  void SetPortCount(int u, int port_count) {
    if (port_counts_->size() <= u) {
      port_counts_->resize(u + 1);
    }
    port_counts_->at(u) = port_count;
  }

  void SetPortCapacity(int u, int port_capacity) {
    if (port_capacities_->size() <= u) {
      port_capacities_->resize(u + 1);
    }
    port_capacities_->at(u) = port_capacity;
  }

  virtual ~IPGraph() { adj_list_.reset(); }

 private:
  unique_ptr<std::vector<std::vector<IPEdgeEndpoint>>> adj_list_;
  int node_count_, edge_count_, total_port_count_;
  unique_ptr<std::vector<int>> port_capacities_;
  unique_ptr<std::vector<int>> port_counts_;
};

// Represents the OTN network.
class OTNGraph {
 public:
  OTNGraph() {
    adj_list_ = unique_ptr<std::vector<std::vector<OTNEdgeEndpoint>>>(
        new std::vector<std::vector<OTNEdgeEndpoint>>());
    interface_info_ = unique_ptr<std::vector<interface_metadata>>(
        new std::vector<interface_metadata>());
    interfaces_installed_ = 
      unique_ptr<std::vector<std::vector<int>>>(new std::vector<std::vector<int>>());
    node_count_ = edge_count_ = 0;
  }

  typedef struct interface_metadata {
    long capacity;
    int cost;
    interface_metadata() : capacity(0), cost(INF) {}
    interface_metadata(int cap, int c): capacity(cap), cost(c) {}
  } interface_metadata;

  // Accessor methods.
  int node_count() { return node_count_; }
  int edge_count() { return edge_count_; }
  std::vector<std::vector<OTNEdgeEndpoint>>* adj_list() {
    return static_cast<std::vector<std::vector<OTNEdgeEndpoint>>*>(
        adj_list_.get());
  }
  std::vector<interface_metadata>* interface_info() {
    return static_cast<std::vector<interface_metadata>*>(interface_info_.get());
  }

  std::vector<std::vector<int>>* interfaces_installed() {
    return static_cast<std::vector<std::vector<int>>*>(
        interfaces_installed_.get());
  }

  int GetInterfaceCount(int node_id, int k) const {
    return this->interfaces_installed_->at(node_id)[k];
  }

  void SetNodeCount(int node_count) {
    this->node_count_ = node_count;
    this->adj_list_->resize(node_count);
  }

  // u and v are 0-based identifiers of an edge endpoint. An edge is
  // bi-directional, i.e., calling Graph::AddEdge with u = 1, v = 3 will add
  // both (1, 3) and (3, 1) in the graph.
  bool AddEdge(int u, int v, int k, int src_idx, int dst_idx, long bw, 
      int cost) {
    if (adj_list_->size() < u + 1) adj_list_->resize(u + 1);
    if (adj_list_->size() < v + 1) adj_list_->resize(v + 1);
    auto& neighbors = adj_list_->at(u);
    for (auto end_point : neighbors) {
      if (end_point.node_id() == v) return false;
    }
    adj_list_->at(u).push_back(
        OTNEdgeEndpoint(v, k, src_idx, dst_idx, bw, cost));
    adj_list_->at(v).push_back(
        OTNEdgeEndpoint(u, k, dst_idx, src_idx, bw, cost));
    ++edge_count_;
    node_count_ = adj_list_->size();
    return true;
  }

  int GetInterfaceCost(int interface_type) const {
    if (interface_type >= interface_info_->size()) return INF;
    return interface_info_->at(interface_type).cost;
  }

  inline int GetNodeDegree(int u) const { return adj_list_->at(u).size(); }
  virtual ~OTNGraph() { 
    adj_list_.reset(); 
    interface_info_.reset();
  }

 private:
  unique_ptr<std::vector<std::vector<OTNEdgeEndpoint>>> adj_list_;
  unique_ptr<std::vector<interface_metadata>> interface_info_;
  unique_ptr<std::vector<std::vector<int>>> interfaces_installed_;
  int node_count_, edge_count_;
};

// Represents the DWDM optical network.
class DWDMGraph {
 public:
  DWDMGraph() {
    adj_list_ = unique_ptr<std::vector<std::vector<DWDMEdgeEndpoint>>>(
        new std::vector<std::vector<DWDMEdgeEndpoint>>());
    node_count_ = edge_count_ = num_wavelengths_ = wavelength_capacity_ = 0;
  }

  // Accessor methods.
  int node_count() const { return node_count_; }
  int edge_count() const { return edge_count_; }
  int& num_wavelengths() { return num_wavelengths_; }
  long& wavelength_capacity() { return wavelength_capacity_; }
  std::vector<std::vector<DWDMEdgeEndpoint>>* adj_list() {
    return static_cast<std::vector<std::vector<DWDMEdgeEndpoint>>*>(
        adj_list_.get());
  }

  void SetNodeCount(int node_count) {
    this->node_count_ = node_count;
    this->adj_list_->resize(node_count);
  }

  // u and v are 0-based identifiers of an edge endpoint. An edge is
  // bi-directional, i.e., calling Graph::AddEdge with u = 1, v = 3 will add
  // both (1, 3) and (3, 1) in the graph.
  bool AddEdge(int u, int v, const std::vector<long>& lambda_res_bw, int cost) {
    if (adj_list_->size() < u + 1) adj_list_->resize(u + 1);
    if (adj_list_->size() < v + 1) adj_list_->resize(v + 1);
    auto& neighbors = adj_list_->at(u);
    for (auto end_point : neighbors) {
      if (end_point.node_id() == v) return false;
    }
    adj_list_->at(u).push_back(DWDMEdgeEndpoint(v, lambda_res_bw, cost));
    adj_list_->at(v).push_back(DWDMEdgeEndpoint(u, lambda_res_bw, cost));
    ++edge_count_;
    node_count_ = adj_list_->size();
    return true;
  }

  bool IsWavelengthAvailable(int u, int v, int lambda, long bw_threshold = 0) {
    auto& neighbors = adj_list_->at(u);
    for (auto& end_point : neighbors) {
      if (end_point.node_id() == v) 
        return end_point.lambda_residual_bandwidth()[lambda] >= bw_threshold;
    }
    return false;
  }

  std::vector<long>& GetWavelengthResidualBw(int u, int v) {
    auto& neighbors = adj_list_->at(u);
    for (auto& end_point : neighbors) {
      if (end_point.node_id() == v) 
        return end_point.lambda_residual_bandwidth();
    }
  }

  inline int GetNodeDegree(int u) const { return adj_list_->at(u).size(); }

  virtual ~DWDMGraph() { adj_list_.reset(); }

 private:
  unique_ptr<std::vector<std::vector<DWDMEdgeEndpoint>>> adj_list_;
  int node_count_, edge_count_, num_wavelengths_;
  long wavelength_capacity_;
};

// Type definition for convenience.

// Datastructure for an IP link. It contains the endpoints (first, second) and
// order to distinguish between multiple links between the same pair of
// endpoints.
typedef struct ip_edge_t {
  int first, second, order;
  ip_edge_t(int f, int s, int o = 0) : first(f), second(s), order(o) {}
  bool operator<(const ip_edge_t& e) const {
    if (first != e.first) return first < e.first;
    if (second != e.second) return second < e.second;
    return order < e.order;
  }
  bool operator==(const ip_edge_t& e) const {
    return first == e.first && second == e.second && order == e.order;
  }
} ip_edge_t;

// Type definitions for an IP path and an IP link to an IP path mapping. The
// latter is used for virtual link mapping.
typedef std::vector<ip_edge_t> ip_path_t;
typedef std::map<ip_edge_t, ip_path_t> ip_edge_map_t;

// Datastructure for an OTN link. first and second represent the endpoint OTN
// nodes. Moreover, we also identify the specific interface type and interface
// instance by interface_type and interface_instance, respectively.
typedef struct otn_edge_t {
  int first, second, interface_type, src_idx, dst_idx;
  otn_edge_t(int f, int s, int k, int j, int l)
      : first(f), second(s), interface_type(k), src_idx(j), dst_idx(l) {}
  bool operator<(const otn_edge_t& e) const {
    if (first != e.first) return first < e.first;
    if (second != e.second) return second < e.second;
    if (interface_type != e.interface_type) return interface_type < e.interface_type;
    if (src_idx != e.src_idx) return src_idx < e.src_idx;
    return dst_idx < e.dst_idx;
  }
  bool operator==(const otn_edge_t& e) const {
    return first == e.first && second == e.second &&
           interface_type == e.interface_type && 
           src_idx == e.src_idx && dst_idx == e.dst_idx;
  }
} otn_edge_t;

// The following two are used for representing the mapping between an IP link
// and an OTN path. An IP link is mapped to a sequence of OTN nodes and interfaces.
typedef std::vector<otn_edge_t> otn_path_t;
typedef std::map<ip_edge_t, otn_path_t> otn_edge_map_t;

// Datastructure to represent a DWDM connection between two optical nodes first,
// second over a wavelength lambda.
typedef struct dwdm_edge_t {
  int first, second, lambda;
  dwdm_edge_t(int f, int s, int l) : first(f), second(s), lambda(l) {}
  bool operator<(const dwdm_edge_t& e) const {
    if (first != e.first) return first < e.first;
    if (second != e.second) return second < e.second;
    return lambda < e.lambda;
  }
  bool operator==(const dwdm_edge_t& e) const {
    return first == e.first && second == e.second && lambda == e.lambda;
  }
} dwdm_edge_t;

// The following type defenitions are used to represent the mapping between an
// OTN connection (i.e., a connection between two interfaces of a specific type)
// over a DWDM path with continous wavelength.
typedef std::vector<dwdm_edge_t> dwdm_path_t;
typedef std::map<otn_edge_t, dwdm_path_t> dwdm_edge_map_t;

// The following template class can be used to represent different types of
// mapping, e.g., VN -> IP mapping, IP->OTN Mapping, OTN->DWDM Mapping. For the
// appropriate type of mapping instantiate an object with appropriate EdgeMapT.
// For instance the following will create an IP to OTN mapping:
//      OverlayMapping<otn_edge_map_t> ip_otn_mapping;
template <class EdgeMapT>
class OverlayMapping {
 public:
  std::vector<int> node_map;
  EdgeMapT edge_map;
  OverlayMapping() {}
  OverlayMapping(const std::vector<int>& nmap, const EdgeMapT& emap)
      : node_map(nmap), edge_map(emap) {}
};

#endif  // DATASTRUCTURE_H_
