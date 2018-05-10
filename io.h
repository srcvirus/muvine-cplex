#ifndef IO_H_
#define IO_H_

#include "datastructure.h"
#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <memory>
#include <string>

using std::unique_ptr;

typedef std::vector<std::vector<std::string>> csv_vector_t;
typedef unique_ptr<csv_vector_t> csv_vector_ptr_t;

unique_ptr<std::map<std::string, std::string>> ParseArgs(int argc,
                                                         char* argv[]) {
  unique_ptr<std::map<std::string, std::string>> arg_map(
      new std::map<std::string, std::string>());
  for (int i = 1; i < argc; ++i) {
    char* key = strtok(argv[i], "=");
    char* value = strtok(NULL, "=");
    DEBUG(" [%s] => [%s]\n", key, value);
    arg_map->insert(std::make_pair(key, value));
  }
  return std::move(arg_map);
}

csv_vector_ptr_t ReadCSVFile(const char* filename,
                             const char* delim = ",\n\r") {
  DEBUG("[Parsing %s]\n", filename);
  FILE* file_ptr = fopen(filename, "r");
  if (!file_ptr) {
    DEBUG("Invalid file %s\n", filename);
    return nullptr;
  }
  const static int kBufferSize = 1024;
  char line_buffer[kBufferSize];
  csv_vector_ptr_t ret_vector(new csv_vector_t());
  std::vector<std::string> current_line;
  int row_number = 0;
  while (fgets(line_buffer, kBufferSize, file_ptr)) {
    DEBUG("Read %d characters\n", strlen(line_buffer));
    if (strlen(line_buffer) <= 0) continue;
    if (line_buffer[0] == '\n' || line_buffer[0] == '\r') continue;
    current_line.clear();
    char* token = strtok(line_buffer, delim);
    current_line.push_back(token);
    while ((token = strtok(NULL, delim))) {
      current_line.push_back(token);
    }
    ret_vector->push_back(current_line);
  }
  fclose(file_ptr);
  DEBUG("Parsed %d lines\n", static_cast<int>(ret_vector->size()));
  return std::move(ret_vector);
}

unique_ptr<configuration> ParseConfigurationFile(const char* filename) {
  csv_vector_ptr_t csv_vector = ReadCSVFile(filename, "=\n\r");
  if (csv_vector.get() == NULL) return nullptr;
  unique_ptr<configuration> config(new configuration());
  for (int i = 0; i < csv_vector->size(); ++i) {
    const auto& row = csv_vector->at(i);
    string key = row[0];
    string value = row[1];
    if (key == "otn_topology_file") {
      config->otn_topology_file = value;
    } else if (key == "ip_topology_file") {
      config->ip_topology_file = value;
    } else if (key == "dwdm_topology_file") {
      config->dwdm_topology_file = value;
    } else if (key == "vn_topology_file") {
      config->vn_topology_file = value;
    } else if (key == "vn_location_file") {
      config->vn_location_file = value;
    } else if (key == "ip_link_mapping_file") {
      config->ip_link_mapping_file = value;
    } else if (key == "ip_port_info_file") {
      config->ip_port_info_file = value;
    } else if (key == "ip_node_mapping_file") {
      config->ip_node_mapping_file = value;
    } else if (key == "otn_link_mapping_file") {
      config->otn_link_mapping_file = value;
    } else if (key == "otn_node_mapping_file") {
      config->otn_node_mapping_file = value;
    }
  }
  return std::move(config);
}

unique_ptr<IPGraph> InitializeIPGraphFromFile(const char* filename) {
  int node_count = 0, edge_count = 0;
  csv_vector_ptr_t csv_vector = ReadCSVFile(filename);
  if (csv_vector.get() == NULL) {
    return nullptr;
  }
  unique_ptr<IPGraph> graph(new IPGraph());
  node_count = std::stoi(csv_vector->at(0)[0].c_str());
  edge_count = std::stoi(csv_vector->at(0)[1].c_str());
  graph->SetNodeCount(node_count);
  for (int i = 1; i < csv_vector->size(); ++i) {
    const auto& row = csv_vector->at(i);

    // Each line has the following format:
    // SourceID, DestinationID, Cost, Bandwidth.
    int u = std::stoi(row[0]);
    int v = std::stoi(row[1]);
    int cost = std::stoi(row[2]);
    long bw = std::stol(row[3]);

    DEBUG("Line[%d]: u = %d, v = %d, cost = %d, bw = %ld\n", i, u, v, cost, bw);
    graph->AddEdge(u, v, bw, cost);
  }
  return std::move(graph);
}

unique_ptr<OTNGraph> InitializeOTNGraphFromFile(const char* filename) {
  int node_count = 0, edge_count = 0, num_module_types = 0;
  csv_vector_ptr_t csv_vector = ReadCSVFile(filename);
  if (csv_vector.get() == NULL) return nullptr;
  unique_ptr<OTNGraph> graph(new OTNGraph());

  // The first line is formatted as:
  // num_nodes,num_links,k,<k entries for module capacities>,<k entries for
  // module cost>.
  int row_offset = 0;
  node_count = std::stoi(csv_vector->at(0)[0].c_str());
  edge_count = std::stoi(csv_vector->at(0)[1].c_str());
  graph->SetNodeCount(node_count);
  num_module_types = std::stoi(csv_vector->at(0)[2].c_str());
  for (int i = 0; i < num_module_types; ++i) {
    int capacity_k = std::stoi(csv_vector->at(row_offset)[3 + i].c_str());
    int cost_k = std::stoi(csv_vector->at(row_offset)[3 + i + num_module_types].c_str());
    graph->interface_info()->emplace_back(capacity_k, cost_k);
  }
  ++row_offset;

  graph->interfaces_installed()->resize(
      node_count, std::vector<int>(num_module_types, 0));
  // The first line is followed by num_node lines, each line is formatted as
  // <node_id>,<num_intfs_type_0>,<num_intfs_type1>,...,<num_intfs_type_k>.
  for (int i = row_offset; i < row_offset + node_count; ++i) {
    const auto& row = csv_vector->at(i);
    int node_id = std::stoi(row[0].c_str());
    for (int k = 1; k < row.size(); ++k) {
      int num_intfs_k = std::stoi(row[k].c_str());
      graph->interfaces_installed()->at(node_id)[k - 1] = num_intfs_k;
    }
  }
  row_offset += node_count;

  for (int i = row_offset; i < csv_vector->size(); ++i) {
    const auto& row = csv_vector->at(i);

    // Each line has the following format:
    // SourceID,DestID,k,src_intfs_idx,dst_intfs_idx,bw,cost.
    int p = std::stoi(row[0].c_str());
    int q = std::stoi(row[1].c_str());
    int k = std::stoi(row[2].c_str());
    int j = std::stoi(row[3].c_str());
    int l = std::stoi(row[4].c_str());
    long bw = std::stol(row[5].c_str());
    int c = std::stoi(row[6].c_str());
    graph->AddEdge(p, q, k, j, l, bw, c);
  }
  return std::move(graph);
}

unique_ptr<DWDMGraph> InitializeDWDMGraphFromFile(const char* filename) {
  int node_count = 0, edge_count = 0, num_wavelengths = 0;
  long wavelength_capacity = 0;
  csv_vector_ptr_t csv_vector = ReadCSVFile(filename);
  if (csv_vector.get() == NULL) return nullptr;
  unique_ptr<DWDMGraph> graph(new DWDMGraph());
  node_count = std::stoi(csv_vector->at(0)[0].c_str());
  edge_count = std::stoi(csv_vector->at(0)[1].c_str());

  // The first line is formatted as:
  // num_nodes,num_link,num_wl_per_fiber,wl_capacity.
  num_wavelengths = std::stoi(csv_vector->at(0)[2].c_str());
  wavelength_capacity = std::stoi(csv_vector->at(0)[3].c_str());
  std::vector<long> lambda_residual_capacity(num_wavelengths, 
      wavelength_capacity);
  graph->num_wavelengths() = num_wavelengths;
  graph->wavelength_capacity() = wavelength_capacity;
  graph->SetNodeCount(node_count);

  for (int i = 1; i < csv_vector->size(); ++i) {
    // Each line has the following format:
    // SourceID,DestID,Unit cost of using a wavelength on a fiber,a vector
    // contaiing residual capacity of each wavelength on the fiber.
    auto row = csv_vector->at(i);
    int a = std::stoi(row[0].c_str());
    int b = std::stoi(row[1].c_str());
    int c = std::stoi(row[2].c_str());
    for (int j = 3; j < 3 + num_wavelengths; ++j) {
      lambda_residual_capacity[j - 3] = std::stoi(row[j].c_str());
    }
    graph->AddEdge(a, b, lambda_residual_capacity, c);
  }
  return std::move(graph);
}

unique_ptr<std::vector<std::vector<int>>> InitializeVNLocationsFromFile(
    const char* filename, int num_virtual_nodes) {
  DEBUG("Parsing %s\n", filename);
  unique_ptr<std::vector<std::vector<int>>> ret_vector(
      new std::vector<std::vector<int>>(num_virtual_nodes));
  csv_vector_ptr_t csv_vector = ReadCSVFile(filename);
  if (csv_vector.get() == NULL) {
    return nullptr;
  }
  DEBUG("Parsing %s successful\n", filename);
  for (int i = 0; i < csv_vector->size(); ++i) {
    const std::vector<std::string>& row = csv_vector->at(i);
    int vnode_id = atoi(row[0].c_str());
    for (int j = 1; j < row.size(); ++j) {
      ret_vector->at(vnode_id).push_back(atoi(row[j].c_str()));
    }
  }
  return std::move(ret_vector);
}

unique_ptr<OverlayMapping<otn_edge_map_t>> InitializeIPOTNMappingFromFile(
    const char* nmap_file, const char* emap_file) {
  DEBUG("Reading node embedding from: %s\n", nmap_file);
  DEBUG("Reading edge embedding from: %s\n", emap_file);
  unique_ptr<OverlayMapping<otn_edge_map_t>> mapping(
      new OverlayMapping<otn_edge_map_t>());

  // Initialize node mapping.
  csv_vector_ptr_t nmap_csv_vector = ReadCSVFile(nmap_file);
  if (nmap_csv_vector.get() == NULL) {
    return nullptr;
  }
  for (int i = 0; i < nmap_csv_vector->size(); ++i) {
    const std::vector<std::string>& row = nmap_csv_vector->at(i);
    int vnode = std::stoi(row[0].c_str());
    int vnode_map = std::stoi(row[1].c_str());
    if (vnode > static_cast<int>(mapping->node_map.size()) - 1)
      mapping->node_map.resize(vnode + 1);
    mapping->node_map[vnode] = vnode_map;
  }

  // Initialize edge mapping.
  csv_vector_ptr_t emap_csv_vector = ReadCSVFile(emap_file);
  if (emap_csv_vector.get() == NULL) {
    return nullptr;
  }
  for (int i = 0; i < emap_csv_vector->size(); ++i) {
    const std::vector<std::string>& row = emap_csv_vector->at(i);
    int u = std::stoi(row[0].c_str());
    int v = std::stoi(row[1].c_str());
    int order = std::stoi(row[2].c_str());
    int p = std::stoi(row[3].c_str());
    int q = std::stoi(row[4].c_str());
    int k = std::stoi(row[5].c_str());
    int j = std::stoi(row[6].c_str());
    int l = std::stoi(row[7].c_str());
    ip_edge_t overlay_link(u, v, order);
    otn_edge_t underlay_link(p, q, k, j, l);
    if (mapping->edge_map.find(overlay_link) == mapping->edge_map.end()) {
      mapping->edge_map[overlay_link] = otn_path_t();
    }
    mapping->edge_map[overlay_link].push_back(underlay_link);
    DEBUG("Current embedding path length of (%d, %d, %d) is %u\n", u, v, order,
          mapping->edge_map[overlay_link].size());
  }
  DEBUG("Embedding of %d links read successfully\n", mapping->edge_map.size());
  return std::move(mapping);
}

unique_ptr<OverlayMapping<dwdm_edge_map_t>> InitializeOTNDWDMMappingFromFile(
    const char* nmap_file, const char* emap_file) {
  DEBUG("Reading node embedding from: %s\n", nmap_file);
  DEBUG("Reading edge embedding from: %s\n", emap_file);
  unique_ptr<OverlayMapping<dwdm_edge_map_t>> mapping(
      new OverlayMapping<dwdm_edge_map_t>());
  csv_vector_ptr_t nmap_csv_vector = ReadCSVFile(nmap_file);
  if (nmap_csv_vector.get() == NULL) return nullptr;
  for (int i = 0; i < nmap_csv_vector->size(); ++i) {
    const auto& row = nmap_csv_vector->at(i);
    int vnode = std::stoi(row[0].c_str());
    int vnode_map = std::stoi(row[1].c_str());
    if (vnode > static_cast<int>(mapping->node_map.size()) - 1)
      mapping->node_map.resize(vnode + 1);
    mapping->node_map[vnode] = vnode_map;
  }
  csv_vector_ptr_t emap_csv_vector = ReadCSVFile(emap_file);
  for (int i = 0; i < emap_csv_vector->size(); ++i) {
    const auto& row = emap_csv_vector->at(i);
    int p = std::stoi(row[0].c_str());
    int q = std::stoi(row[1].c_str());
    int k = std::stoi(row[2].c_str());
    int j = std::stoi(row[3].c_str());
    int l = std::stoi(row[4].c_str());
    int a = std::stoi(row[5].c_str());
    int b = std::stoi(row[6].c_str());
    int ll = std::stoi(row[7].c_str());
    otn_edge_t overlay_link(p, q, k, j, l);
    dwdm_edge_t underlay_link(a, b, l);
    if (mapping->edge_map.find(overlay_link) == mapping->edge_map.end()) {
      mapping->edge_map[overlay_link] = dwdm_path_t();
    }
    mapping->edge_map[overlay_link].push_back(underlay_link);
    DEBUG("Current embedding path length of (%d, %d, %d, %d, %d) is %u\n", p, q, k,
          j, l, mapping->edge_map[overlay_link].size());
  }
  DEBUG("Embedding of %d links read successfully\n", mapping->edge_map.size());
  return std::move(mapping);
}

unique_ptr<std::vector<std::vector<int>>> InitializePortInfoFromFile(
    const char* port_info_file) {
  unique_ptr<std::vector<std::vector<int>>> port_info(
      new std::vector<std::vector<int>>());
  csv_vector_ptr_t csv_vector = ReadCSVFile(port_info_file);
  for (int i = 0; i < csv_vector->size(); ++i) {
    const std::vector<std::string>& row = csv_vector->at(i);
    int u = std::stoi(row[0].c_str());
    int num_ports = std::stoi(row[1].c_str());
    int port_capacity = std::stoi(row[2].c_str());
    std::vector<int> v;
    v.push_back(u);
    v.push_back(num_ports);
    v.push_back(port_capacity);
    port_info->push_back(v);
  }
  return std::move(port_info);
}

#endif  // IO_H_
