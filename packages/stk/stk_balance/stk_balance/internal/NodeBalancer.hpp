#ifndef NODEBALANCER_HPP
#define NODEBALANCER_HPP

#include <set>
#include <map>
#include <vector>

#include "stk_mesh/base/Types.hpp"

namespace stk {
namespace balance {
namespace internal {

class NodeBalancer {

public:
  NodeBalancer(stk::mesh::BulkData& bulk);

  bool balance_node_entities(const double targetLoadBalance,
                             const unsigned maxIterations);

private:
  void getInterfaceDescription(std::set<int>& neighborProcessors,
                               std::map<stk::mesh::Entity, std::vector<int>>& interfaceNodesAndProcessors);

  void getGlobalLoadImbalance(double &loadFactor,
                              int& numLocallyOwnedNodes);

  void exchangeLocalSizes(const std::set<int>& neighborProcessors,
                          int& numLocallyOwnedNodes,
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
                          std::map<int, int>& numLocallyOwnedByRank);
#else
                          std::map<>& numLocallyOwnedByRank);
#endif

  void changeOwnersOfNodes(const std::map<stk::mesh::Entity, std::vector<int> >& interfaceNodesAndProcessors,
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
                           std::map<int, int>& numLocallyOwnedByRank,
#else
                           std::map<>& numLocallyOwnedByRank,
#endif
                           int numLocallyOwnedNodes);

  stk::mesh::BulkData & m_bulkData;
  const stk::mesh::MetaData & m_metaData;
};

}
}
}

#endif // NODEBALANCER_HPP
