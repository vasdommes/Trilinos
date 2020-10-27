#ifndef _ZOLTAN2_IANGRAPHADAPTER_HPP_
#define _ZOLTAN2_IANGRAPHADAPTER_HPP_

#include <Zoltan2_GraphAdapter.hpp>
#include "dist_graph.h"
#include <iostream>
namespace Zoltan2{

template <typename User, typename UserCoord=User>
  class IanGraphAdapter :public GraphAdapter<User> {
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::node_t   node_t;
  //this typedef is slow, use the comment for production
  typedef typename InputTraits<User>::offset_t offset_t;
  //typedef uint64_t offset_t;//lno_t offset_t;
  typedef User user_t;
  typedef UserCoord userCoord_t;
  typedef GraphAdapter<User, UserCoord> base_adapter_t;
#endif

  virtual ~IanGraphAdapter(){
    delete [] dist_graph.out_edges;
    delete [] dist_graph.out_offsets;
    delete [] dist_graph.local_unmap;
    delete [] dist_graph.ghost_unmap;
    delete [] dist_graph.ghost_tasks;
  }
  //take a graph as an argument
  IanGraphAdapter(color_dist_graph_t* dist_graph_) : GraphAdapter<User>(){
    //std::cout<<"Starting IanGraphAdapter constructor\n";
    std::cout<<"size of offset_t = "<<sizeof(offset_t)<<"\n";
    //dist_graph = new pdist_graph_t;
    dist_graph.n = dist_graph_->n;
    dist_graph.m = dist_graph_->m;
    dist_graph.n_local = dist_graph_->n_local;
    dist_graph.m_local = dist_graph_->m_local;
    dist_graph.n_offset = dist_graph_->n_offset;
    if(dist_graph.n_offset == dist_graph.n_local) dist_graph.n_offset++;
    dist_graph.n_ghost = dist_graph_->n_ghost;
    dist_graph.n_total = dist_graph_->n_total;
    dist_graph.out_edges = new gno_t[dist_graph.m_local];
    std::cout<<"allocated out_edges, "<<dist_graph.m_local<<" entries\n";
    for(int i = 0; i < dist_graph.m_local; i++) {
      if(dist_graph_->out_edges[i] < dist_graph_->n_local){
        dist_graph.out_edges[i] = dist_graph_->local_unmap[dist_graph_->out_edges[i]];
      } else {
        dist_graph.out_edges[i] = dist_graph_->ghost_unmap[dist_graph_->out_edges[i] - dist_graph_->n_local];
      }
    }
    dist_graph.out_offsets = new offset_t[dist_graph.n_local+1];
    std::cout<<"allocated out_offsets, "<<dist_graph.n_local+1<<" entries\n";
    for(int i = 0; i < dist_graph.n_local+1; i++) dist_graph.out_offsets[i] = dist_graph_->out_offsets[i];
    dist_graph.out_offsets[dist_graph.n_offset-1] = dist_graph.m_local; 
    std::cout<<"last offsets index = "<<dist_graph.out_offsets[dist_graph.n_local]<<"\n";
    dist_graph.local_unmap = new gno_t[dist_graph.n_local];
    std::cout<<"allocated local_unmap, "<<dist_graph.n_local<<" entries\n";
    for(int i = 0; i < dist_graph.n_local; i++) dist_graph.local_unmap[i] = dist_graph_->local_unmap[i];
    dist_graph.ghost_unmap = new gno_t[dist_graph.n_ghost];
    std::cout<<"allocated ghost_unmap, "<<dist_graph.n_ghost<<" entries\n";
    for(int i = 0; i < dist_graph.n_ghost; i++) dist_graph.ghost_unmap[i] = dist_graph_->ghost_unmap[i];
    dist_graph.ghost_tasks = new gno_t[dist_graph.n_ghost];
    std::cout<<"allocated ghost_tasks, "<<dist_graph.n_ghost<<" entries\n";
    for(int i = 0; i < dist_graph.n_ghost; i++) dist_graph.ghost_tasks[i] = dist_graph_->ghost_tasks[i];

    /*std::cout<<"IanGraphAdapter: n_local="<<dist_graph.n_local<<" m_local="<<dist_graph.m_local<<" n_offset="<<dist_graph.n_offset<<"\n";
    std::cout<<"IanGraphAdapter: out_offsets: ";
    for(int i = 0; i < dist_graph.n_offset; i++) std::cout<<dist_graph.out_offsets[i]<<" ";
    std::cout<<"\n";
    std::cout<<"IanGraphAdapter: out_edges: ";
    for(int i = 0; i < dist_graph.m_local; i++) std::cout<<dist_graph.out_edges[i]<<" ";
    std::cout<<"\n";
    std::cout<<"Ending IanGraphAdapter constructor\n";*/
  }

  size_t getLocalNumVertices() const {
    //return local # verts using my datastructures
    //std::cout<<"IanGraphAdapter::getLocalNumVertices returning "<<dist_graph.n_local<<"\n";
    return (size_t) dist_graph.n_local;
  }

  size_t getLocalNumEdges() const {
    //return local # edges using my datastructures
    //std::cout<<"IanGraphAdapter::getLocalNumEdges returning "<<dist_graph.m_local<<"\n";
    return (size_t) dist_graph.m_local;
  }

  void getVertexIDsView(const gno_t *&vertexIds) const{
    //return set vertex IDs equal to a pointer to global IDs
    //std::cout<<"IanGraphAdapter::getVertexIDsView setting argument to dist_graph.local_unmap\n";
    vertexIds = dist_graph.local_unmap;
  }
 
  int getNumWeightsPerEdge() const{
    return 0;
  } 
  void getEdgesView(const offset_t *&offsets,
                           const gno_t *&adjIds)  const{
    /*gno_t* global_adjs = new gno_t[dist_graph.m_local];
    for(int i = 0; i < dist_graph.m_local; i++) {
      if(dist_graph.out_edges[i] < dist_graph.n_local){
        global_adjs[i] = dist_graph.local_unmap[dist_graph.out_edges[i]];
      } else {
        global_adjs[i] = dist_graph.ghost_unmap[dist_graph.out_edges[i] - dist_graph.n_local];
      }
    }
    adjIds = global_adjs;
    offsets = dist_graph.out_offsets;*/
    //std::cout<<"IanGraphAdapter::getEdgesView setting arguments to dist_graph.out_offsets and dist_graph.out_edges\n";
    offsets = dist_graph.out_offsets;
    adjIds = dist_graph.out_edges;
  }

private:
  //add graph datastructure here
  //color_dist_graph_t* dist_graph;
  struct pdist_graph_t{
    gno_t n;
    gno_t m;
    
    gno_t n_local;
    gno_t m_local;
    
    gno_t n_offset;
    gno_t n_ghost;
    gno_t n_total;
    
    gno_t* out_edges;
    offset_t* out_offsets;
    
    gno_t* local_unmap;
    gno_t* ghost_unmap;
    gno_t* ghost_tasks;
  } dist_graph;
  
};

}
#endif
