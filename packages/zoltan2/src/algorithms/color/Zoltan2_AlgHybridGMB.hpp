#ifndef _ZOLTAN2_ALGHYBRIDGMB_HPP_
#define _ZOLTAN2_ALGHYBRIDGMB_HPP_

#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <queue>
#include <sys/time.h>
#include <algorithm>

#include "Zoltan2_Algorithm.hpp"
#include "Zoltan2_GraphModel.hpp"
#include "Zoltan2_ColoringSolution.hpp"
#include "Zoltan2_Util.hpp"
#include "Zoltan2_TPLTraits.hpp"

#include "Tpetra_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"

#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosGraph_Distance1ColorHandle.hpp"
#include <stdlib.h>

//////////////////////////////////////////////
//! \file Zoltan2_AlgHybridGMB.hpp
//! \brief A hybrid version of the framework proposed by Gebremedhin, Manne, 
//!        and Boman

namespace Zoltan2 {

template <typename Adapter>
class AlgHybridGMB : public Algorithm<Adapter>
{
  public:
  
    using lno_t = typename Adapter::lno_t;
    using gno_t = typename Adapter::gno_t;
    using offset_t = typename Adapter::offset_t;//lno_t;//typename Adapter::offset_t;
    using scalar_t = typename Adapter::scalar_t;
    using base_adapter_t = typename Adapter::base_adapter_t;
    using map_t = Tpetra::Map<lno_t, gno_t>;
    using femv_scalar_t = int;
    using femv_t = Tpetra::FEMultiVector<femv_scalar_t, lno_t, gno_t>;
    using device_type = Tpetra::Map<>::device_type;//Kokkos::Device<Kokkos::Cuda,Kokkos::Cuda::memory_space>;
    using execution_space = Tpetra::Map<>::execution_space;//Kokkos::Cuda;
    using memory_space = Tpetra::Map<>::memory_space;//Kokkos::Cuda::memory_space;
    double timer() {
      struct timeval tp;
      gettimeofday(&tp, NULL);
      return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);
    }
    
  private:

    void buildModel(modelFlag_t &flags); 

    //function to invoke KokkosKernels distance-1 coloring    
    template <class ExecutionSpace, typename TempMemorySpace, 
              typename MemorySpace>
    void colorInterior(const size_t nVtx, 
                       Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace,MemorySpace> > adjs_view,
                       Kokkos::View<offset_t*, Kokkos::Device<ExecutionSpace, MemorySpace> > offset_view, 
                       Teuchos::ArrayRCP<int> colors,
                       Teuchos::RCP<femv_t> femv,
		       Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace, MemorySpace> > vertex_list,
		       size_t vertex_list_size = 0,
                       bool recolor=false){
      

      using KernelHandle =  KokkosKernels::Experimental::KokkosKernelsHandle
          <offset_t, lno_t, lno_t, ExecutionSpace, MemorySpace, 
           MemorySpace>;
      using lno_row_view_t = Kokkos::View<offset_t*, Kokkos::Device<ExecutionSpace, MemorySpace>>;
      using lno_nnz_view_t = Kokkos::View<lno_t*, Kokkos::Device<ExecutionSpace, MemorySpace>>;

      KernelHandle kh;

      if(recolor){
        kh.create_graph_coloring_handle(KokkosGraph::COLORING_VBBIT);
      } else {
        kh.create_graph_coloring_handle(KokkosGraph::COLORING_EB);  
      }
      if(vertex_list_size != 0){
        kh.get_graph_coloring_handle()->set_vertex_list(vertex_list,vertex_list_size);
      }
      //kh.create_graph_coloring_handle(KokkosGraph::COLORING_EB);
      kh.set_shmem_size(16128);
      kh.set_verbose(true);
      //kh.get_graph_coloring_handle()->set_eb_num_initial_colors(10);	
      //set the initial coloring of the kh.get_graph_coloring_handle() to be
      //the data view from the femv.
      Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<MemorySpace>();
      Kokkos::View<int*, Tpetra::Map<>::device_type >  sv = subview(femvColors, Kokkos::ALL, 0);
      //Kokkos::View<int*,Kokkos::Device<ExecutionSpace, MemorySpace> > color_view("Colors",sv.size());
      //Kokkos::deep_copy(color_view,sv);
      kh.get_graph_coloring_handle()->set_vertex_colors(sv);
      kh.get_graph_coloring_handle()->set_tictoc(true);
      KokkosGraph::Experimental::graph_color_symbolic<KernelHandle, lno_row_view_t, lno_nnz_view_t>
                                                     (&kh, nVtx, nVtx, offset_view, adjs_view);
       
      numColors = kh.get_graph_coloring_handle()->get_num_colors();
      //Kokkos::deep_copy(sv,color_view);
      std::cout<<"\nKokkosKernels Coloring: "<<kh.get_graph_coloring_handle()->get_overall_coloring_time()<<" iterations: "<<kh.get_graph_coloring_handle()->get_num_phases()<<"\n\n";
    }
    
    double doOwnedToGhosts(RCP<const map_t> mapOwnedPlusGhosts,
                         size_t nVtx,
			 //Teuchos::ArrayView<const offset_t> offsets,
			 //Teuchos::ArrayView<const lno_t> adjs,
                         //const std::vector<lno_t>& verts_to_send,
			 Kokkos::View<bool*,device_type> send_flags,
			 //Kokkos::View<lno_t[1],device_type> verts_to_send_size,
                         //ArrayView<int> owners,
                         Kokkos::View<int*, device_type>& colors,
                         gno_t& recv, gno_t& send){
      
      int nprocs = comm->getSize();
      int* sendcnts = new int[nprocs];
      int* recvcnts = new int[nprocs];
      for(int i = 0; i < nprocs; i++){
        sendcnts[i] = 0;
        recvcnts[i] = 0;
      }
      
      for(size_t i = 0; i < exportLIDs.size(); i++){
        if(send_flags(exportLIDs[i])){
	  sendcnts[exportPIDs[i]] += 2;
	}
      }

      /*for(size_t i=0; i < verts_to_send_size(0); i++){
	bool used_proc[nprocs];
	for(int x = 0; x < nprocs; x++) used_proc[x] = false;
        for(offset_t j = offsets[verts_to_send(i)]; j < offsets[verts_to_send(i)+1]; j++){
	  lno_t nbor = adjs[j];
	  if(nbor >= nVtx){
	    //one for the GID, one for the color, two overall.
	    if(owners[nbor-nVtx] != comm->getRank() && owners[nbor-nVtx] != -1){
	      if(!used_proc[owners[nbor-nVtx]]){
		//std::cout<<comm->getRank()<<": sending vertex "<<mapOwnedPlusGhosts->getGlobalElement(verts_to_send[i])<<"to proc "<<owners[nbor-nVtx]<<"\n";
	        sendcnts[owners[nbor-nVtx]]+=2;
		used_proc[owners[nbor-nVtx]] = true;
	      }
            }
	  }
	}
      }*/
      int status = MPI_Alltoall(sendcnts,1,MPI_INT,recvcnts,1,MPI_INT,MPI_COMM_WORLD);
      
      int* sdispls = new int[nprocs+1];
      int* rdispls = new int[nprocs+1];
      //construct sdispls (for building sendbuf), and sum the total sendcount
      sdispls[0] = 0;
      rdispls[0] = 0;
      gno_t sendsize = 0;
      gno_t recvsize = 0;
      int* sentcount = new int[nprocs];
      for(int i = 1; i < comm->getSize()+1; i++){
        sdispls[i] = sdispls[i-1] + sendcnts[i-1];
        rdispls[i] = rdispls[i-1] + recvcnts[i-1];
        sendsize += sendcnts[i-1];
        recvsize += recvcnts[i-1];
        sentcount[i-1] = 0;
        std::cout<<comm->getRank()<<": sending "<<sendcnts[i-1]<<" GIDs to proc "<<i-1<<"\n";
        std::cout<<comm->getRank()<<": receiving "<<recvcnts[i-1]<<" GIDs from proc "<<i-1<<"\n";
      }
      recv = recvsize;
      send = sendsize;
      int* sendbuf = new int[sendsize];
      int* recvbuf = new int[recvsize];
      
      for(size_t i = 0; i < exportLIDs.size(); i++){
        if(send_flags(exportLIDs[i])){
	  size_t idx = sdispls[exportPIDs[i]] + sentcount[exportPIDs[i]];
	  sentcount[exportPIDs[i]] += 2;
	  sendbuf[idx++] = mapOwnedPlusGhosts->getGlobalElement(exportLIDs[i]);
	  sendbuf[idx] = colors(exportLIDs[i]);
	}
      }

      /*for(size_t i = 0; i < verts_to_send_size(0); i++){
	bool used_proc[nprocs];
	for(int x = 0; x < nprocs; x++) used_proc[x] = false;
	lno_t curr_vert = verts_to_send(i);
        for(offset_t j = offsets[verts_to_send(i)]; j < offsets[verts_to_send(i)+1]; j++){
	  lno_t nbor = adjs[j];
	  if(nbor >= nVtx){
	    if(owners[nbor-nVtx] != comm->getRank() && owners[nbor-nVtx] != -1){
	      if(!used_proc[owners[nbor-nVtx]]){
	        //grab the last used index
	        int idx = sdispls[owners[nbor-nVtx]] + sentcount[owners[nbor-nVtx]];
	        sentcount[owners[nbor-nVtx]]+=2;
	        //build up the sendbuf
	        sendbuf[idx++] = mapOwnedPlusGhosts->getGlobalElement(curr_vert);
	        sendbuf[idx] = colors(curr_vert);
		used_proc[owners[nbor-nVtx]] = true;
	      }
	    }
	  }
	}
      }*/
      comm->barrier();
      double comm_total = 0.0;
      double comm_temp = timer();
      status = MPI_Alltoallv(sendbuf, sendcnts, sdispls, MPI_INT, recvbuf,recvcnts,rdispls,MPI_INT,MPI_COMM_WORLD);
      comm_total += timer() - comm_temp;
      

      for(int i = 0; i < comm->getSize(); i++){
        for(int j = rdispls[i]; j < rdispls[i+1]; j+=2){
          lno_t lid = mapOwnedPlusGhosts->getLocalElement(recvbuf[j]);
	  if(lid < nVtx) std::cout<<comm->getRank()<<": received a locally owned vertex, somehow\n";
	  colors(lid) = recvbuf[j+1];
        }
      }
       
      return comm_total;
    }
    
    RCP<const base_adapter_t> adapter;
    RCP<GraphModel<base_adapter_t> > model;
    RCP<Teuchos::ParameterList> pl;
    RCP<Environment> env;
    RCP<const Teuchos::Comm<int> > comm;
    int numColors;
    Teuchos::ArrayView<const lno_t> exportLIDs;
    Teuchos::ArrayView<const int> exportPIDs; 
  public:
    //constructor for the  hybrid distributed distance-1 algorithm
    AlgHybridGMB(
      const RCP<const base_adapter_t> &adapter_, 
      const RCP<Teuchos::ParameterList> &pl_,
      const RCP<Environment> &env_,
      const RCP<const Teuchos::Comm<int> > &comm_)
    : adapter(adapter_), pl(pl_), env(env_), comm(comm_) {
      std::cout<<comm->getRank()<<": inside coloring constructor\n";
      numColors = 4;
      modelFlag_t flags;
      flags.reset();
      buildModel(flags);
      std::cout<<comm->getRank()<<": done constructing coloring class\n";
    }


    //Main entry point for graph coloring
    void color( const RCP<ColoringSolution<Adapter> > &solution ) {
      std::cout<<comm->getRank()<<": inside coloring function\n";
      int rank = comm->getRank(); 
      //this will color the global graph in a manner similar to Zoltan
      
      //get vertex GIDs in a locally indexed array (stolen from Ice-Sheet 
      //interface)
      std::cout<<comm->getRank()<<": getting owned vtxIDs\n";
      ArrayView<const gno_t> vtxIDs;
      ArrayView<StridedData<lno_t, scalar_t> > vwgts;
      size_t nVtx = model->getVertexList(vtxIDs, vwgts);
      //we do not use weights at this point
      std::cout<<comm->getRank()<<": getting edge list\n";
      //get edge information from the model
      ArrayView<const gno_t> adjs;
    //  ArrayView<const typename Adapter::offset_t> offsets_inter;
      ArrayView<const offset_t> offsets;
      ArrayView<StridedData<lno_t, scalar_t> > ewgts;
      size_t nEdge = model->getEdgeList(adjs, offsets, ewgts);
      //again, weights are not used
      
      RCP<const map_t> mapOwned;
      RCP<const map_t> mapWithCopies;
      
      std::vector<gno_t> finalGIDs;
      std::vector<offset_t> finalOffset_vec;
      std::vector<lno_t> finalAdjs_vec;      

      lno_t nInterior = 0;
      std::vector<lno_t> reorderToLocal;
      for(size_t i = 0;  i< nVtx; i++) reorderToLocal.push_back(i);
      std::cout<<comm->getRank()<<": Setting up local datastructures\n";
      //Set up a typical local mapping here.
      std::unordered_map<gno_t,lno_t> globalToLocal;
      std::vector<gno_t> ownedPlusGhosts;
      for (gno_t i = 0; i < vtxIDs.size(); i++){
        if(vtxIDs[i] < 0) std::cout<<comm->getRank()<<": found a negative GID\n";
        globalToLocal[vtxIDs[i]] = i;
        ownedPlusGhosts.push_back(vtxIDs[i]);
      }
      gno_t nghosts = 0;
      for (int i = 0; i < adjs.size(); i++){
        if(globalToLocal.count(adjs[i]) == 0){
          //new unique ghost found
          if(adjs[i] < 0) std::cout<<comm->getRank()<<": found a negative adjacency\n";
          ownedPlusGhosts.push_back(adjs[i]);
          globalToLocal[adjs[i]] = vtxIDs.size() + nghosts;
          nghosts++;
            
        }
      }
      std::cout<<comm->getRank()<<": vec.max_size() = "<<finalAdjs_vec.max_size()<<", adjs.size() = "<<adjs.size()<<"\n";  
      finalAdjs_vec.resize(adjs.size()); 
      for(size_t i = 0; i < finalAdjs_vec.size();i++){
        finalAdjs_vec[i] = globalToLocal[adjs[i]];
      }
      for(int i = 0; i < offsets.size(); i++) finalOffset_vec.push_back(offsets[i]);
      finalGIDs = ownedPlusGhosts;
      

      std::cout<<comm->getRank()<<": creating Tpetra Maps\n";
      Tpetra::global_size_t dummy = Teuchos::OrdinalTraits
                                             <Tpetra::global_size_t>::invalid();
      mapOwned = rcp(new map_t(dummy, vtxIDs, 0, comm));

      dummy = Teuchos::OrdinalTraits <Tpetra::global_size_t>::invalid();
      mapWithCopies = rcp(new map_t(dummy, 
                                Teuchos::arrayViewFromVector(ownedPlusGhosts),
                                0, comm)); 
                                      
      //create the FEMultiVector for the distributed communication.
      //We also use the views from this datastructure as arguments to
      //KokkosKernels coloring functions.
      std::cout<<comm->getRank()<<": creating FEMultiVector\n";
      typedef Tpetra::Import<lno_t, gno_t> import_t;
      Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned, 
                                                            mapWithCopies));
      Teuchos::RCP<femv_t> femv = rcp(new femv_t(mapOwned, 
                                                    importer, 1, true));
      //Get color array to fill
      ArrayRCP<int> colors = solution->getColorsRCP();
      for(size_t i=0; i<nVtx; i++){
        colors[i] = 0;
      } 
      
      //Create random numbers seeded on global IDs so that we don't
      //need to communicate for consistency. These numbers determine
      //which vertex gets recolored in the event of a conflict.
      //taken directly from the Zoltan coloring implementation 
      std::vector<int> rand(finalGIDs.size());
      for(size_t i = 0; i < finalGIDs.size(); i++){
        std::srand(finalGIDs[i]);
        rand[i] = std::rand();
      }

      std::vector<int> ghostOwners(finalGIDs.size() - nVtx);
      std::vector<gno_t> ghosts(finalGIDs.size() - nVtx);
      for(int i = nVtx; i < finalGIDs.size(); i++) ghosts[i-nVtx] = finalGIDs[i];
      ArrayView<int> owners = Teuchos::arrayViewFromVector(ghostOwners);
      ArrayView<const gno_t> ghostGIDs = Teuchos::arrayViewFromVector(ghosts);
      Tpetra::LookupStatus ls = mapOwned->getRemoteIndexList(ghostGIDs, owners);
      
      for(int i = 0; i < finalOffset_vec.size()-1; i++){
        std::sort(finalAdjs_vec.begin()+finalOffset_vec[i],finalAdjs_vec.begin()+finalOffset_vec[i+1]);
      }
      
      ArrayView<const offset_t> finalOffsets = Teuchos::arrayViewFromVector(finalOffset_vec);
      ArrayView<const lno_t> finalAdjs = Teuchos::arrayViewFromVector(finalAdjs_vec);
      //sort adjacencies, for KokkosKernels comparison
      /*for(int i = 0; i < 100; i++){
        comm->barrier();
        femv->switchActiveMultiVector();
        double comm_tmp = timer();
        femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
        double comm_time = timer() - comm_tmp;
        std::cout<<comm->getRank()<<": Comm "<<i<<" time: "<<comm_time<<"\n";
      }*/
      
      exportLIDs = importer->getExportLIDs();
      exportPIDs = importer->getExportPIDs();

      for(int i = 0; i < comm->getSize(); i++){
	if(comm->getRank() == i){
	  for(int j = 0; j < exportLIDs.size(); j++){
	    std::cout<<i<<": send lid "<<exportLIDs[j]<<" to proc "<<exportPIDs[j]<<"\n";
	  }
	}
        comm->barrier();
      }

      /*for(int i = 0; i < comm->getSize(); i++){
	if(comm->getRank() == i){
	  for(int j = 0; j < exportPIDs.size(); j++){
	    std::cout<<i<<": "<<exportPIDs[j]<<"\n";
	  }
	}
        comm->barrier();
      }*/
      // call coloring function
      hybridGMB(nVtx, nInterior, finalAdjs, finalOffsets,colors,femv,finalGIDs,rand,owners,mapWithCopies);
      
      //copy colors to the output array.
      for(int i = 0; i < colors.size(); i++){
        colors[reorderToLocal[i]] = femv->getData(0)[i];
      }
      /*for(lno_t i = 0; i < nVtx; i++){
        std::cout<<comm->getRank()<<": global vert "<< finalGIDs[i] <<" is color "<< colors[reorderToLocal[i]]<<"\n";
      }*/
     
      comm->barrier();
    }
     
    void hybridGMB(const size_t nVtx,lno_t nInterior, Teuchos::ArrayView<const lno_t> adjs, 
                   Teuchos::ArrayView<const offset_t> offsets, 
                   Teuchos::ArrayRCP<int> colors, Teuchos::RCP<femv_t> femv,
                   std::vector<gno_t> reorderGIDs,
                   std::vector<int> rand,
                   ArrayView<int> owners,
                   RCP<const map_t> mapOwnedPlusGhosts){
      std::cout<<comm->getRank()<<": inside coloring algorithm\n";
      double total_time = 0.0;
      double interior_time = 0.0;
      double comm_time = 0.0;
      double comp_time = 0.0;
      double recoloring_time = 0.0;
      double conflict_detection = 0.0;
      //for(int i = nVtx; i < reorderGIDs.size(); i++){
        //ghosts[i-nVtx] = reorderGIDs[i];
        //std::cout<<comm->getRank()<<": ghosts["<<i-nVtx<<"] = "<<reorderGIDs[i]<<", Owned by proc"<<owners[i-nVtx]<<"\n";
      //}
      //make views out of arrayViews
      offset_t local_max_degree = 0;
      offset_t global_max_degree = 0;
      for(int i = 0; i < nVtx; i++){
        offset_t curr_degree = offsets[i+1] - offsets[i];
        if(curr_degree > local_max_degree){
          local_max_degree = curr_degree;
        }
      }
      Teuchos::reduceAll<int, offset_t>(*comm,Teuchos::REDUCE_MAX,1, &local_max_degree, &global_max_degree);
      if(comm->getRank() == 0) std::cout<<"Input has max degree "<<global_max_degree<<"\n";
      std::cout<<comm->getRank()<<": creating Kokkos Views\n"; 
      /*Kokkos::View<offset_t*, Tpetra::Map<>::device_type> host_offsets("Host Offset view", offsets.size());
      for(int i = 0; i < offsets.size(); i++){
        host_offsets(i) = offsets[i];
      }
      Kokkos::View<lno_t*, Tpetra::Map<>::device_type> host_adjs("Host Adjacencies view", adjs.size());
      for(int i = 0; i < adjs.size(); i++){
        host_adjs(i) = adjs[i];
      }
      std::cout<<comm->getRank()<<": creating recoloring graph\n";*/
     
      Kokkos::View<offset_t*, device_type> dist_degrees("Owned+Ghost degree view",rand.size());
      typename Kokkos::View<offset_t*, device_type>::HostMirror dist_degrees_host = Kokkos::create_mirror(dist_degrees);
      //set degree counts for ghosts
      for(int i = 0; i < adjs.size(); i++){
        if(adjs[i] < nVtx) continue;
        dist_degrees_host(adjs[i])++;
      }
      //set degree counts for owned verts
      for(int i = 0; i < offsets.size()-1; i++){
        dist_degrees_host(i) = offsets[i+1] - offsets[i];
      }
      
      Kokkos::View<offset_t*, device_type> dist_offsets("Owned+Ghost Offset view", rand.size()+1);
      typename Kokkos::View<offset_t*, device_type>::HostMirror dist_offsets_host = Kokkos::create_mirror(dist_offsets);

      //set offsets and total # of adjacencies
      dist_offsets_host(0) = 0;
      uint64_t total_adjs = 0;
      for(size_t i = 1; i < rand.size()+1; i++){
        dist_offsets_host(i) = dist_degrees_host(i-1) + dist_offsets_host(i-1);
        total_adjs+= dist_degrees_host(i-1);
      }

      Kokkos::View<lno_t*, device_type> dist_adjs("Owned+Ghost adjacency view", total_adjs);
      typename Kokkos::View<lno_t*, device_type>::HostMirror dist_adjs_host = Kokkos::create_mirror(dist_adjs);
      //now, use the degree view as a counter
      for(size_t i = 0; i < rand.size(); i++){
        dist_degrees_host(i) = 0;
      }
      for(int i = 0; i < adjs.size(); i++) dist_adjs_host(i) = adjs[i];
      if(comm->getSize() > 1){
        for(size_t i = 0; i < nVtx; i++){
          for(size_t j = offsets[i]; j < offsets[i+1]; j++){
            //if the adjacency is a ghost
            if( (size_t)adjs[j] >= nVtx){
              //add the symmetric edge to its adjacency list (already accounted for by offsets)
              dist_adjs_host(dist_offsets_host(adjs[j]) + dist_degrees_host(adjs[j])) = i;
              dist_degrees_host(adjs[j])++;
            }
          }
      	}
      }
      
      /*std::cout<<comm->getRank()<<": writing graph adjacency list out to file\n";
      std::ofstream adj_file("graph.adj");
      for(int i = 0; i < rand.size(); i++){
        for(int j = dist_offsets_host(i); j < dist_offsets_host(i+1); j++){
          adj_file<<dist_adjs_host(j)<<" ";
        }
        adj_file<<"\n";
      }
      adj_file.close();
      std::cout<<comm->getRank()<<": Done writing to file\n";*/
      
      std::cout<<comm->getRank()<<": copying host mirrors to device views\n";
      Kokkos::deep_copy(dist_degrees, dist_degrees_host);
      Kokkos::deep_copy(dist_offsets, dist_offsets_host);
      Kokkos::deep_copy(dist_adjs, dist_adjs_host);
      std::cout<<comm->getRank()<<": done copying to device\n";
      
      std::string kokkos_only_interior = pl->get<std::string>("Kokkos_only_interior","false");
      size_t kokkosVerts = nVtx;
      if(kokkos_only_interior == "true" && comm->getSize() != 1) {
        kokkosVerts = nInterior;
      }
      std::cout<<comm->getRank()<<": done creating views, coloring interior\n";
      interior_time = timer();
      //call the KokkosKernels coloring function with the Tpetra default spaces.
      /*this->colorInterior<Tpetra::Map<>::execution_space,
                          Tpetra::Map<>::memory_space,
                          Tpetra::Map<>::memory_space>*/
      bool use_vbbit = (global_max_degree < 6000);
      this->colorInterior<execution_space, memory_space,memory_space>
                 (kokkosVerts, dist_adjs, dist_offsets, colors, femv,dist_adjs,0,use_vbbit);
      interior_time = timer() - interior_time;
      total_time = interior_time;
      comp_time = interior_time;
      //This is the Kokkos version of two queues. These will attempt to be used in parallel.
      /*Kokkos::View<lno_t*, device_type> recoloringQueue("recoloringQueue",nVtx);
      Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA(const int& i){
        recoloringQueue(i) = -1;
      });
      Kokkos::View<lno_t*, device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > recoloringQueue_atomic=recoloringQueue;*/


      Kokkos::View<gno_t[1], device_type> recoloringSize("Recoloring Queue Size");
      recoloringSize(0) = 0;
      Kokkos::View<gno_t[1], device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > recoloringSize_atomic = recoloringSize; 
      Kokkos::View<int*,device_type> rand_dev("randVec",rand.size());
      typename Kokkos::View<int*, device_type>::HostMirror rand_host = Kokkos::create_mirror(rand_dev);
      for(size_t i = 0; i < rand.size(); i++){
        rand_host(i) = rand[i];
      }
      Kokkos::View<gno_t*, device_type> gid_dev("GIDs",reorderGIDs.size());
      typename Kokkos::View<gno_t*,device_type>::HostMirror gid_host = Kokkos::create_mirror(gid_dev);
      for(size_t i = 0; i < reorderGIDs.size(); i++){
        gid_host(i) = reorderGIDs[i];
      }
      Kokkos::deep_copy(rand_dev,rand_host);
      Kokkos::deep_copy(gid_dev, gid_host);
      std::cout<<comm->getRank()<<": done creating recoloring datastructures, begin initial recoloring\n";
      offset_t boundary_size = 0;
      for(offset_t i = 0; i < nVtx; i++){
        for(offset_t j = offsets[i]; j < offsets[i+1]; j++){
          if(adjs[j] >= nVtx) {
            boundary_size++;
            break;
          }
        }
      }
      std::cout<<comm->getRank()<<": creating send views\n";
      Kokkos::View<bool*, device_type> send_flags("send flags",nVtx);
      Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA(const int& i){
	send_flags(i) = false;	      
      });
      Kokkos::View<lno_t*, device_type> verts_to_send_view("verts to send",boundary_size);
      Kokkos::parallel_for(boundary_size, KOKKOS_LAMBDA(const int& i){
        verts_to_send_view(i) = -1;
      });
      Kokkos::View<lno_t*, device_type, Kokkos::MemoryTraits<Kokkos::Atomic>> verts_to_send_atomic = verts_to_send_view;
      
      Kokkos::View<int[1], device_type> verts_to_send_size("verts to send size");
      verts_to_send_size(0) = 0;
      Kokkos::View<int[1], device_type, Kokkos::MemoryTraits<Kokkos::Atomic> > verts_to_send_size_atomic = verts_to_send_size;
      std::cout<<comm->getRank()<<": Done creating send views, initializing...\n";
      std::cout<<comm->getRank()<<": boundary_size = "<<boundary_size<<" verts_to_send_size_atomic(0) = "<<verts_to_send_size_atomic(0)<<"\n";
      Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA(const int&i){
        for(offset_t j = dist_offsets(i); j < dist_offsets(i+1); j++){
	  if(dist_adjs(j) >= nVtx){
	    verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
	    send_flags(i) = true;
	    break;
	  }
	} 
      });
      Kokkos::fence();
      offset_t max_uncolored_degree = 0;
      Kokkos::View<int*, device_type> ghost_colors("ghost color backups", rand.size()-nVtx);
      //typename Kokkos::View<int*, device_type>::HostMirror ghost_colors_host = Kokkos::create_mirror(ghost_colors);
      std::cout<<comm->getRank()<<": Done initializing\n";
      //bootstrap distributed coloring, add conflicting vertices to the recoloring queue.
      gno_t sentPerRound[100];
      gno_t recvPerRound[100];
      if(comm->getSize() > 1){
        comm->barrier();
        Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
        Kokkos::View<int*, device_type> femv_colors = subview(femvColors, Kokkos::ALL, 0);
        //femv->switchActiveMultiVector();
        //double comm_temp = timer();
        for(int i = 0; i < verts_to_send_size(0); i++){
          if(verts_to_send_view(i) < 0 || verts_to_send_view(i) >= nVtx) std::cout<<comm->getRank()<<": found a vert_to_send that is "<<verts_to_send_view(i)<<"\n";
        }
        std::cout<<comm->getRank()<<": going to communicate\n";
        gno_t recv, sent;
        comm_time = doOwnedToGhosts(mapOwnedPlusGhosts,nVtx,/*offsets,adjs,verts_to_send_view,verts_to_send_size, owners*/send_flags,femv_colors,recv,sent);
        sentPerRound[0] = sent;
        recvPerRound[0] = recv;
        std::cout<<comm->getRank()<<": done communicating\n";
	verts_to_send_size(0) = 0;
        //femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
        //comm_time = timer() - comm_temp;
        total_time += comm_time;
        //femv->switchActiveMultiVector();
        //get a subview of the colors:
        
        Kokkos::parallel_for(rand.size()-nVtx,KOKKOS_LAMBDA(const int& i){
          ghost_colors(i) = femv_colors(i+nVtx);
        });
        //Kokkos::deep_copy(ghost_colors_host, ghost_colors);

        comm->barrier();
        Kokkos::View<int, device_type, Kokkos::MemoryTraits<Kokkos::Atomic>> boolean_checks("check yoself");
        boolean_checks() = 0;
        double temp = timer();
        /*Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA (const int& i){
          int currColor = femv_colors(i);
          for(offset_t j = dist_offsets(i); j < dist_offsets(i+1); j++){
            int nborColor = femv_colors(dist_adjs(j));
            if(currColor == nborColor ){
              if(rand_dev(i) > rand_dev(dist_adjs(j))){
                femv_colors(i) = 0;
                recoloringSize_atomic(0)++;
		verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
                break;
              }else if(rand_dev(dist_adjs(j)) > rand_dev(i)){
                femv_colors(dist_adjs(j)) = 0;
                recoloringSize_atomic(0)++;
              } else {
                if (gid_dev(i) >= gid_dev(dist_adjs(j))){ //&& rand_dev(i) == rand_dev(dist_adjs(j))){
                  femv_colors(i) = 0;
                  recoloringSize_atomic(0)++;
		  verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
                  break;
                } else {
                  femv_colors(dist_adjs(j)) = 0;
                  recoloringSize_atomic(0)++;
                }
              }
            }
          }
        });*/
        Kokkos::parallel_reduce(rand.size()-nVtx, KOKKOS_LAMBDA(const int& i,gno_t& recoloring_size){
          lno_t lid = i+nVtx;
          int currColor = femv_colors(lid);
          for(offset_t j = dist_offsets(lid); j < dist_offsets(lid+1); j++){
            int nborColor = femv_colors(dist_adjs(j));
            //boolean_checks()++;
            if(currColor == nborColor ){
              //boolean_checks()++;
              if(rand_dev(lid) > rand_dev(dist_adjs(j))){
                //boolean_checks()+=2;
                femv_colors(lid) = 0;
                recoloring_size++;
                //recoloringSize_atomic(0)++;
                //verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
                break;
              }else if (rand_dev(dist_adjs(j)) > rand_dev(lid)){
                //boolean_checks()+=2;
                femv_colors(dist_adjs(j)) = 0;
                recoloring_size++;
                //recoloringSize_atomic(0)++;
              } else {
                //boolean_checks()++;
                if (gid_dev(lid) >= gid_dev(dist_adjs(j))){// && rand_dev(i) == rand_dev(dist_adjs(j))){
                  //boolean_checks()+=2;
                  femv_colors(lid) = 0;
                  recoloring_size++;
                  //recoloringSize_atomic(0)++;
                  //verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
                  break;
                } else {
                  //boolean_checks()+=2;
                  femv_colors(dist_adjs(j)) = 0;
                  recoloring_size++;
                  //recoloringSize_atomic(0)++;
                }
              }
            }
          }
        },recoloringSize(0));
        Kokkos::fence();
        Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA(const int& i){
          if(femv_colors(i) == 0){
            //boolean_checks()++;
            verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
	    send_flags(i) = true;
          } else {
	    send_flags(i) = false;
	  }
        });
        //ensure the parallel_for has completed before continuing.
        Kokkos::fence();
        std::cout<<comm->getRank()<<": boolean_checks = "<<boolean_checks()<<"\n";
        conflict_detection += timer() - temp;
        total_time += conflict_detection;
        comp_time += conflict_detection;
      }
      std::cout<<comm->getRank()<<": done initial recoloring, begin recoloring loop\n";
      double totalPerRound[100];
      double commPerRound[100];
      double compPerRound[100];
      double recoloringPerRound[100];
      double conflictDetectionPerRound[100];
      double serialRecoloringPerRound[100];
      int vertsPerRound[100];
      bool done = false; //We're only done when all processors are done
      if(comm->getSize() == 1) done = true;
      totalPerRound[0] = interior_time + comm_time + conflict_detection;
      recoloringPerRound[0] = 0;
      commPerRound[0] = comm_time;
      compPerRound[0] = interior_time + conflict_detection;
      conflictDetectionPerRound[0] = conflict_detection;
      recoloringPerRound[0] = 0;
      vertsPerRound[0] = 0;
      //sentPerRound[0] = 0;
      //recvPerRound[0] = 0;
      int distributedRounds = 1; //this is the same across all processors
      //typename Kokkos::View<offset_t*, device_type>::HostMirror verts_to_send_host = Kokkos::create_mirror(verts_to_send_view);
      //bool serial = false;
      //while the queue is not empty
      while(recoloringSize(0) > 0 || !done){
        //get a subview of the colors:
        Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
        Kokkos::View<int*, device_type> femv_colors = subview(femvColors, Kokkos::ALL, 0);
        //color everything in the recoloring queue, put everything on conflict queue
        if(distributedRounds < 100) {
          int localVertsToRecolor = 0;
          for(int i = 0; i < nVtx; i ++){
            if(femv_colors(i) == 0) localVertsToRecolor++;
          }
          vertsPerRound[distributedRounds] = localVertsToRecolor;//recoloringSize(0);
        }
        std::cout<<comm->getRank()<<": starting to recolor\n";
        for(int i = 0; i < nVtx; i++){
          offset_t degree = offsets[i+1] - offsets[i];
          if(femv_colors(i) == 0 && degree > max_uncolored_degree) max_uncolored_degree = degree;
        }
        std::cout<<comm->getRank()<<": max_uncolored_degree = "<<max_uncolored_degree<<"\n";
	Kokkos::View<lno_t*, device_type> verts_to_recolor("verts_to_recolor", boundary_size);
        //typename Kokkos::View<offset_t*, device_type>::HostMirror verts_to_recolor_host = Kokkos::create_mirror(verts_to_recolor);        
	Kokkos::deep_copy(verts_to_recolor, verts_to_send_view);
        /*if(!serial && verts_to_send_size(0) <= 5000){
          serial = true;
        }
        if(!serial){
          Kokkos::deep_copy(verts_to_recolor_host, verts_to_send_view);
          Kokkos::deep_copy(verts_to_send_host, verts_to_send_view);
        }*/
        //Kokkos::deep_copy(verts_to_recolor_host, verts_to_send_view);
        //std::vector<int> old_colors;
        /*for(int i = 0; i < femv_colors.size(); i++){
          old_colors.push_back(femv_colors(i));
        }


        comm->barrier();
        double recolor_serial = timer();
        if(verts_to_send_size(0) > 0){
          this->colorInterior<Kokkos::DefaultHostExecutionSpace,
                              Kokkos::CudaUVMSpace::memory_space,
                              Kokkos::CudaUVMSpace::memory_space>
                             (femv_colors.size(), dist_adjs_host, dist_offsets_host, colors, femv, verts_to_recolor_host, verts_to_send_size(0),use_vbbit);
        }
        serialRecoloringPerRound[distributedRounds] = timer() - recolor_serial;

        for(int i = 0; i < femv_colors.size(); i++){
          femv_colors(i) = old_colors[i];
        }*/

        comm->barrier();
        double recolor_temp = timer();
        //use KokkosKernels to recolor the conflicting vertices.  
	if(verts_to_send_size(0) > 0){
          //if(/*verts_to_send_size(0) > 5000*/true){
            this->colorInterior<execution_space,
                                memory_space,
                                memory_space>
                               (femv_colors.size(),dist_adjs,dist_offsets,colors,femv,verts_to_recolor,verts_to_send_size(0),use_vbbit);
          /*} else {
            this->colorInterior<Kokkos::DefaultHostExecutionSpace,
                                Kokkos::CudaUVMSpace::memory_space,
                                Kokkos::CudaUVMSpace::memory_space>
                               (femv_colors.size(), dist_adjs_host, dist_offsets_host, colors, femv, verts_to_recolor_host, verts_to_send_size(0),use_vbbit);    
          }*/
	}
        recoloringPerRound[distributedRounds] = timer() - recolor_temp;
        recoloring_time += recoloringPerRound[distributedRounds];
        total_time += recoloringPerRound[distributedRounds];
        comp_time += recoloringPerRound[distributedRounds];
        compPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
        totalPerRound[distributedRounds] = recoloringPerRound[distributedRounds];
	
        
        std::cout<<comm->getRank()<<": done recoloring\n";
        recoloringSize(0) = 0;

        //replace old ghost colors
        //if(/*verts_to_send_size(0) > 5000*/true){
          Kokkos::parallel_for(rand.size()-nVtx, KOKKOS_LAMBDA(const int& i){
            femv_colors(i+nVtx) = ghost_colors(i);
          });
          Kokkos::fence();
        /*} else {
          for(int i = 0; i < ghost_colors_host.size(); i++){
            femv_colors(i+nVtx) = ghost_colors_host(i);
          }
        }*/
        //communicate
        typename Kokkos::View<lno_t*, device_type>::HostMirror verts_to_send_host = Kokkos::create_mirror(verts_to_send_view);
        Kokkos::deep_copy(verts_to_send_host, verts_to_send_view);
        comm->barrier();
        //femv->switchActiveMultiVector();
        double comm_temp = timer();
        //femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);
	//comm_time = doGhostUpdate(mapOwnedPlusGhosts,nVtx,inst_to_req,owners,femv_colors);
	gno_t sent,recv;
        commPerRound[distributedRounds] = doOwnedToGhosts(mapOwnedPlusGhosts,nVtx,/*offsets,adjs,verts_to_send_host, verts_to_send_size, owners*/send_flags,femv_colors,recv,sent); 	
        std::cout<<comm->getRank()<<": total sent in round "<<distributedRounds<<" = "<<sent<<"\n";
        std::cout<<comm->getRank()<<": total recv in round "<<distributedRounds<<" = "<<recv<<"\n";
        sentPerRound[distributedRounds] = sent;
        recvPerRound[distributedRounds] = recv;
        //commPerRound[distributedRounds] = doOwnedToGhosts(mapOwnedPlusGhosts,nVtx,offsets,adjs,verts_to_send,owners,femv_colors); 
	//commPerRound[distributedRounds] += doGhostUpdate(mapOwnedPlusGhosts,nVtx,verts_to_req,owners,femv_colors);
        //commPerRound[distributedRounds] = timer() - comm_temp;
        comm_time += commPerRound[distributedRounds];
        totalPerRound[distributedRounds] += commPerRound[distributedRounds];
        total_time += commPerRound[distributedRounds];
        //femv->switchActiveMultiVector();
        //detect conflicts in parallel. For a detected conflict,
        //reset the vertex-to-be-recolored's color to 0, in order to
        //allow KokkosKernels to recolor correctly.
        //if(/*verts_to_send_size(0) > 5000*/true){
          Kokkos::parallel_for(rand.size()-nVtx, KOKKOS_LAMBDA(const int& i){
            ghost_colors(i) = femv_colors(i+nVtx);
          });
          Kokkos::fence();
          //Kokkos::deep_copy(ghost_colors_host, ghost_colors);
        /*} else {
          for(int i = 0; i < ghost_colors_host.size(); i++){
            ghost_colors_host(i) = femv_colors(i+nVtx);
          }
        }*/
        //int prev_verts_to_send = verts_to_send_size(0);
	verts_to_send_size(0) = 0;
        comm->barrier();
        double detection_temp = timer();
        /*Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA (const int& i){
          int currColor = femv_colors(i);
          for(offset_t j = dist_offsets(i); j < dist_offsets(i+1); j++){
            int nborColor = femv_colors(dist_adjs(j));
            if(currColor == nborColor ){
              if(rand_dev(i) > rand_dev(dist_adjs(j))){
                femv_colors(i) = 0;
                recoloringSize_atomic(0)++;
		verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
                break;
              }else if (rand_dev(dist_adjs(j)) > rand_dev(i)){
                femv_colors(dist_adjs(j)) = 0;
                recoloringSize_atomic(0)++;
              } else {
                if (gid_dev(i) >= gid_dev(dist_adjs(j))){// && rand_dev(i) == rand_dev(dist_adjs(j))){
                  femv_colors(i) = 0;
                  recoloringSize_atomic(0)++;
		  verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
                  break;
                } else {
                  femv_colors(dist_adjs(j)) = 0;
                  recoloringSize_atomic(0)++;
                }
              }
            }
          }
        });*/
        //if verts to recolor > 5000
        //if(/*prev_verts_to_send > 5000*/true){
        /*  Kokkos::parallel_for(rand.size()-nVtx, KOKKOS_LAMBDA(const int& i){
            lno_t lid = i+nVtx;
            int currColor = femv_colors(lid);
            for(offset_t j = dist_offsets(lid); j < dist_offsets(lid+1); j++){
              int nborColor = femv_colors(dist_adjs(j));
              if(currColor == nborColor ){
                if(rand_dev(lid) > rand_dev(dist_adjs(j))){
                  femv_colors(lid) = 0;
                  recoloringSize_atomic(0)++;
                  //verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
                  break;
                }else if (rand_dev(dist_adjs(j)) > rand_dev(lid)){
                  femv_colors(dist_adjs(j)) = 0;
                  recoloringSize_atomic(0)++;
                } else {
                  if (gid_dev(lid) >= gid_dev(dist_adjs(j))){// && rand_dev(i) == rand_dev(dist_adjs(j))){
                    femv_colors(lid) = 0;
                    recoloringSize_atomic(0)++;
                    //verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
                    break;
                  } else {
                    femv_colors(dist_adjs(j)) = 0;
                    recoloringSize_atomic(0)++;
                  }
                }
              }
            }
          });
          Kokkos::fence();*/
        Kokkos::parallel_reduce(rand.size()-nVtx, KOKKOS_LAMBDA(const int& i,gno_t& recoloring_size){
          lno_t lid = i+nVtx;
          int currColor = femv_colors(lid);
          for(offset_t j = dist_offsets(lid); j < dist_offsets(lid+1); j++){
            int nborColor = femv_colors(dist_adjs(j));
            //boolean_checks()++;
            if(currColor == nborColor ){
              //boolean_checks()++;
              if(rand_dev(lid) > rand_dev(dist_adjs(j))){
                //boolean_checks()+=2;
                femv_colors(lid) = 0;
                recoloring_size++;
                //recoloringSize_atomic(0)++;
                //verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
                break;
              }else if (rand_dev(dist_adjs(j)) > rand_dev(lid)){
                //boolean_checks()+=2;
                femv_colors(dist_adjs(j)) = 0;
                recoloring_size++;
                //recoloringSize_atomic(0)++;
              } else {
                //boolean_checks()++;
                if (gid_dev(lid) >= gid_dev(dist_adjs(j))){// && rand_dev(i) == rand_dev(dist_adjs(j))){
                  //boolean_checks()+=2;
                  femv_colors(lid) = 0;
                  recoloring_size++;
                  //recoloringSize_atomic(0)++;
                  //verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
                  break;
                } else {
                  //boolean_checks()+=2;
                  femv_colors(dist_adjs(j)) = 0;
                  recoloring_size++;
                  //recoloringSize_atomic(0)++;
                }
              }
            }
          }
        },recoloringSize(0));
          Kokkos::parallel_for(nVtx, KOKKOS_LAMBDA(const int& i){
            if(femv_colors(i) == 0){
              verts_to_send_atomic(verts_to_send_size_atomic(0)++) = i;
            }
          });
        
          //For Cuda, this fence is necessary to ensure the Kokkos::parallel_for is finished
          //before continuing with the coloring. 
          Kokkos::fence();
        /*} else {
          for(int i = nVtx; i < rand.size(); i++){
            lno_t lid = i;
            int currColor = femv_colors(lid);
            for(offset_t j = dist_offsets_host(lid); j < dist_offsets_host(lid+1); j++){
              int nborColor = femv_colors(dist_adjs_host(j));
              if(currColor == nborColor){
                if(rand[lid] > rand[dist_adjs_host(j)]){
                  femv_colors(lid) = 0;
                  recoloringSize_atomic(0)++;
                  verts_to_send_host(verts_to_send_size(0)++) = lid;
                  break;
                } else if (rand[dist_adjs_host(j)] > rand[lid]){
                  femv_colors(dist_adjs_host(j)) = 0;
                  recoloringSize_atomic(0)++;
                } else {
                  if(gid_host(lid) >= gid_host(dist_adjs_host(j))){
                    femv_colors(lid) = 0;
                    recoloringSize_atomic(0)++;
                    verts_to_send_host(verts_to_send_size(0)++) = lid;
                    break;
                  } else {
                    femv_colors(dist_adjs_host(j)) = 0;
                    recoloringSize_atomic(0)++;
                  }
                }
              }
            }
          }
        }*/
        conflictDetectionPerRound[distributedRounds] = timer() - detection_temp;
        conflict_detection += conflictDetectionPerRound[distributedRounds];
        compPerRound[distributedRounds] += conflictDetectionPerRound[distributedRounds];
        totalPerRound[distributedRounds] += conflictDetectionPerRound[distributedRounds];
        total_time += conflictDetectionPerRound[distributedRounds];
        comp_time += conflictDetectionPerRound[distributedRounds];
        
        //do a reduction to determine if we're done
        int globalDone = 0;
        int localDone = recoloringSize(0);
        Teuchos::reduceAll<int, int>(*comm,Teuchos::REDUCE_SUM,1, &localDone, &globalDone);
        //std::cout<<comm->getRank()<<": globaldone="<<globalDone<<"\n";
        //We're only allowed to stop once everyone has no work to do.
        //collectives will hang if one process exits. 
        distributedRounds++;
        done = !globalDone;
      }
      
      std::cout<<comm->getRank()<<": done recoloring loop, computing statistics\n";
      int localBoundaryVertices = 0;
      for(int i = 0; i < nVtx; i++){
        for(int j = offsets[i]; j < offsets[i+1]; j++){
          if(adjs[j] >= nVtx){
            localBoundaryVertices++;
            break;
          }
        }
      }
      Kokkos::View<int**, Kokkos::LayoutLeft> femvColors = femv->template getLocalView<memory_space>();
      Kokkos::View<int*, device_type> femv_colors = subview(femvColors, Kokkos::ALL, 0);
      //print how many rounds of speculating/correcting happened (this should be the same for all ranks):
      if(comm->getRank()==0) printf("did %d rounds of distributed coloring\n", distributedRounds);
      int totalVertsPerRound[100];
      int totalBoundarySize = 0;
      double finalTotalPerRound[100];
      double maxRecoloringPerRound[100];
      double finalSerialRecoloringPerRound[100];
      double minRecoloringPerRound[100];
      double finalCommPerRound[100];
      double finalCompPerRound[100];
      double finalConflictDetectionPerRound[100];
      gno_t finalRecvPerRound[100];
      gno_t finalSentPerRound[100];
      for(int i = 0; i < 100; i++) {
        totalVertsPerRound[i] = 0;
        finalTotalPerRound[i] = 0.0;
        maxRecoloringPerRound[i] = 0.0;
        minRecoloringPerRound[i] = 0.0;
        finalCommPerRound[i] = 0.0;
        finalCompPerRound[i] = 0.0;
        finalConflictDetectionPerRound[i] = 0.0;
        finalSentPerRound[i] = 0;
        finalRecvPerRound[i] = 0;
      }
      Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM,1, &localBoundaryVertices,&totalBoundarySize);
      Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM,100,vertsPerRound,totalVertsPerRound);
      Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,totalPerRound,finalTotalPerRound);
      Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,recoloringPerRound,maxRecoloringPerRound);
      Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MIN,100,recoloringPerRound,minRecoloringPerRound);
      Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,serialRecoloringPerRound,finalSerialRecoloringPerRound);
      Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,commPerRound,finalCommPerRound);
      Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,compPerRound,finalCompPerRound);
      Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,100,conflictDetectionPerRound, finalConflictDetectionPerRound);
      Teuchos::reduceAll<int,gno_t> (*comm, Teuchos::REDUCE_SUM,100,recvPerRound, finalRecvPerRound);
      Teuchos::reduceAll<int,gno_t> (*comm, Teuchos::REDUCE_SUM,100,sentPerRound, finalSentPerRound);
      printf("Rank %d: boundary size: %d\n",comm->getRank(),localBoundaryVertices);
      for(int i = 0; i < std::min(distributedRounds,100); i++){
        printf("Rank %d: recolor %d vertices in round %d\n",comm->getRank(),vertsPerRound[i],i);
        if(comm->getRank()==0) printf("recolored %d vertices in round %d\n",totalVertsPerRound[i],i);
        if(comm->getRank()==0) printf("total time in round %d: %f\n",i,finalTotalPerRound[i]);
        if(comm->getRank()==0) printf("recoloring time in round %d: %f\n",i,maxRecoloringPerRound[i]);
        if(comm->getRank()==0) printf("serial recoloring time in round %d: %f\n",i,finalSerialRecoloringPerRound[i]);
        if(comm->getRank()==0) printf("min recoloring time in round %d: %f\n",i,minRecoloringPerRound[i]);
        if(comm->getRank()==0) printf("conflict detection time in round %d: %f\n",i,finalConflictDetectionPerRound[i]);
        if(comm->getRank()==0) printf("comm time in round %d: %f\n",i,finalCommPerRound[i]);
        if(comm->getRank()==0) printf("total sent in round %d: %d\n",i,finalSentPerRound[i]);
        if(comm->getRank()==0) printf("total recv in round %d: %d\n",i,finalRecvPerRound[i]);
        if(comm->getRank()==0) printf("comp time in round %d: %f\n",i,finalCompPerRound[i]);
      }
      double global_total_time = 0.0;
      double global_recoloring_time=0.0;
      double global_min_recoloring_time=0.0;
      double global_conflict_detection=0.0;
      double global_comm_time=0.0;
      double global_comp_time=0.0;
      double global_interior_time = 0.0;
      Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,1,&total_time,&global_total_time);
      Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,1,&recoloring_time,&global_recoloring_time);
      Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MIN,1,&recoloring_time,&global_min_recoloring_time);
      Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,1,&conflict_detection,&global_conflict_detection);
      Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,1,&comm_time,&global_comm_time);
      Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,1,&comp_time,&global_comp_time);
      Teuchos::reduceAll<int,double>(*comm, Teuchos::REDUCE_MAX,1,&interior_time,&global_interior_time);
      comm->barrier();
      fflush(stdout);
      if(comm->getRank()==0){
        printf("Boundary size: %d\n",totalBoundarySize);
        printf("Total Time: %f\n",global_total_time);
        printf("Interior Time: %f\n",global_interior_time);
        printf("Recoloring Time: %f\n",global_recoloring_time);
        printf("Min Recoloring Time: %f\n",global_min_recoloring_time);
        printf("Conflict Detection Time: %f\n",global_conflict_detection);
        printf("Comm Time: %f\n",global_comm_time);
        printf("Comp Time: %f\n",global_comp_time);
      }
      std::cout<<comm->getRank()<<": exiting coloring\n";
    }
};

template <typename Adapter>
void AlgHybridGMB<Adapter>::buildModel(modelFlag_t &flags){
  flags.set(REMOVE_SELF_EDGES);
    
  this->env->debug(DETAILED_STATUS, "   building graph model");
  std::cout<<comm->getRank()<<": starting to construct graph model\n";
  this->model = rcp(new GraphModel<base_adapter_t>(this->adapter, this->env,
                                                   this->comm, flags));
  std::cout<<comm->getRank()<<": done constructing graph model\n";
  this->env->debug(DETAILED_STATUS, "   graph model built");
}


}//end namespace Zoltan2
#endif
