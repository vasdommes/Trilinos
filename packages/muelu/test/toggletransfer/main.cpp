// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <MueLu_Exceptions.hpp>
#include <MueLu_ParameterListInterpreter.hpp>

#include "MueLu_SemiCoarsenPFactory.hpp" // for semi-coarsening constants
#include <MueLu_TestHelpers.hpp>
#include <MueLu_MultiPhys.hpp>

/**********************************************************************************/
/* CREATE INITAL MATRIX                                                           */
/**********************************************************************************/
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#ifdef HAVE_MUELU_TPETRA
#include <BelosTpetraAdapter.hpp>
#endif
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#endif

#include <unistd.h>
/**********************************************************************************/

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  Teuchos::oblackholestream blackhole;

  bool success = true;
  bool verbose = true;
  try {
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    int nx, ny;

    // current problem read from files is on a 25 x 25 mesh
    //
    // The block system looks like
    //
    //     [ E    0   0    0    0] 
    //     [ 0   D1   0   L1    0]
    //     [ 0    0  D2    0   L2]
    //     [ 0   L1   0   D3    0]
    //     [ 0    0  L2    0   D4]
    //
    // E  (1250x1250)  corresponds to a 2D lin elasticity
    // L1 ( 625x 625) is a variable coef Poisson problem
    // L2 ( 625x 625) is a variable coef Poisson problem
    // D1-D4 are 625x625 diagonal matrices
    
    nx = 25; ny = 25;

    Teuchos::RCP<const Map>   mapScalar = MapFactory::Build(lib, nx*ny, 0, comm);
    Teuchos::RCP<const Map>   map2DofsPerNode = Xpetra::MapFactory<LO,GO,Node>::Build(mapScalar, 2);

    // build map for 5x5 block system  and then read in this matrix

    GO  maxGID = map2DofsPerNode->getMaxAllGlobalIndex();
    GO  minGID = map2DofsPerNode->getMinAllGlobalIndex();

    Teuchos::ArrayView<const GO> mapGIDs  = map2DofsPerNode->getLocalElementList();
    Teuchos::Array<GO>  globalsMultiPhysA(mapGIDs.size()*3); 

    for (int ii = 0; ii < mapGIDs.size(); ii++) { 
      globalsMultiPhysA[ii                    ] = mapGIDs[ii]; 
      globalsMultiPhysA[ii +  mapGIDs.size()  ] = mapGIDs[ii] +   (maxGID - minGID + 1);
      globalsMultiPhysA[ii +  mapGIDs.size()*2] = mapGIDs[ii] + 2*(maxGID - minGID + 1);
    }
    Teuchos::RCP<const Map>   mapMultiPhysA = Xpetra::MapFactory<LO,GO,Node>::Build(lib, Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), globalsMultiPhysA, map2DofsPerNode->getIndexBase(), comm);
    Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >     multiPhysA = Xpetra::IO<SC,LO,GO,Node>::Read("multiPhys.mat",mapMultiPhysA); 


    // To solve the block system a block diagonal auxiliary operator is used in conjunction with 5 invocations of Hierarchy Setup
    //
    //     [ aux4    0     0      0      0] 
    //     [ 0    aux1     0      0      0]
    //     [ 0       0  aux3      0      0]
    //     [ 0       0     0   aux2      0]
    //     [ 0       0     0      0   aux1]
    //
    // Here aux4 correspond to 2 distance Laplacians (one corresponding to even dofs and the other to odd dofs).
    // aux1 is a single distance Laplacian. aux2 is L1 and aux3 is L2.
    
    Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >     aux1 = Xpetra::IO<SC,LO,GO,Node>::Read("aux1.mat",mapScalar);  
    Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >     aux2 = Xpetra::IO<SC,LO,GO,Node>::Read("aux2.mat",mapScalar);  
    Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >     aux3 = Xpetra::IO<SC,LO,GO,Node>::Read("aux3.mat",mapScalar); 
    Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >     aux4 = Xpetra::IO<SC,LO,GO,Node>::Read("aux4.mat",map2DofsPerNode); 

   Teuchos::ArrayRCP<Teuchos::RCP<Matrix>> arrayOfAuxMatrices(5);
    arrayOfAuxMatrices[0] = aux4;
    arrayOfAuxMatrices[1] = aux1;
    arrayOfAuxMatrices[2] = aux3;
    arrayOfAuxMatrices[3] = aux2;
    arrayOfAuxMatrices[4] = aux1;

    // We also pass in 5 sets of coordinates and null spaces

    typedef Teuchos::ScalarTraits<SC> STS;
    typedef typename STS::magnitudeType real_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::coordinateType coordinate_type;
    typedef Xpetra::MultiVector<real_type,LO,GO,NO> RealValuedMultiVector;
    typedef Xpetra::MultiVectorFactory<coordinate_type,LO,GO,NO> RealValuedMultiVectorFactory;

    Teuchos::ArrayRCP<Teuchos::RCP<RealValuedMultiVector>> arrayOfCoords(5);
    Teuchos::ArrayRCP<Teuchos::RCP<MultiVector>>           arrayOfNullspaces(5);

    Teuchos::RCP<RealValuedMultiVector> Coords = RealValuedMultiVectorFactory::Build(mapScalar, 2);
    Teuchos::ArrayRCP<real_type> x = Coords->getDataNonConst(0);
    Teuchos::ArrayRCP<real_type> y = Coords->getDataNonConst(1);
    Teuchos::ArrayView<const GO> coordGIDs = mapScalar->getLocalElementList();
    Scalar meanX = 0.0, meanY = 0.0;
    for (GO p = 0; p < coordGIDs.size(); p += 1) { 
       GlobalOrdinal ind = coordGIDs[p];
       size_t i = ind % nx, j = ind / nx;
       x[p] = (i+1);
       y[p] = (j+1);
       meanX += x[p];
       meanY += y[p];
    }
    meanX = meanX/((Scalar) coordGIDs.size());
    meanY = meanY/((Scalar) coordGIDs.size());

    arrayOfCoords[0] = Coords; 
    arrayOfCoords[1] = Coords; 
    arrayOfCoords[2] = Coords; 
    arrayOfCoords[3] = Coords; 
    arrayOfCoords[4] = Coords; 

    Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>> scalarNullSpace = MultiVectorFactory::Build(mapScalar, 1); 
    Scalar one2(1.0);
    scalarNullSpace->putScalar(one2);

    Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>> elasNullSpace = MultiVectorFactory::Build(map2DofsPerNode, 3); 
    Scalar zero2(0.0);

    elasNullSpace->putScalar(zero2);
    Teuchos::ArrayRCP<Scalar> nsValues0, nsValues1, nsValues2;

    nsValues0 = elasNullSpace->getDataNonConst(0);
    nsValues1 = elasNullSpace->getDataNonConst(1);
    nsValues2 = elasNullSpace->getDataNonConst(2);
    for (int j=0; j<2*coordGIDs.size(); j+=2) {
      // translation
      nsValues0[j+0] = 1.0;
      nsValues1[j+1] = 1.0;
      // rotate around z-axis (x-y plane)
      nsValues2[j+0] = -(y[j/2]-meanY);
      nsValues2[j+1] =  (x[j/2]-meanX);
    }

    arrayOfNullspaces[0] = elasNullSpace;
    arrayOfNullspaces[1] = scalarNullSpace;
    arrayOfNullspaces[2] = scalarNullSpace;
    arrayOfNullspaces[3] = scalarNullSpace;
    arrayOfNullspaces[4] = scalarNullSpace;

    Teuchos::ParameterList  bigList;
    Teuchos::updateParametersFromXmlFileAndBroadcast("bigList.xml", Teuchos::Ptr<Teuchos::ParameterList>(&bigList), *comm);
    Teuchos::RCP<Operator> preconditioner = rcp(new MueLu::MultiPhys<SC,LO,GO,NO>(multiPhysA,arrayOfAuxMatrices, arrayOfNullspaces, arrayOfCoords, 5, bigList, true));

    Teuchos::RCP<Teuchos::TimeMonitor> globalTimeMonitor = Teuchos::rcp (new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Timings: Global Time")));

    // =========================================================================
    // System solution (Ax = b)
    // =========================================================================
    comm->barrier();
    typedef Teuchos::ScalarTraits<SC> STS;
    SC zero = STS::zero(), one = STS::one();

    Teuchos::RCP<Vector> X = VectorFactory::Build(multiPhysA->getRowMap());
    Teuchos::RCP<Vector> B = VectorFactory::Build(multiPhysA->getRowMap());

    {
      // set seed for reproducibility
      Utilities::SetRandomSeed(*comm);
      X->randomize();
      multiPhysA->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

      Teuchos::Array<typename STS::magnitudeType> norms(1);
      B->norm2(norms);
      B->scale(one/norms[0]);
      X->putScalar(zero);
    }

    comm->barrier();

    // Belos linear problem
    typedef MultiVector          MV;
    typedef Belos::OperatorT<MV> OP;


    Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(multiPhysA));
    Teuchos::RCP<OP> belosPrecOp = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(preconditioner)); // Turns a Xpetra::Matrix object into a Belos operator
    Teuchos::RCP<Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
    belosProblem->setRightPrec(belosPrecOp);
    bool set = belosProblem->setProblem();
    if (set == false) {
      std::cout << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }
    Teuchos::RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList();
    belosList->set("Maximum Iterations",    100); 
    belosList->set("Convergence Tolerance", 1e-6);
    belosList->set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList->set("Output Frequency",      1);
    belosList->set("Output Style",          Belos::Brief);
    belosList->set("Implicit Residual Scaling", "None");

    Belos::SolverFactory<SC, MV, OP> solverFactory;
    Teuchos::RCP<Belos::SolverManager<SC, MV, OP> > solver = solverFactory.create("cg", belosList);
    solver->setProblem(belosProblem);
    Belos::ReturnType retStatus = Belos::Unconverged;
    retStatus = solver->solve();

    int iters = solver -> getNumIters();
    success = (iters<50 && retStatus == Belos::Converged);
    if (success)
      std::cout << "SUCCESS! Belos converged in " << iters << " iterations." << std::endl;
    else
      std::cout << "FAILURE! Belos did not converge fast enough." << std::endl;

    // Timer final summaries
    globalTimeMonitor = Teuchos::null; // stop this timer before summary
    Teuchos::TimeMonitor::summarize();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}
