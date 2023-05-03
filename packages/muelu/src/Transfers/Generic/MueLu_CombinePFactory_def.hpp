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
#ifndef MUELU_COMBINEPFACTORY_DEF_HPP
#define MUELU_COMBINEPFACTORY_DEF_HPP

#include <stdlib.h>
#include <iomanip>


// #include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include <Xpetra_IO.hpp>

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_CombinePFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> CombinePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->setEntry("combine: numBlks",ParameterEntry(1));
    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void CombinePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  DeclareInput(Level& fineLevel, Level& /* coarseLevel */) const {
//    Input(fineLevel, "subblock");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void CombinePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel,
                                Level& coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void CombinePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel,
                                Level& coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    const ParameterList& pL = GetParameterList();
    const LO nBlks = as<LO>(pL.get<int> ("combine: numBlks"));

    RCP<Matrix>                A             = Get< RCP<Matrix> >               (fineLevel, "A");

    // Record all matrices that each define a block in block diagonal comboP
    // matrix used for PDE/multiblock interpolation.  Additionally, count 
    // total number of local rows, nonzeros, coarseDofs, and colDofs.

    Teuchos::ArrayRCP<RCP<Matrix>>       arrayOfMatrices(nBlks);
    size_t nComboRowMap = 0, nnzCombo = 0,  nComboColMap = 0, nComboDomMap = 0;

coarseLevel.print(std::cout,Teuchos::VERB_EXTREME);
    Teuchos::ArrayRCP<size_t> DomMapSizePerBlk(nBlks);
    Teuchos::ArrayRCP<size_t> ColMapSizePerBlk(nBlks);
    for(int j = 0; j < nBlks; j++) {
      std::string blockName = "Psubblock" + Teuchos::toString(j);
      if (coarseLevel.IsAvailable(blockName, NoFactory::get())) {
        arrayOfMatrices[j] = coarseLevel.Get< RCP<Matrix> >(blockName, NoFactory::get());
        nComboRowMap       += Teuchos::as<size_t>((arrayOfMatrices[j])->getRowMap()->getLocalNumElements());
        DomMapSizePerBlk[j] = Teuchos::as<size_t>((arrayOfMatrices[j])->getDomainMap()->getLocalNumElements());
        ColMapSizePerBlk[j] = Teuchos::as<size_t>((arrayOfMatrices[j])->getColMap()->getLocalNumElements());
        nComboDomMap       += DomMapSizePerBlk[j];
        nComboColMap       += ColMapSizePerBlk[j];
        nnzCombo         += Teuchos::as<size_t>((arrayOfMatrices[j])->getLocalNumEntries());
        TEUCHOS_TEST_FOR_EXCEPTION((arrayOfMatrices[j])->getDomainMap()->getIndexBase() != 0, Exceptions::RuntimeError, "interpolation subblocks must use 0 indexbase");
      }
      else arrayOfMatrices[j] = Teuchos::null;
    }
printf("after the combo do loop %d %d\n",(int)  A->getRowMap()->getLocalNumElements(), (int) nComboRowMap);
    TEUCHOS_TEST_FOR_EXCEPTION(nComboRowMap != A->getRowMap()->getLocalNumElements(), Exceptions::RuntimeError, "sum of subblock rows != #row's Afine");

    // build up csr arrays for combo block diagonal P
    Teuchos::ArrayRCP<size_t> comboPRowPtr(nComboRowMap+1);
    Teuchos::ArrayRCP<LocalOrdinal> comboPCols(nnzCombo);
    Teuchos::ArrayRCP<Scalar>       comboPVals(nnzCombo);

FILE *fp;
fp = fopen("comboArrays.m","w");
    size_t nnzCnt = 0, nrowCntFromPrevBlks = 0,ncolCntFromPrevBlks = 0;
    for(int j = 0; j < nBlks; j++) {
      // grab csr pointers for individual blocks of P
      if  (arrayOfMatrices[j] != Teuchos::null) {
        Teuchos::ArrayRCP<const size_t>     subblockRowPtr((arrayOfMatrices[j])->getLocalNumRows());
        Teuchos::ArrayRCP<const LocalOrdinal> subblockCols((arrayOfMatrices[j])->getLocalNumEntries());
        Teuchos::ArrayRCP<const Scalar>         subblockVals((arrayOfMatrices[j])->getLocalNumEntries());
        Teuchos::RCP<CrsMatrixWrap> subblockwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>((arrayOfMatrices[j]));
        Teuchos::RCP<CrsMatrix>     subblockcrs  = subblockwrap->getCrsMatrix();
        subblockcrs->getAllValues(subblockRowPtr, subblockCols, subblockVals);

        // copy jth block into csr arrays of comboP 
        
        size_t nGhost = 0;
        for (decltype(subblockRowPtr.size()) i = 0; i < subblockRowPtr.size() - 1; i++) {
          size_t rowLength = subblockRowPtr[i+1] - subblockRowPtr[i];
          comboPRowPtr[ nrowCntFromPrevBlks+i ] = nnzCnt;
fprintf(fp,"comboPRowPtr(%d) = %d;\n",  (int) (nrowCntFromPrevBlks+i+1),(int) (nnzCnt+1)  );
          for (size_t k = 0; k < rowLength; k++) {
            if  ( (int) subblockCols[k+subblockRowPtr[i]] < (int) DomMapSizePerBlk[j]) comboPCols[nnzCnt  ] = subblockCols[k+subblockRowPtr[i]] + ncolCntFromPrevBlks;
            else  {   // stick it after all nonGhostts
              comboPCols[nnzCnt  ] = subblockCols[k+subblockRowPtr[i]]-DomMapSizePerBlk[j] + nGhost + nComboDomMap; 
              printf("NOT HERE\n"); 
              if ( (int) comboPCols[nnzCnt] >= (int) nComboColMap) { printf("ERROR\n"); exit(1); }
            } 
            comboPVals[nnzCnt++] = subblockVals[k+subblockRowPtr[i]];
fprintf(fp,"comboPCols(%d) = %d;   comboPVals(%d) = %20.13e;\n",(int) nnzCnt, (int) (comboPCols[nnzCnt-1]+1), (int) nnzCnt, (double) comboPVals[nnzCnt-1]);
          }
        }
        nrowCntFromPrevBlks += Teuchos::as<size_t>((arrayOfMatrices[j])->getRowMap()->getLocalNumElements());
        ncolCntFromPrevBlks += DomMapSizePerBlk[j];
        nGhost += (ColMapSizePerBlk[j] - DomMapSizePerBlk[j]);
      }
    }
    comboPRowPtr[nComboRowMap] = nnzCnt; 
fprintf(fp,"comboPRowPtr(%d) = %d;\n",(int) nComboRowMap+1, (int) nnzCnt+1);
fclose(fp);


    // Come up with global IDs for the coarse grid maps. We assume that each xxx
    // block has a minimum GID of 0.  Since MueLu is generally creating these 
    // GIDS, this assumption is probably correct, but we'll check it.


    Teuchos::Array<GlobalOrdinal> comboDomainMapGIDs(nComboDomMap);
    Teuchos::Array<GlobalOrdinal> comboColMapGIDs(nComboColMap);

    GlobalOrdinal offset = 0; 
    size_t        domainMapIndex = 0, colMapIndex = 0;
    size_t nGhost = 0;
    for(int j = 0; j < nBlks; j++) {
      if  (arrayOfMatrices[j] != Teuchos::null) {
        TEUCHOS_TEST_FOR_EXCEPTION(arrayOfMatrices[j]->getDomainMap()->getMinAllGlobalIndex() < 0,Exceptions::RuntimeError, "Global ID assumption for domainMap not met within subblock");
        GlobalOrdinal maxGID = (arrayOfMatrices[j])->getDomainMap()->getMaxAllGlobalIndex() + 1;
        for(size_t c = 0; c < DomMapSizePerBlk[j]; ++c) {
          comboColMapGIDs[   domainMapIndex  ] =  offset + (arrayOfMatrices[j])->getColMap()->getGlobalElement(c);
          comboDomainMapGIDs[domainMapIndex++] =  offset + (arrayOfMatrices[j])->getDomainMap()->getGlobalElement(c);
        }
        
        for(size_t c = 0; c <  ColMapSizePerBlk[j] - DomMapSizePerBlk[j]; ++c) {
//        for(size_t c = 0; c <  ColMapSizePerBlk[j]; ++c) {
//          comboColMapGIDs[colMapIndex++] =  offset + (arrayOfMatrices[j])->getColMap()->getGlobalElement(c);
          comboColMapGIDs[nGhost + nComboDomMap + colMapIndex++] =  offset + (arrayOfMatrices[j])->getColMap()->getGlobalElement(c+DomMapSizePerBlk[j]);
        }
        nGhost += (ColMapSizePerBlk[j] - DomMapSizePerBlk[j]);
        offset += maxGID+1;
      }
    }

    RCP<const Map> coarseDomainMap = Xpetra::MapFactory<LO,GO,NO>::Build(A->getDomainMap()->lib(),Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),comboDomainMapGIDs, 0, A->getDomainMap()->getComm());
    RCP<const Map> coarseColMap    = Xpetra::MapFactory<LO,GO,NO>::Build(A->getDomainMap()->lib(),Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),comboColMapGIDs,    0, A->getDomainMap()->getComm());


    Teuchos::RCP<CrsMatrix> comboPCrs = CrsMatrixFactory::Build(A->getRowMap(), coarseColMap,nnzCombo);

    for (size_t i = 0; i < nComboRowMap; i++) {
      comboPCrs->insertLocalValues(i, comboPCols.view(comboPRowPtr[i], comboPRowPtr[i+1] - comboPRowPtr[i]),
                                    comboPVals.view(comboPRowPtr[i], comboPRowPtr[i+1] - comboPRowPtr[i]));
    }
    comboPCrs->fillComplete(coarseDomainMap, A->getRowMap());

    Teuchos::RCP<Matrix> comboP = Teuchos::rcp(new CrsMatrixWrap(comboPCrs));

    Set(coarseLevel, "P", comboP);
Xpetra::IO<SC,LO,GO,Node>::Write( "combo",*comboP);
  }


} //namespace MueLu

#define MUELU_COMBINEPFACTORY_SHORT
#endif // MUELU_COMBINEPFACTORY_DEF_HPP
