//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

#ifndef _FROSCH_HARMONICCOARSEOPERATOR_DECL_HPP
#define _FROSCH_HARMONICCOARSEOPERATOR_DECL_HPP

#include <FROSch_CoarseOperator_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
              class LO = int,
              class GO = DefaultGlobalOrdinal,
#endif
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
    class HarmonicCoarseOperator : public CoarseOperator<SC,LO,GO,NO> {
#else
    class HarmonicCoarseOperator : public CoarseOperator<SC,NO> {
#endif

    protected:

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using XMapPtr                   = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr              = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr;
        using XMapPtrVecPtr             = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtrVecPtr;
        using XMapPtrVecPtr2D           = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtrVecPtr2D;
        using ConstXMapPtrVecPtr2D      = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtrVecPtr2D;
#else
        using LO = typename Tpetra::Map<>::local_ordinal_type;
        using GO = typename Tpetra::Map<>::global_ordinal_type;
        using XMapPtr                   = typename SchwarzOperator<SC,NO>::XMapPtr;
        using ConstXMapPtr              = typename SchwarzOperator<SC,NO>::ConstXMapPtr;
        using XMapPtrVecPtr             = typename SchwarzOperator<SC,NO>::XMapPtrVecPtr;
        using XMapPtrVecPtr2D           = typename SchwarzOperator<SC,NO>::XMapPtrVecPtr2D;
        using ConstXMapPtrVecPtr2D      = typename SchwarzOperator<SC,NO>::ConstXMapPtrVecPtr2D;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using XMatrixPtr                = typename SchwarzOperator<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr           = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMatrixPtr;
#else
        using XMatrixPtr                = typename SchwarzOperator<SC,NO>::XMatrixPtr;
        using ConstXMatrixPtr           = typename SchwarzOperator<SC,NO>::ConstXMatrixPtr;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using XMultiVectorPtr           = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtr;
        using ConstXMultiVectorPtr      = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMultiVectorPtr;
        using XMultiVectorPtrVecPtr     = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtrVecPtr;
#else
        using XMultiVectorPtr           = typename SchwarzOperator<SC,NO>::XMultiVectorPtr;
        using ConstXMultiVectorPtr      = typename SchwarzOperator<SC,NO>::ConstXMultiVectorPtr;
        using XMultiVectorPtrVecPtr     = typename SchwarzOperator<SC,NO>::XMultiVectorPtrVecPtr;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using XCrsGraph                 = typename SchwarzOperator<SC,LO,GO,NO>::XCrsGraph;
        using GraphPtr                  = typename SchwarzOperator<SC,LO,GO,NO>::GraphPtr;
        using ConstXCrsGraphPtr         = typename SchwarzOperator<SC,LO,GO,NO>::ConstXCrsGraphPtr;
#else
        using XCrsGraph                 = typename SchwarzOperator<SC,NO>::XCrsGraph;
        using GraphPtr                  = typename SchwarzOperator<SC,NO>::GraphPtr;
        using ConstXCrsGraphPtr         = typename SchwarzOperator<SC,NO>::ConstXCrsGraphPtr;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using ParameterListPtr          = typename SchwarzOperator<SC,LO,GO,NO>::ParameterListPtr;
#else
        using ParameterListPtr          = typename SchwarzOperator<SC,NO>::ParameterListPtr;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using TSerialDenseMatrixPtr     = typename SchwarzOperator<SC,LO,GO,NO>::TSerialDenseMatrixPtr;
#else
        using TSerialDenseMatrixPtr     = typename SchwarzOperator<SC,NO>::TSerialDenseMatrixPtr;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using TSerialQRDenseSolverPtr   = typename SchwarzOperator<SC,LO,GO,NO>::TSerialQRDenseSolverPtr;
#else
        using TSerialQRDenseSolverPtr   = typename SchwarzOperator<SC,NO>::TSerialQRDenseSolverPtr;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using CoarseSpacePtr            = typename SchwarzOperator<SC,LO,GO,NO>::CoarseSpacePtr;
        using CoarseSpacePtrVecPtr      = typename SchwarzOperator<SC,LO,GO,NO>::CoarseSpacePtrVecPtr;
#else
        using CoarseSpacePtr            = typename SchwarzOperator<SC,NO>::CoarseSpacePtr;
        using CoarseSpacePtrVecPtr      = typename SchwarzOperator<SC,NO>::CoarseSpacePtrVecPtr;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using EntitySetPtr              = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetPtr;
        using EntitySetConstPtr         = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetConstPtr;
        using EntitySetPtrVecPtr        = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetPtrVecPtr;
        using EntitySetPtrConstVecPtr   = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetPtrConstVecPtr;
#else
        using EntitySetPtr              = typename SchwarzOperator<SC,NO>::EntitySetPtr;
        using EntitySetConstPtr         = typename SchwarzOperator<SC,NO>::EntitySetConstPtr;
        using EntitySetPtrVecPtr        = typename SchwarzOperator<SC,NO>::EntitySetPtrVecPtr;
        using EntitySetPtrConstVecPtr   = typename SchwarzOperator<SC,NO>::EntitySetPtrConstVecPtr;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using InterfaceEntityPtr        = typename SchwarzOperator<SC,LO,GO,NO>::InterfaceEntityPtr;
        using InterfaceEntityPtrVec     = typename SchwarzOperator<SC,LO,GO,NO>::InterfaceEntityPtrVec;
        using InterfaceEntityPtrVecPtr  = typename SchwarzOperator<SC,LO,GO,NO>::InterfaceEntityPtrVecPtr;
#else
        using InterfaceEntityPtr        = typename SchwarzOperator<SC,NO>::InterfaceEntityPtr;
        using InterfaceEntityPtrVec     = typename SchwarzOperator<SC,NO>::InterfaceEntityPtrVec;
        using InterfaceEntityPtrVecPtr  = typename SchwarzOperator<SC,NO>::InterfaceEntityPtrVecPtr;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using SubdomainSolverPtr        = typename SchwarzOperator<SC,LO,GO,NO>::SubdomainSolverPtr;
#else
        using SubdomainSolverPtr        = typename SchwarzOperator<SC,NO>::SubdomainSolverPtr;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using UN                        = typename SchwarzOperator<SC,LO,GO,NO>::UN;
        using UNVec                     = typename SchwarzOperator<SC,LO,GO,NO>::UNVec;
        using UNVecPtr                  = typename SchwarzOperator<SC,LO,GO,NO>::UNVecPtr;
        using ConstUNVecView            = typename SchwarzOperator<SC,LO,GO,NO>::ConstUNVecView;
#else
        using UN                        = typename SchwarzOperator<SC,NO>::UN;
        using UNVec                     = typename SchwarzOperator<SC,NO>::UNVec;
        using UNVecPtr                  = typename SchwarzOperator<SC,NO>::UNVecPtr;
        using ConstUNVecView            = typename SchwarzOperator<SC,NO>::ConstUNVecView;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using IntVec                    = typename SchwarzOperator<SC,LO,GO,NO>::IntVec;
        using IntVec2D                  = typename SchwarzOperator<SC,LO,GO,NO>::IntVec2D;
#else
        using IntVec                    = typename SchwarzOperator<SC,NO>::IntVec;
        using IntVec2D                  = typename SchwarzOperator<SC,NO>::IntVec2D;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using LOVec                     = typename SchwarzOperator<SC,LO,GO,NO>::LOVec;
        using LOVecPtr                  = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr;
        using ConstLOVecView            = typename SchwarzOperator<SC,LO,GO,NO>::ConstLOVecView;
        using LOVecPtr2D                = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr2D;
#else
        using LOVec                     = typename SchwarzOperator<SC,NO>::LOVec;
        using LOVecPtr                  = typename SchwarzOperator<SC,NO>::LOVecPtr;
        using ConstLOVecView            = typename SchwarzOperator<SC,NO>::ConstLOVecView;
        using LOVecPtr2D                = typename SchwarzOperator<SC,NO>::LOVecPtr2D;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using GOVec                     = typename SchwarzOperator<SC,LO,GO,NO>::GOVec;
        using GOVecView                 = typename SchwarzOperator<SC,LO,GO,NO>::GOVecView;
        using GOVec2D                   = typename SchwarzOperator<SC,LO,GO,NO>::GOVec2D;
#else
        using GOVec                     = typename SchwarzOperator<SC,NO>::GOVec;
        using GOVecView                 = typename SchwarzOperator<SC,NO>::GOVecView;
        using GOVec2D                   = typename SchwarzOperator<SC,NO>::GOVec2D;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using SCVec                     = typename SchwarzOperator<SC,LO,GO,NO>::SCVec;
        using SCVecPtr                  = typename SchwarzOperator<SC,LO,GO,NO>::SCVecPtr;
        using ConstSCVecPtr             = typename SchwarzOperator<SC,LO,GO,NO>::ConstSCVecPtr;
        using ConstSCVecView            = typename SchwarzOperator<SC,LO,GO,NO>::ConstSCVecView;
#else
        using SCVec                     = typename SchwarzOperator<SC,NO>::SCVec;
        using SCVecPtr                  = typename SchwarzOperator<SC,NO>::SCVecPtr;
        using ConstSCVecPtr             = typename SchwarzOperator<SC,NO>::ConstSCVecPtr;
        using ConstSCVecView            = typename SchwarzOperator<SC,NO>::ConstSCVecView;
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        using ConstBoolVecPtr           = typename SchwarzOperator<SC,LO,GO,NO>::ConstBoolVecPtr;
#else
        using ConstBoolVecPtr           = typename SchwarzOperator<SC,NO>::ConstBoolVecPtr;
#endif

    public:

        HarmonicCoarseOperator(ConstXMatrixPtr k,
                               ParameterListPtr parameterList);

        virtual int initialize() = 0;

        XMapPtr assembleCoarseMap();

        XMapPtr computeCoarseSpace(CoarseSpacePtr coarseSpace);

    protected:

        int intializeCoarseMap();

        int assembleInterfaceCoarseSpace();

        int addZeroCoarseSpaceBlock(ConstXMapPtr dofsMap);

        int computeVolumeFunctions(UN blockId,
                                   UN dimension,
                                   ConstXMapPtr nodesMap,
                                   ConstXMultiVectorPtr nodeList,
                                   EntitySetConstPtr interior);

        virtual XMultiVectorPtrVecPtr computeTranslations(UN blockId,
                                                          EntitySetConstPtr entitySet);

        virtual XMultiVectorPtrVecPtr computeRotations(UN blockId,
                                                       UN dimension,
                                                       ConstXMultiVectorPtr nodeList,
                                                       EntitySetConstPtr entitySet,
                                                       UN discardRotations = 0);

        virtual LOVecPtr detectLinearDependencies(GOVecView indicesGammaDofsAll,
                                                        ConstXMapPtr rowMap,
                                                        ConstXMapPtr rangeMap,
                                                        ConstXMapPtr repeatedMap,
                                                        SC treshold);

        virtual XMultiVectorPtr computeExtensions(ConstXMapPtr localMap,
                                                  GOVecView indicesGammaDofsAll,
                                                  GOVecView indicesIDofsAll,
                                                  XMatrixPtr kII,
                                                  XMatrixPtr kIGamma);

        virtual ConstXMatrixPtr removeCouplingBetweenDofs(ConstXMatrixPtr matrix,
                                                          ConstXMapPtr map,
                                                          TwoDArray<int> &couplingIDsToRemove);

        virtual int buildElementNodeList();

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
        virtual int buildGlobalGraph(Teuchos::RCP<DDInterface<SC,LO,GO,NO> > theDDInterface_);
#else
        virtual int buildGlobalGraph(Teuchos::RCP<DDInterface<SC,NO> > theDDInterface_);
#endif

        virtual int buildCoarseGraph();


        SubdomainSolverPtr ExtensionSolver_;

        CoarseSpacePtrVecPtr InterfaceCoarseSpaces_ = CoarseSpacePtrVecPtr(0);
        CoarseSpacePtr AssembledInterfaceCoarseSpace_;

        UNVecPtr Dimensions_ = UNVecPtr(0);
        UNVecPtr DofsPerNode_ = UNVecPtr(0);

        LOVecPtr2D GammaDofs_ = LOVecPtr2D(0);
        LOVecPtr2D IDofs_ = LOVecPtr2D(0);

        ConstXMapPtrVecPtr2D DofsMaps_ = ConstXMapPtrVecPtr2D(0); // notwendig??

        ConstXMapPtr KRowMap_; // AH: What is this?

        UN NumberOfBlocks_ = 0;

        UN MaxNumNeigh_ = 0;
    };

}

#endif
