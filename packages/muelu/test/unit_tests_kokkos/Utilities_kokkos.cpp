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
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_UseDefaultTypes.hpp"

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Utilities_kokkos.hpp"
#include "MueLu_Version.hpp"

namespace MueLuTests {


#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Utilities_kokkos, CuthillMcKee, Scalar, LocalOrdinal, GlobalOrdinal, Node)
#else
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Utilities_kokkos, CuthillMcKee, Scalar, Node)
#endif
  {
#   include "MueLu_UseShortNames.hpp"
#ifndef TPETRA_ENABLE_TEMPLATE_ORDINALS
    using LocalOrdinal = typename Tpetra::Map<>::local_ordinal_type;
    using GlobalOrdinal = typename Tpetra::Map<>::global_ordinal_type;
#endif
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;
    
    using Teuchos::RCP;
    using Teuchos::rcp;
    using real_type = typename Teuchos::ScalarTraits<SC>::coordinateType;
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
    using RealValuedMultiVector = Xpetra::MultiVector<real_type,LO,GO,NO>;
#else
    using RealValuedMultiVector = Xpetra::MultiVector<real_type,NO>;
#endif

    // Build the problem
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
    RCP<Matrix> A = MueLuTests::TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(2001);
#else
    RCP<Matrix> A = MueLuTests::TestHelpers_kokkos::TestFactory<SC, NO>::Build1DPoisson(2001);
#endif
    RCP<const Map> map = A->getMap();
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
    RCP<RealValuedMultiVector> coordinates = Xpetra::MultiVectorFactory<real_type, LO, GO, NO>::Build(map, 1);
#else
    RCP<RealValuedMultiVector> coordinates = Xpetra::MultiVectorFactory<real_type,NO>::Build(map, 1);
#endif
    RCP<MultiVector> nullspace = MultiVectorFactory::Build(map, 1);
    nullspace->putScalar(Teuchos::ScalarTraits<SC>::one());

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
    MueLu::Utilities_kokkos<SC,LO,GO,NO>::Transpose(*A);// compile test
#else
    MueLu::Utilities_kokkos<SC,NO>::Transpose(*A);// compile test
#endif

    // CM Test
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
    RCP<Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> > ordering
      = MueLu::Utilities_kokkos<SC,LO,GO,NO>::CuthillMcKee(*A);
#else
    RCP<Xpetra::Vector<LocalOrdinal,Node> > ordering
      = MueLu::Utilities_kokkos<SC,NO>::CuthillMcKee(*A);
#endif


#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
    RCP<Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> > ordering2
      = MueLu::Utilities_kokkos<SC,LO,GO,NO>::ReverseCuthillMcKee(*A);
#else
    RCP<Xpetra::Vector<LocalOrdinal,Node> > ordering2
      = MueLu::Utilities_kokkos<SC,NO>::ReverseCuthillMcKee(*A);
#endif


    TEST_EQUALITY(1,1);
  }


#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
#define MUELU_ETI_GROUP(SC,LO,GO,NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Utilities_kokkos, CuthillMcKee, SC, LO, GO, NO)
#else
#define MUELU_ETI_GROUP(SC,NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(Utilities_kokkos, CuthillMcKee, SC, NO)
#endif

#include <MueLu_ETI_4arg.hpp>


} //namespace MueLuTests
