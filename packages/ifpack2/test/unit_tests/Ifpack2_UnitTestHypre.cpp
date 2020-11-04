/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_Hypre.hpp>


#include "Teuchos_Array.hpp"
#include "Teuchos_Comm.hpp"
#include <string>
#include <stdio.h>
#include <map>
#include <HYPRE.h>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;
using Teuchos::RCP;
using Teuchos::rcp;


#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Hypre, Construct, Scalar, LocalOrdinal, GlobalOrdinal, Node)
#else
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Ifpack2Hypre, Construct, Scalar, Node)
#endif
{
#ifndef TPETRA_ENABLE_TEMPLATE_ORDINALS
  using LocalOrdinal = typename Tpetra::Map<>::local_ordinal_type;
  using GlobalOrdinal = typename Tpetra::Map<>::global_ordinal_type;
#endif
  global_size_t num_rows_per_proc = 10;
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
#else
  const RCP<const Tpetra::Map<Node> > rowmap = tif_utest::create_tpetra_map<Node>(num_rows_per_proc);
#endif
  
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);
#else
  RCP<const Tpetra::CrsMatrix<Scalar,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,Node>(rowmap);
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec(crsmatrix);
#else
  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,Node> > prec(crsmatrix);
#endif

  prec.initialize();
  TEST_EQUALITY(0,0);
}


#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Hypre, BoomerAMG, Scalar, LocalOrdinal, GlobalOrdinal, Node)
#else
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Ifpack2Hypre, BoomerAMG, Scalar, Node)
#endif
{
#ifndef TPETRA_ENABLE_TEMPLATE_ORDINALS
  using LocalOrdinal = typename Tpetra::Map<>::local_ordinal_type;
  using GlobalOrdinal = typename Tpetra::Map<>::global_ordinal_type;
#endif
  const GlobalOrdinal num_rows_per_proc = 1000;
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  using multivector_type = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
#else
  using multivector_type = Tpetra::MultiVector<Scalar,Node>;
#endif
  const double tol = 1e-7;

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  auto rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
#else
  auto rowmap = tif_utest::create_tpetra_map<Node>(num_rows_per_proc);
#endif
  int NumProc = rowmap->getComm()->getSize();
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  auto A = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap,-Teuchos::ScalarTraits<Scalar>::one());
#else
  auto A = tif_utest::create_test_matrix<Scalar,Node>(rowmap,-Teuchos::ScalarTraits<Scalar>::one());
#endif

  // Create the parameter list
  Teuchos::ParameterList list("Preconditioner List");
  Teuchos::ParameterList precList = list.sublist("hypre: Preconditioner functions");
  precList.set("HYPRE_BoomerAMGSetPrintLevel", 0);// print amg solution info
  precList.set("HYPRE_BoomerAMGSetCoarsenType", 6);
  precList.set("HYPRE_BoomerAMGSetRelaxType", 6);  //Sym G.S./Jacobi hybrid
  precList.set("HYPRE_BoomerAMGSetNumSweeps", 1);
  precList.set("HYPRE_BoomerAMGSetTol", 0.0);      // conv. tolerance zero
  precList.set("HYPRE_BoomerAMGSetMaxIter", 1);   //do only one iteration!
  list.set("hypre: Solver", "PCG");
  list.set("hypre: Preconditioner", "BoomerAMG");
  list.set("hypre: SolveOrPrecondition", "Solver");
  list.set("hypre: SetPreconditioner", true);
  
  // Create the preconditioner
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > preconditioner(A);
#else
  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,Node> > preconditioner(A);
#endif
  preconditioner.setParameters(list);
  preconditioner.compute();
    
  int numVec = 2;
  multivector_type X(preconditioner.getRangeMap(), numVec);
  multivector_type KnownX(preconditioner.getRangeMap(), numVec);
  multivector_type B(preconditioner.getDomainMap(), numVec);
  KnownX.putScalar(Teuchos::ScalarTraits<Scalar>::one());
  preconditioner.applyMat(KnownX, B);

  preconditioner.compute();
  preconditioner.apply(B, X);

  auto v1v = X.get1dView();
  auto v2v = KnownX.get1dView();
  TEST_COMPARE_FLOATING_ARRAYS(v1v,v2v,tol*10*pow(10.0,NumProc));
}

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Hypre, BoomerAMGNonContiguous, Scalar, LocalOrdinal, GlobalOrdinal, Node)
#else
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Ifpack2Hypre, BoomerAMGNonContiguous, Scalar, Node)
#endif
{
#ifndef TPETRA_ENABLE_TEMPLATE_ORDINALS
  using LocalOrdinal = typename Tpetra::Map<>::local_ordinal_type;
  using GlobalOrdinal = typename Tpetra::Map<>::global_ordinal_type;
#endif
  const GlobalOrdinal num_rows_per_proc = 1000;
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  using multivector_type = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
#else
  using multivector_type = Tpetra::MultiVector<Scalar,Node>;
#endif
  const double tol = 1e-7;

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  auto rowmap = tif_utest::create_odd_even_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
#else
  auto rowmap = tif_utest::create_odd_even_map<Node>(num_rows_per_proc);
#endif
  int NumProc = rowmap->getComm()->getSize();
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  auto A = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap,-Teuchos::ScalarTraits<Scalar>::one());
#else
  auto A = tif_utest::create_test_matrix<Scalar,Node>(rowmap,-Teuchos::ScalarTraits<Scalar>::one());
#endif

  // Create the parameter list
  Teuchos::ParameterList list("Preconditioner List");
  Teuchos::ParameterList precList = list.sublist("hypre: Preconditioner functions");
  precList.set("HYPRE_BoomerAMGSetPrintLevel", 0);// print amg solution info
  precList.set("HYPRE_BoomerAMGSetCoarsenType", 6);
  precList.set("HYPRE_BoomerAMGSetRelaxType", 6);  //Sym G.S./Jacobi hybrid
  precList.set("HYPRE_BoomerAMGSetNumSweeps", 1);
  precList.set("HYPRE_BoomerAMGSetTol", 0.0);      // conv. tolerance zero
  precList.set("HYPRE_BoomerAMGSetMaxIter", 1);   //do only one iteration!
  list.set("hypre: Solver", "PCG");
  list.set("hypre: Preconditioner", "BoomerAMG");
  list.set("hypre: SolveOrPrecondition", "Solver");
  list.set("hypre: SetPreconditioner", true);
  

  // Create the preconditioner
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > preconditioner(A);
#else
  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,Node> > preconditioner(A);
#endif
  preconditioner.setParameters(list);
  preconditioner.compute();
    
  int numVec = 2;
  multivector_type X(preconditioner.getRangeMap(), numVec);
  multivector_type KnownX(preconditioner.getRangeMap(), numVec);
  multivector_type B(preconditioner.getDomainMap(), numVec);
  KnownX.putScalar(Teuchos::ScalarTraits<Scalar>::one());
  preconditioner.applyMat(KnownX, B);

  preconditioner.compute();
  preconditioner.apply(B, X);

  auto v1v = X.get1dView();
  auto v2v = KnownX.get1dView();
  TEST_COMPARE_FLOATING_ARRAYS(v1v,v2v,tol*10*pow(10.0,NumProc));
}


// Tests the hypre interface's ability to work with both a preconditioner and linear
// solver via ApplyInverse
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Hypre, Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  using multivector_type = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
#else
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Ifpack2Hypre, Apply, Scalar, Node) {
  using LocalOrdinal = typename Tpetra::Map<>::local_ordinal_type;
  using GlobalOrdinal = typename Tpetra::Map<>::global_ordinal_type;
  using multivector_type = Tpetra::MultiVector<Scalar,Node>;
#endif
  const double tol = 1E-7;

  global_size_t num_rows_per_proc = 10;
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
#else
  const RCP<const Tpetra::Map<Node> > rowmap = tif_utest::create_tpetra_map<Node>(num_rows_per_proc);
#endif
  
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);
#else
  RCP<const Tpetra::CrsMatrix<Scalar,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,Node>(rowmap);
#endif
  int NumProc = rowmap->getComm()->getSize();

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > preconditioner(crsmatrix);
#else
  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,Node> > preconditioner(crsmatrix);
#endif
  preconditioner.initialize();

  //
  // Create the solution vector and RHS
  //
  int numVec = 2;
  multivector_type X(preconditioner.getRangeMap(), numVec);
  multivector_type KnownX(preconditioner.getRangeMap(), numVec);
  multivector_type B(preconditioner.getDomainMap(), numVec);
  KnownX.putScalar(Teuchos::ScalarTraits<Scalar>::one());
  preconditioner.compute();
  preconditioner.applyMat(KnownX, B);

  Teuchos::ParameterList list("New List");
  Teuchos::ParameterList & solverList = list.sublist("hypre: Solver functions");
  solverList.set("HYPRE_PCGSetMaxIter", 1000);               // max iterations
  solverList.set("HYPRE_PCGSetTol", 1e-9);                   // conv. tolerance
  solverList.set("HYPRE_PCGSetPrintLevel", 1);               // print solve info
  solverList.set("HYPRE_PCGSetLogging", 0);                  // needed to get run info later
  list.set("hypre: Solver functions",solverList);
  list.set("hypre: SolveOrPrecondition", "Solver");
  list.set("hypre: Solver", "PCG");
  list.set("hypre: Preconditioner", "ParaSails");
  list.set("hypre: SetPreconditioner", true);
  
  preconditioner.compute();
  preconditioner.apply(B, X);

  auto v1v = X.get1dView();
  auto v2v = KnownX.get1dView();
  TEST_COMPARE_FLOATING_ARRAYS(v1v,v2v,tol*10*pow(10.0,NumProc));
}



#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Hypre, DiagonalMatrixInOrder, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  using multivector_type = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
#else
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Ifpack2Hypre, DiagonalMatrixInOrder, Scalar, Node) {
  using LocalOrdinal = typename Tpetra::Map<>::local_ordinal_type;
  using GlobalOrdinal = typename Tpetra::Map<>::global_ordinal_type;
  using multivector_type = Tpetra::MultiVector<Scalar,Node>;
#endif
  const double tol = 1E-7;

  LocalOrdinal num_rows_per_proc = 10;
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = tif_utest::create_diagonal_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_diagonal_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(crsgraph);
#else
  RCP<Tpetra::CrsGraph<Node> > crsgraph = tif_utest::create_diagonal_graph<Node>(num_rows_per_proc);
  RCP<Tpetra::CrsMatrix<Scalar,Node> > crsmatrix = tif_utest::create_diagonal_matrix<Scalar,Node>(crsgraph);
#endif

  int NumProc = crsmatrix->getMap()->getComm()->getSize();

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > preconditioner(crsmatrix);
#else
  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,Node> > preconditioner(crsmatrix);
#endif
  preconditioner.initialize();


  //
  // Create the initial guess and RHS
  //
  int numVec = 2;
  multivector_type X(preconditioner.getRangeMap(), numVec);
  multivector_type KnownX(preconditioner.getRangeMap(), numVec);
  multivector_type B(preconditioner.getDomainMap(), numVec);
  KnownX.putScalar(Teuchos::ScalarTraits<Scalar>::one());
  preconditioner.compute();
  preconditioner.applyMat(KnownX, B);

  //
  // Create the parameter list
  //
  Teuchos::ParameterList list("Preconditioner List");
  Teuchos::ParameterList & solverList = list.sublist("hypre: Solver functions");
  solverList.set("HYPRE_PCGSetMaxIter", 1000);               // max iterations
  solverList.set("HYPRE_PCGSetTol", tol);                    // conv. tolerance
  solverList.set("HYPRE_PCGSetPrintLevel", 0);               // print solve info
  solverList.set("HYPRE_PCGSetLogging", 1);                  // needed to get run info later
  list.set("hypre: Solver", "PCG");
  list.set("hypre: SolveOrPrecondition", "Solver");
  list.set("hypre: SetPreconditioner", false);

  //
  // Create the preconditioner (which is actually a PCG solver)
  //
  preconditioner.setParameters(list);
  preconditioner.compute();

  //
  // Solve the linear system
  //
  preconditioner.apply(B,X);
  auto v1v = X.get1dView();
  auto v2v = KnownX.get1dView();
  TEST_COMPARE_FLOATING_ARRAYS(v1v,v2v,tol*10*pow(10.0,NumProc));
}

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Hypre, DiagonalMatrixNonContiguous, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  using multivector_type = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
#else
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Ifpack2Hypre, DiagonalMatrixNonContiguous, Scalar, Node) {
  using LocalOrdinal = typename Tpetra::Map<>::local_ordinal_type;
  using GlobalOrdinal = typename Tpetra::Map<>::global_ordinal_type;
  using multivector_type = Tpetra::MultiVector<Scalar,Node>;
#endif
  const double tol = 1E-7;

  LocalOrdinal num_rows_per_proc = 10;
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = tif_utest::create_odd_even_diagonal_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_diagonal_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(crsgraph);
#else
  RCP<Tpetra::CrsGraph<Node> > crsgraph = tif_utest::create_odd_even_diagonal_graph<Node>(num_rows_per_proc);
RCP<Tpetra::CrsMatrix<Scalar,Node> > crsmatrix = tif_utest::create_diagonal_matrix<Scalar,Node>(crsgraph);
#endif

  int NumProc = crsmatrix->getMap()->getComm()->getSize();

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > preconditioner(crsmatrix);
#else
  Ifpack2::Hypre<Tpetra::RowMatrix<Scalar,Node> > preconditioner(crsmatrix);
#endif
  preconditioner.initialize();


  //
  // Create the initial guess and RHS
  //
  int numVec = 2;
  multivector_type X(preconditioner.getRangeMap(), numVec);
  multivector_type KnownX(preconditioner.getRangeMap(), numVec);
  multivector_type B(preconditioner.getDomainMap(), numVec);
  KnownX.putScalar(Teuchos::ScalarTraits<Scalar>::one());
  preconditioner.compute();
  preconditioner.applyMat(KnownX, B);

  //
  // Create the parameter list
  //
  Teuchos::ParameterList list("Preconditioner List");
  Teuchos::ParameterList & solverList = list.sublist("hypre: Solver functions");
  solverList.set("HYPRE_PCGSetMaxIter", 1000);               // max iterations
  solverList.set("HYPRE_PCGSetTol", tol);                    // conv. tolerance
  solverList.set("HYPRE_PCGSetPrintLevel", 0);               // print solve info
  solverList.set("HYPRE_PCGSetLogging", 1);                  // needed to get run info later
  list.set("hypre: Solver", "PCG");
  list.set("hypre: SolveOrPrecondition", "Solver");
  list.set("hypre: SetPreconditioner", false);

  //
  // Create the preconditioner (which is actually a PCG solver)
  //
  preconditioner.setParameters(list);
  preconditioner.compute();

  //
  // Solve the linear system
  //
  preconditioner.apply(B,X);
  auto v1v = X.get1dView();
  auto v2v = KnownX.get1dView();
  TEST_COMPARE_FLOATING_ARRAYS(v1v,v2v,tol*10*pow(10.0,NumProc));
}


#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
#define UNIT_TEST_GROUP_SC_LO_GO_NO(Scalar,LocalOrdinal,GlobalOrdinal,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Ifpack2Hypre, Construct, Scalar, LocalOrdinal,GlobalOrdinal,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Ifpack2Hypre, Apply, Scalar, LocalOrdinal,GlobalOrdinal,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Ifpack2Hypre, DiagonalMatrixInOrder, Scalar, LocalOrdinal,GlobalOrdinal,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Ifpack2Hypre, DiagonalMatrixNonContiguous, Scalar, LocalOrdinal,GlobalOrdinal,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Ifpack2Hypre, BoomerAMG, Scalar, LocalOrdinal,GlobalOrdinal,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Ifpack2Hypre, BoomerAMGNonContiguous, Scalar, LocalOrdinal,GlobalOrdinal,Node) 
#else
#define UNIT_TEST_GROUP_SC_LO_GO_NO(Scalar,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Ifpack2Hypre, Construct, Scalar,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Ifpack2Hypre, Apply, Scalar,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Ifpack2Hypre, DiagonalMatrixInOrder, Scalar,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Ifpack2Hypre, DiagonalMatrixNonContiguous, Scalar,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Ifpack2Hypre, BoomerAMG, Scalar,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Ifpack2Hypre, BoomerAMGNonContiguous, Scalar,Node) 
#endif


#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLGN_REAL( UNIT_TEST_GROUP_SC_LO_GO_NO )

} // namespace (anonymous)
