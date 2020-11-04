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
#ifndef MUELU_STRATIMIKOSSMOOTHER_DEF_HPP
#define MUELU_STRATIMIKOSSMOOTHER_DEF_HPP

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_THYRA)

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_StratimikosSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include "Teuchos_AbstractFactoryStd.hpp"
#if defined(HAVE_MUELU_THYRATPETRAADAPTERS) && defined(HAVE_MUELU_IFPACK2)
#include <Thyra_Ifpack2PreconditionerFactory.hpp>
#endif


namespace MueLu {

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::StratimikosSmoother(const std::string type, const Teuchos::ParameterList& paramList)
#else
  template <class Node>
  StratimikosSmoother<double, Node>::StratimikosSmoother(const std::string type, const Teuchos::ParameterList& paramList)
#endif
  : type_(type)
  {
    bool isSupported = type == "STRATIMIKOS";
    this->declareConstructionOutcome(!isSupported, "Stratimikos does not provide the smoother '" + type_ + "'.");
    if (isSupported)
      SetParameterList(paramList);
  }

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const Teuchos::ParameterList& paramList) {
#else
  template <class Node>
  void StratimikosSmoother<double, Node>::SetParameterList(const Teuchos::ParameterList& paramList) {
#endif
    Factory::SetParameterList(paramList);
  }

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
#else
  template <class Node>
  void StratimikosSmoother<double, Node>::DeclareInput(Level& currentLevel) const {
#endif
    this->Input(currentLevel, "A");
  }

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
#else
  template <class Node>
  void StratimikosSmoother<double, Node>::Setup(Level& currentLevel) {
#endif
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);

    A_ = Factory::Get< RCP<Matrix> >(currentLevel, "A");
    SetupStratimikos(currentLevel);
    SmootherPrototype::IsSetup(true);
    this->GetOStream(Statistics1) << description() << std::endl;
  }

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::SetupStratimikos(Level& /* currentLevel */) {
#else
  template <class Node>
  void StratimikosSmoother<double, Node>::SetupStratimikos(Level& /* currentLevel */) {
#endif

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
    RCP<const Thyra::LinearOpBase<Scalar> > thyraA = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A_)->getCrsMatrix());
#else
    RCP<const Thyra::LinearOpBase<Scalar> > thyraA = Xpetra::ThyraUtils<Scalar,Node>::toThyra(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A_)->getCrsMatrix());
#endif

    // Build Stratimikos solver
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    typedef Thyra::PreconditionerFactoryBase<Scalar> Base;
#if defined(HAVE_MUELU_THYRATPETRAADAPTERS) && defined(HAVE_MUELU_IFPACK2)
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
    typedef Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Impl;
#else
    typedef Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<Scalar, Node> > Impl;
#endif
    linearSolverBuilder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
#endif
    linearSolverBuilder.setParameterList(rcpFromRef(const_cast<ParameterList&>(this->GetParameterList())));


    // Build a new "solver factory" according to the previously specified parameter list.
    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);
    solver_ = Thyra::linearOpWithSolve(*solverFactory, thyraA);
  }

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
#else
  template <class Node>
  void StratimikosSmoother<double, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
#endif
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::StratimikosSmoother::Apply(): Setup() has not been called");

    // Apply
    if (InitialGuessIsZero) {
      X.putScalar(0.0);
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
      RCP<      Thyra::VectorBase<Scalar> >thyraX = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyraVector(X.getVectorNonConst(0)));
      RCP<const Thyra::VectorBase<Scalar> >thyraB = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyraVector(B.getVector(0));
#else
      RCP<      Thyra::VectorBase<Scalar> >thyraX = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(Xpetra::ThyraUtils<Scalar,Node>::toThyraVector(X.getVectorNonConst(0)));
      RCP<const Thyra::VectorBase<Scalar> >thyraB = Xpetra::ThyraUtils<Scalar,Node>::toThyraVector(B.getVector(0));
#endif
      Thyra::SolveStatus<Scalar> status = Thyra::solve<Scalar>(*solver_, Thyra::NOTRANS, *thyraB, thyraX.ptr());

    } else {
      typedef Teuchos::ScalarTraits<Scalar> TST;
      RCP<MultiVector> Residual   = Utilities::Residual(*A_, X, B);
      RCP<MultiVector> Correction = MultiVectorFactory::Build(A_->getDomainMap(), X.getNumVectors());

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
      RCP<      Thyra::VectorBase<Scalar> >thyraX = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyraVector(Correction->getVectorNonConst(0)));
      RCP<const Thyra::VectorBase<Scalar> >thyraB = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyraVector(Residual->getVector(0));
#else
      RCP<      Thyra::VectorBase<Scalar> >thyraX = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(Xpetra::ThyraUtils<Scalar,Node>::toThyraVector(Correction->getVectorNonConst(0)));
      RCP<const Thyra::VectorBase<Scalar> >thyraB = Xpetra::ThyraUtils<Scalar,Node>::toThyraVector(Residual->getVector(0));
#endif
      Thyra::SolveStatus<Scalar> status = Thyra::solve<Scalar>(*solver_, Thyra::NOTRANS, *thyraB, thyraX.ptr());

      X.update(TST::one(), *Correction, TST::one());
    }

  }

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<MueLu::SmootherPrototype<double, LocalOrdinal, GlobalOrdinal, Node> > StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
#else
  template <class Node>
  RCP<MueLu::SmootherPrototype<double, Node> > StratimikosSmoother<double, Node>::Copy() const {
#endif
    RCP<StratimikosSmoother> smoother = rcp(new StratimikosSmoother(*this) );
    smoother->SetParameterList(this->GetParameterList());
    return smoother;
  }

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::description() const {
#else
  template <class Node>
  std::string StratimikosSmoother<double, Node>::description() const {
#endif
    std::ostringstream out;
    if (SmootherPrototype::IsSetup()) {
      out << solver_->description();
    } else {
      out << "STRATIMIKOS {type = " << type_ << "}";
    }
    return out.str();
  }

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
#else
  template <class Node>
  void StratimikosSmoother<double, Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
#endif
    MUELU_DESCRIBE;

    if (verbLevel & Parameters1) {
      out0 << "Parameter list: " << std::endl;
      Teuchos::OSTab tab2(out);
      out << this->GetParameterList();
    }

    if (verbLevel & External)
      if (solver_ != Teuchos::null) {
        Teuchos::OSTab tab2(out);
        out << *solver_ << std::endl;
      }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "RCP<solver_>: " << solver_ << std::endl;
    }
  }

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
#else
  template <class Node>
  size_t StratimikosSmoother<double, Node>::getNodeSmootherComplexity() const {
#endif
    return Teuchos::OrdinalTraits<size_t>::invalid();
  }


} // namespace MueLu

#endif // HAVE_MUELU_STRATIMIKOS
#endif // MUELU_STRATIMIKOSSMOOTHER_DEF_HPP
