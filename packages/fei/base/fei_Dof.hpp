/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_Dof_hpp_
#define _fei_Dof_hpp_

#include <fei_macros.hpp>

namespace fei {

/** Dof - mesh-degree-of-freedom.
 *
 * A mesh-dof is the triple (rank, id, field).
 * - Rank is an integer type used to label different 'kinds' of mesh-dofs,
 *   e.g. node vs face, etc.
 * - Id identifies a particular instance of a rank, e.g., node 97.
 * - Field is an integer type used to label solution fields such as
 *   temperature or displacement, etc.
 *
 * Thus if the user chooses to give nodes a rank of 0 and the temperature
 * field a label of 4, then the temperature at node 97 can be represented
 * as the dof (0, 97, 4).
 *
 * Notes:
 *
 * 1. The Dof class is templated on the two integer types LocalOrdinal
 *    and GlobalOrdinal. Ranks and Fields have type LocalOrdinal, while
 *    ids have type GlobalOrdinal. The distinction is somewhat arbitrary,
 *    but the assumption is that there may be billions (or more?) ids, but
 *    probably far fewer distinct ranks or fields. So in extreme cases a
 *    user may wish to use a different type (larger) for ids than for the
 *    rank or field.
 */
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
template<class LocalOrdinal, class GlobalOrdinal>
#endif
class Dof {
 public:
#ifndef TPETRA_ENABLE_TEMPLATE_ORDINALS
  using LocalOrdinal = typename Tpetra::Map<>::local_ordinal_type;
  using GlobalOrdinal = typename Tpetra::Map<>::global_ordinal_type;
#endif
  /** constructor */
  Dof(LocalOrdinal rank, GlobalOrdinal id, LocalOrdinal field)
   : m_rank(rank), m_id(id), m_field(field) {}

  /** destructor */
  ~Dof(){}

  LocalOrdinal rank() const { return m_rank; }
  GlobalOrdinal id() const { return m_id; }
  LocalOrdinal field() const { return m_field; }

 private:
  LocalOrdinal m_rank;
  GlobalOrdinal m_id;
  LocalOrdinal m_field;
};//class Dof

/** Less operator which will order Dofs by rank first, then id then field.
 */
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
template<class LocalOrdinal, class GlobalOrdinal>
#endif
struct less_rank_id_field {
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  bool operator()(const Dof<LocalOrdinal,GlobalOrdinal>& dof1,
                 const Dof<LocalOrdinal,GlobalOrdinal>& dof2) const
#else
  using LocalOrdinal = typename Tpetra::Map<>::local_ordinal_type;
  using GlobalOrdinal = typename Tpetra::Map<>::global_ordinal_type;
  bool operator()(const Dof<>& dof1,
                 const Dof<>& dof2) const
#endif
  {
    if (dof1.rank()==dof2.rank()) {
      if (dof1.id() == dof2.id()) return dof1.field() < dof2.field();
      else return dof1.id() < dof2.id();
    }
    else {
      return dof1.rank() < dof2.rank();
    }
  }

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  bool operator()(const Dof<LocalOrdinal,GlobalOrdinal>* dof1,
                 const Dof<LocalOrdinal,GlobalOrdinal>* dof2) const
#else
  bool operator()(const Dof<>* dof1,
                 const Dof<>* dof2) const
#endif
  {
    if (dof1->rank()==dof2->rank()) {
      if (dof1->id() == dof2->id()) return dof1->field() < dof2->field();
      else return dof1->id() < dof2->id();
    }
    else {
      return dof1->rank() < dof2->rank();
    }
  }
};

/** Less operator which will order Dofs by field first, then rank then id.
 */
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
template<class LocalOrdinal, class GlobalOrdinal>
#endif
struct less_field_rank_id {
#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  bool operator()(const Dof<LocalOrdinal,GlobalOrdinal>& dof1,
                 const Dof<LocalOrdinal,GlobalOrdinal>& dof2) const
#else
  using LocalOrdinal = typename Tpetra::Map<>::local_ordinal_type;
  using GlobalOrdinal = typename Tpetra::Map<>::global_ordinal_type;
  bool operator()(const Dof<>& dof1,
                 const Dof<>& dof2) const
#endif
  {
    if (dof1.field()==dof2.field()) {
      if (dof1.rank() == dof2.rank()) return dof1.id() < dof2.id();
      else return dof1.rank() < dof2.rank();
    }
    else {
      return dof1.field() < dof2.field();
    }
  }

#ifdef TPETRA_ENABLE_TEMPLATE_ORDINALS
  bool operator()(const Dof<LocalOrdinal,GlobalOrdinal>* dof1,
                 const Dof<LocalOrdinal,GlobalOrdinal>* dof2) const
#else
  bool operator()(const Dof<>* dof1,
                 const Dof<>* dof2) const
#endif
  {
    if (dof1->field()==dof2->field()) {
      if (dof1->rank() == dof2->rank()) return dof1->id() < dof2->id();
      else return dof1->rank() < dof2->rank();
    }
    else {
      return dof1->field() < dof2->field();
    }
  }
};

}//namespace fei
#endif

