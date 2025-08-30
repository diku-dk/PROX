#ifndef SPARSE_COMPRESSED_VECTOR_ITERATOR_H
#define SPARSE_COMPRESSED_VECTOR_ITERATOR_H

#include <sparse_if_then_else.h>
#include <sparse_property_maps.h>  // for col and row functions

#include <iterator>     // for iterator_traits

namespace sparse
{
  /**
   * Compressed Vector Row Iterator Class.
   *
   * @tparam is_const     Compile time boolean indicating it the iterator should const or non-const.
   * @tparam V            The vector type the iterator should work on.
   */
  template <bool is_const, typename V>
  class CompressedVectorRowIterator
    {
    private:

// 2010-04-19 mrtn: VS 2005 view the template friend prototype as a function declaration, causing compile error c2668
#if _MSC_VER != 1400
      template<bool b, typename T> friend size_t row(CompressedVectorRowIterator<b, T> const&);
      template<bool b, typename T> friend size_t col(CompressedVectorRowIterator<b, T> const&);
#else
      friend size_t row<is_const,V>(CompressedVectorRowIterator<is_const,V> const&);
      friend size_t col<is_const,V>(CompressedVectorRowIterator<is_const,V> const&);
#endif

    private:
        using vector_type = V;

        using vector_type_ptr = typename detail::if_then_else<is_const, const vector_type*, vector_type*>::result;
        using data_container_type = typename vector_type::accessor::data_container_type;
        using data_iterator = typename data_container_type::iterator;
        using const_data_iterator = typename data_container_type::const_iterator;
        using data_iterator_type = typename detail::if_then_else<is_const, const_data_iterator, data_iterator>::result;

    public:
        using iterator_category = typename std::iterator_traits<data_iterator_type>::iterator_category;

        using value_type = typename data_iterator_type::value_type;
        using difference_type = typename data_iterator_type::difference_type;
        using pointer = typename data_iterator_type::pointer;
        using reference = typename data_iterator_type::reference;

    protected:
        vector_type_ptr m_src;
        size_t m_idx;

    public:

      CompressedVectorRowIterator()
      : m_src(0)
      , m_idx(0u)
      {}

      CompressedVectorRowIterator(size_t r, vector_type_ptr m)
      : m_src(m)
      , m_idx(r)
      {}

      CompressedVectorRowIterator(CompressedVectorRowIterator const& orig)
      : m_src(orig.m_src)
      , m_idx(orig.m_idx)
      {}

      CompressedVectorRowIterator& operator=(CompressedVectorRowIterator const & rhs)
      {
        this->m_src = rhs.m_src;
        this->m_idx = rhs.m_idx;
        return *this;
      }

      reference operator*() const
      {
        typedef typename vector_type::accessor A;
        return A::data(*(this->m_src))[m_idx];
      }

      // 2009-07-01 Kenny: Why not ue this-> for accessing members/methods in a template class?
      pointer operator->() const { return &operator*(); }

      // 2009-07-01 Kenny: Why not ue this-> for accessing members/methods in a template class?
      CompressedVectorRowIterator & operator++()
      {
        ++m_idx;
        return *this;
      }

      // 2009-07-01 Kenny: Why not ue this-> for accessing members/methods in a template class?
      CompressedVectorRowIterator operator++(int)
      {
        CompressedVectorRowIterator tmp = *this;
        ++m_idx;
        return tmp;
      }

      // 2009-07-01 Kenny: Why not ue this-> for accessing members/methods in a template class?
      CompressedVectorRowIterator& operator--()
      {
        --m_idx;
        return *this;
      }

      // 2009-07-01 Kenny: Why not ue this-> for accessing members/methods in a template class?
      CompressedVectorRowIterator operator--(int)
      {
        CompressedVectorRowIterator tmp = *this;
        --m_idx;
        return tmp;
      }

      // 2009-07-01 Kenny: Why not ue this-> for accessing members/methods in a template class?
      CompressedVectorRowIterator& operator+=(int n)
      {
        m_idx += n;
        return *this;
      }

      // 2009-07-01 Kenny: Why not ue this-> for accessing members/methods in a template class?
      CompressedVectorRowIterator& operator-=(int n)
      {
        m_idx -= n;
        return *this;
      }

      // 2009-07-01 Kenny: Why not ue this-> for accessing members/methods in a template class?
      CompressedVectorRowIterator operator+(int n)       { return CompressedVectorRowIterator(m_idx + n, m_src); }
      CompressedVectorRowIterator operator-(int n) const { return CompressedVectorRowIterator(m_idx - n, m_src); }

      // 2009-07-01 Kenny: Why not ue this-> for accessing members/methods in a template class?
      difference_type operator-(CompressedVectorRowIterator const & x) const
      {
        // 2009-07-01 Kenny: Static_cast type safe?
        return   static_cast<difference_type>(m_idx) - static_cast<difference_type>(x.m_idx);
      }

      // 2009-07-01 Kenny: Why not ue this-> for accessing members/methods in a template class?
      bool operator != (CompressedVectorRowIterator const & y) const { return !operator == (y); }

      // 2009-07-01 Kenny: Why not ue this-> for accessing members/methods in a template class?
      // 2009-07-01 Kenny: Why not test /assert if src's are the same?
      bool operator < (CompressedVectorRowIterator const & y) const
      {
        assert(m_idx == y.m_idx  || !"rows did not match");
        return m_idx < y.m_idx;
      }

      // 2009-07-01 Kenny: Why not ue this-> for accessing members/methods in a template class?
      // 2009-07-01 Kenny: Why not test /assert if src's are the same?
      bool operator == (CompressedVectorRowIterator const & y) const
      {
        // 2010-04-01 Kenny: Would there be any precendence problems with these logic tests?
        return m_idx == y.m_idx && m_idx == y.m_idx;
      }

    };

  template <bool is_const, typename V>
  inline size_t row(CompressedVectorRowIterator<is_const, V> const& iter) { return iter.m_src->row_idx(iter.m_idx); }

  template <bool is_const, typename V>
  inline size_t col(CompressedVectorRowIterator<is_const, V> const& iter) { return 0u; }

} // namespace sparse

// SPARSE_COMPRESSES_VECTOR_ITERATOR_H
#endif
