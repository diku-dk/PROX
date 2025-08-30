#ifndef SPARSE_VECTOR_H
#define SPARSE_VECTOR_H

#include <sparse_index_iterator.h> // for index iterator
#include <sparse_traits.h>         // for zero_block, volatile_zero_block

#include <cassert>
#include <vector>

namespace sparse
{

  namespace detail
  {
    /**
     * Vector Accessor Class.
     * This class can be used to obtain direct access to the underlying data container type used in the vector class.
     *
     * No guarantees of any sort is given once this accessor is used by an end-user.
     *
     * @tparam V     The vector type that the accessor should grant access to.
     */
    template <typename V>
    class VectorAccessor
      {
      public:
          using vector_type = V;
          using data_container_type = typename V::data_container_type;

          static data_container_type& data(vector_type& src) { return src.m_data; }
          static data_container_type const& data(vector_type const& src) { return src.m_data; }
      };

  }//namespace detail


  /**
   * Block Vector Class.
   *
   * @tparam B    The block type of the vector.
   */
  template <typename B>
  class Vector
    {
    protected:
        using data_container_type = std::vector<B>;

    public:
        using vector_type = Vector<B>;
        using block_type = B;
        using reference = block_type&;
        using const_reference = const block_type&;
        using pointer = block_type*;
        using const_pointer = const block_type*;
        using iterator = IndexIterator<false, vector_type>;
        using const_iterator = IndexIterator<true, vector_type>;

        using accessor = detail::VectorAccessor<vector_type>;

    private:
        friend class detail::VectorAccessor<vector_type>;

    protected:

      data_container_type m_data;

    public:

      Vector(size_t const n=0u)
			: m_data(n, zero_block<block_type>())
      {}

      Vector(vector_type const & orig)
  		: m_data(orig.m_data)
      {}

      ~Vector(){}

      // 2009-07-01 Kenny: Why not use this-> for accessing members/methods in template class?
      vector_type& operator=(vector_type const& orig) { m_data = orig.m_data;  return *this; }

      // 2009-07-01 Kenny: Why not use this-> for accessing members/methods in template class?
      reference operator[](size_t const i)
      {
        assert(i < size() || !"i was too large");
        return m_data[i];
      }

      // 2009-07-01 Kenny: Why not use this-> for accessing members/methods in template class?
      const_reference operator[](size_t const i) const
      {
        assert(i < size() || !"i was too large");
        return m_data[i];
      }

      // 2009-07-01 Kenny: Why not use this-> for accessing members/methods in template class?
      reference       operator()(size_t const i)       { return operator[](i); }
      const_reference operator()(size_t const i) const { return operator[](i); }

      iterator       begin()       { return iterator(0u, this);           }
      const_iterator begin() const { return const_iterator(0u, this);     }
      // 2009-07-01 Kenny: Why not use this-> for accessing members/methods in template class?
      iterator       end()         { return iterator(size(), this);       }
      const_iterator end()   const { return const_iterator(size(), this); }

      void resize(size_t const n, bool preserve = true)
      {
        if (!preserve)
        {
          assert(false || !"preserve is not currently supported");
        }

        // 2009-07-01 Kenny: Why not use this-> for accessing members/methods in template class?
        m_data.resize(n, zero_block<block_type>());
      }

      // 2009-07-01 Kenny: Why not use this-> for accessing members/methods in template class?
      void reserve(size_t const n) { m_data.reserve(n); }

      // 2009-07-01 Kenny: Why not use this-> for accessing members/methods in template class?
      size_t capacity() const { return m_data.capacity(); }

      // 2009-07-01 Kenny: Why not use this-> for accessing members/methods in template class?
      void clear() { m_data.clear(); }

      /**
       * Sets all data in the vector to zero
       */
      void clear_data()
      {
        for( typename data_container_type::iterator iter = m_data.begin(); iter != m_data.end(); ++iter )
        {
          (*iter).clear_data();
        }
      }

      // 2009-07-01 Kenny: Why not use this-> for accessing members/methods in template class?
      size_t size()  const { return m_data.size(); }
      size_t nrows() const { return size();        }
      size_t ncols() const { return 1u;            }

    };

} // namespace sparse

// SPARSE_VECTOR_H
#endif
