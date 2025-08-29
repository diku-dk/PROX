#ifndef OPENTISSUE_CORE_CONTAINERS_GRID_UTIL_ITERATORS_H
#define OPENTISSUE_CORE_CONTAINERS_GRID_UTIL_ITERATORS_H
//
// OpenTissue Template Library
// - A generic toolbox for physics-based modeling and simulation.
// Copyright (C) 2008 Department of Computer Science, University of Copenhagen.
//
// OTTL is licensed under zlib: http://opensource.org/licenses/zlib-license.php
//

#include <iterator>
#include <eigen3/Eigen/Dense>

namespace grid
{
    namespace detail
    {

        //! basic iterator for walking through the WHOLE container
        template <class grid_type_, class reference_type, class pointer_type>
        class Iterator
            : public std::iterator< std::random_access_iterator_tag, typename grid_type_::value_type >
        {
        public:

            typedef grid_type_                                            grid_type;

        private:

            typedef Iterator<grid_type, reference_type, pointer_type>     self_type;

        protected:

            grid_type *   m_grid;
            pointer_type m_pos;

        public:

            grid_type       & get_grid()       { return *m_grid; }
            grid_type const & get_grid() const { return *m_grid; }

            pointer_type const & get_pointer() const { return m_pos; }
            pointer_type       & get_pointer()       { return m_pos; }

        public:

            Iterator()
                : m_grid( nullptr )
                , m_pos( nullptr )
            {}

            Iterator( grid_type * grid, pointer_type pos )
                : m_grid( grid )
                , m_pos( pos )
            {}

            reference_type operator*()
            {
                return* m_pos;
            }

            self_type operator++( int )
            {
                self_type tmp = *this;
                m_pos++;
                return tmp;
            }

            self_type &operator++()
            {
                m_pos++;
                return *this;
            }

            void operator+=( size_t m )
            {
                m_pos += m;
            }

            self_type operator+ ( size_t m ) const
            {
                self_type tmp = *this;
                tmp.m_pos += m;
                return tmp;
            }

            self_type operator- ( size_t m ) const
            {
                self_type tmp = *this;
                tmp.m_pos -= m;
                return tmp;
            }

            self_type operator--( int )
            {
                self_type tmp = *this;
                m_pos--;
                return tmp;
            }

            self_type &operator--()
            {
                m_pos--;
                return *this;
            }

            void operator-=( size_t m )                       {   m_pos -= m;                 }
            int operator-  ( self_type const & other ) const  {  return static_cast<int>(m_pos - other.m_pos);   }
            bool operator< ( self_type const & other ) const  {  return m_pos < other.m_pos;  }
            bool operator<=( self_type const & other ) const  {  return m_pos <= other.m_pos; }
            bool operator> ( self_type const & other ) const  {  return m_pos > other.m_pos;  }
            bool operator>=( self_type const & other ) const  {  return m_pos >= other.m_pos; }

            //Rasmus: Delete old comments...?
            // TODO: henrikd 2005-06-27 - confirm these!
            Eigen::Matrix<size_t, 3, 1> compute_index() const
            {
                long offset = m_pos - m_grid->data();
                Eigen::Matrix<size_t, 3, 1> res;
                size_t rowLen = m_grid->m_nodes.x();
                size_t planeLen = m_grid->m_nodes.x() * m_grid->m_nodes.y();
                res.z() = static_cast<size_t>( offset / static_cast<long>(planeLen) );
                offset = offset % static_cast<long>(planeLen);
                res.y() = static_cast<size_t>( offset / static_cast<long>(rowLen) );
                res.x() = static_cast<size_t>( offset % static_cast<long>(rowLen) );
                return res;
            }

            void operator+=(const Eigen::Matrix<size_t, 3, 1>& idx )
            {
                m_pos += idx.x() + m_grid->m_nodes.x() * idx.y() + m_grid->m_nodes.x() * m_grid->m_nodes.y() * idx.z();
            }

            self_type const & operator=(const Eigen::Matrix<size_t, 3, 1>& idx )
            {
                m_pos = m_grid->data() + idx(0) + m_grid->m_nodes.x() * idx(1) + m_grid->m_nodes.x() * m_grid->m_nodes.y() * idx(2);
                return *this;
            }

            bool operator!=( self_type const & other ) const
            {
                return m_pos != other.m_pos;
            }

            bool operator==( self_type const & other ) const
            {
                return m_pos == other.m_pos;
            }
        };


        //! iterator that keeps track of i,j,k position
        /*! for walking through the WHOLE container */
        template <class grid_type, class reference_type, class pointer_type>
        class IndexIterator
            : public Iterator<grid_type, reference_type, pointer_type>
        {
        private:

            typedef Iterator<grid_type, reference_type, pointer_type>      base_type;
            typedef IndexIterator<grid_type, reference_type, pointer_type> self_type;

            Eigen::Matrix<size_t, 3, 1> m_nodes;

            static size_t const m_out_of_bounds = static_cast<size_t>(-1);

            void update_indices()
            {
                Eigen::Matrix<size_t, 3, 1> v = this->compute_index();
                m_nodes.x() = v.x();
                m_nodes.y() = v.y();
                m_nodes.z() = v.z();
            }

        public:

            IndexIterator()
                : base_type()
                , m_nodes(0,0,0)
            {}

            IndexIterator( base_type const& other )
                : base_type( other )
            {
                operator=( other );
            }

            IndexIterator( grid_type * grid, pointer_type pos )
                : base_type( grid, pos )
            {
                if ( this->m_pos != this->m_grid->data() )
                {
                    base_type::operator=( base_type( this->m_grid, this->m_pos ).compute_index() );
                }
            }

            //Keeping these for backwards compat
            size_t const& i() const { return m_nodes.x(); }
            size_t const& j() const { return m_nodes.y(); }
            size_t const& k() const { return m_nodes.z(); }
            const Eigen::Matrix<size_t, 3, 1>& nodes() const {return m_nodes;}

            Eigen::Matrix<size_t, 3, 1> get_index() const
            {
                return m_nodes;
            }

            template <typename T>
            Eigen::Matrix<T, 3, 1> get_coord() const
            {
                return Eigen::Matrix<T, 3, 1>(
                    m_nodes.x() * this->m_grid->m_dir.x() + this->m_grid->min().x( ),
                    m_nodes.y() * this->m_grid->m_dir.y() + this->m_grid->min().y( ),
                    m_nodes.z() * this->m_grid->m_dir.z() + this->m_grid->min().z( )
                    );
            }

            self_type operator++( int )
            {
                self_type tmp = *this;
                ++( *this );
                return tmp;
            }

            self_type &operator++()
            {
                ++this->m_pos;
                ++m_nodes.x();
                if ( m_nodes.x() >= this->m_grid->I() )
                {
                    m_nodes.x() = 0u;
                    ++m_nodes.y();
                    if (m_nodes.y() >= this->m_grid->m_nodes.y() )
                    {
                        m_nodes.y() = 0u;
                        ++m_nodes.z();
                    }
                }
                return *this;
            }

            self_type operator--( int )
            {
                self_type tmp = *this;
                --( *this );
                return tmp;
            }

            self_type &operator--()
            {
                --this->m_pos;
                --m_nodes.x();
                if ( m_nodes.x() == m_out_of_bounds )
                {
                    m_nodes.x() = this->m_grid->I() - 1;
                    --m_nodes.y();
                    if ( m_nodes.y() == m_out_of_bounds )
                    {
                        m_nodes.y() = this->m_grid->m_nodes.y() - 1;
                        --m_nodes.z();
                    }
                }
                return *this;
            }

            // jump to new location
            self_type const & operator=(const Eigen::Matrix<size_t, 3, 1>& idx )
            {
                m_nodes.x() = idx.x();
                m_nodes.y() = idx.y();
                m_nodes.z() = idx.z();
                this->m_pos = &( *this->m_grid ) ( idx );
                return *this;
            }

            void operator+=( size_t m )
            {
                this->m_pos += m;
                update_indices();
            }

            void operator-=( size_t m )
            {
                this->m_pos -= m;
                update_indices();
            }

            self_type operator+ ( size_t m ) const
            {
                self_type tmp = *this;
                tmp.m_pos += m;
                tmp.update_indices();
                return tmp;
            }

            self_type operator- ( size_t m ) const
            {
                self_type tmp = *this;
                tmp.m_pos -= m;
                tmp.update_indices();
                return tmp;
            }

            self_type const & operator= ( base_type const & other )
            {
                this->m_grid = const_cast<grid_type*>( &( other.get_grid() ) );
                this->m_pos  = other.get_pointer();
                if ( this->m_pos == this->m_grid->data() )
                {
                    m_nodes = Eigen::Matrix<size_t, 3, 1>(0u,0u,0u);
                }
                else
                {
                    update_indices();
                }
                return *this;
            }

        };

    } // namespace detail
} // namespace grid


// OPENTISSUE_CORE_CONTAINERS_GRID_UTIL_ITERATORS_H
#endif
