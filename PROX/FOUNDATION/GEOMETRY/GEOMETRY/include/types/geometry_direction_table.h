#ifndef GEOMETRY_DIRECTION_TABLE_H
#define GEOMETRY_DIRECTION_TABLE_H

#include <eigenhelperall.h>
#include <cstddef> // Needed for size_t
#include <cassert>

namespace geometry
{

template <typename T, size_t N> class DirectionTable
{
protected:
    EigenVector3<T> m_directions[N];

public:
    EigenVector3<T> const& operator()(size_t const& idx) const
    {
      assert( idx < N || !"DirectionTable::operator(): idx was out of range");
      return this->m_directions[idx];
    }

    EigenVector3<T>& operator()(size_t const& idx)
    {
      assert( idx < N || !"DirectionTable::operator(): idx was out of range");
      return this->m_directions[idx];
    }

    size_t size() const { return N; }

  public:

    DirectionTable()
    {}
};

template <typename T, size_t N1, size_t N2>
inline DirectionTable<T, N1 + N2> make_union(DirectionTable<T, N1> const& D1,
                                             DirectionTable<T, N2> const& D2)
{
    DirectionTable<T, N1 + N2> D;

    for (size_t n = 0u; n < N1; ++n) D(n) = D1(n);
    for (size_t n = 0u; n < N2; ++n) D(n + N1) = D2(n);

    return D;
}

template <typename T> inline DirectionTable<T, 3> make3()
{

    DirectionTable<T, 3> D;

    D(0) = EigenVector3<T>(T(1), T(0), T(0));
    D(1) = EigenVector3<T>(T(0), T(1), T(0));
    D(2) = EigenVector3<T>(T(0), T(0), T(1));

    return D;
  }

  template <typename T> inline DirectionTable<T, 4> make4()
  {
      DirectionTable<T, 4> D;

      D(0) = EigenVector3<T>(1, 1, 1).normalized();
      D(1) = EigenVector3<T>(1, 1, -1).normalized();
      D(2) = EigenVector3<T>(1, -1, 1).normalized();
      D(3) = EigenVector3<T>(1, -1, -1).normalized();

      return D;
  }

  template <typename T> inline DirectionTable<T, 6> make6()
  {
      DirectionTable<T, 6> D;

      D(0) = EigenVector3<T>(1, 1, 0).normalized();
      D(1) = EigenVector3<T>(1, -1, 0).normalized();
      D(2) = EigenVector3<T>(1, 0, 1).normalized();
      D(3) = EigenVector3<T>(0, 1, 1).normalized();
      D(4) = EigenVector3<T>(-1, 0, 1).normalized();
      D(5) = EigenVector3<T>(0, -1, 1).normalized();

      return D;
  }

  template <typename T> inline DirectionTable<T, 7> make7()
  {
      return make_union(make3<T>(), make4<T>());
  }

  template <typename T> inline DirectionTable<T, 9> make9()
  {
      return make_union(make3<T>(), make6<T>());
  }

  template <typename T> inline DirectionTable<T, 10> make10()
  {
      return make_union(make4<T>(), make6<T>());
  }

  template <typename T> inline DirectionTable<T, 13> make13()
  {
      return make_union(make7<T>(), make6<T>());
  }

  template <typename T, size_t N> class DirectionTableHelper;

  template <typename T> class DirectionTableHelper<T, 3>
  {
  public:
      static DirectionTable<T, 3> make() { return make3<T>(); }
  };

  template <typename T> class DirectionTableHelper<T, 4>
  {
  public:
      static DirectionTable<T, 4> make() { return make4<T>(); }
  };

  template <typename T> class DirectionTableHelper<T, 6>
  {
  public:
      static DirectionTable<T, 6> make() { return make6<T>(); }
  };

  template <typename T> class DirectionTableHelper<T, 7>
  {
  public:
      static DirectionTable<T, 7> make() { return make7<T>(); }
  };

  template <typename T> class DirectionTableHelper<T, 9>
  {
  public:
      static DirectionTable<T, 9> make() { return make9<T>(); }
  };

  template <typename T> class DirectionTableHelper<T, 10>
  {
  public:
      static DirectionTable<T, 10> make() { return make10<T>(); }
  };

  template <typename T> class DirectionTableHelper<T, 13>
  {
  public:
      static DirectionTable<T, 13> make() { return make13<T>(); }
  };

}// namespace geometry

  namespace geometry
  {

  template<typename T, size_t N>
  class DirectionTableEigen
  {
  protected:

      EigenVector3<T> m_directions[N];

  public:

      const EigenVector3<T>& operator()(size_t const & idx) const
      {
          assert( idx < N || !"DirectionTableEigen::operator(): idx was out of range");
          return this->m_directions[idx];
      }

      EigenVector3<T>& operator()(size_t const & idx)
      {
          assert( idx < N || !"DirectionTableEigen::operator(): idx was out of range");
          return this->m_directions[idx];
      }

      size_t size() const { return N; }

  public:

      DirectionTableEigen()
      {}

  };

  template<typename T, size_t N1, size_t N2>
  inline DirectionTableEigen<T,N1+N2> make_union( DirectionTableEigen<T,N1> const & D1, DirectionTableEigen<T,N2> const & D2)
  {
      DirectionTableEigen<T,N1+N2> D;

      for(size_t n = 0u; n < N1; ++n)
          D(n)    = D1(n);
      for(size_t n = 0u; n < N2; ++n)
          D(n+N1) = D2(n);

      return D;
  }

  template <typename T> inline DirectionTableEigen<T, 3> make3Eigen()
  {
      DirectionTableEigen<T,3> D;

      D(0) =  EigenVector3<T>( 1,  0, 0 );
      D(1) =  EigenVector3<T>( 0, 1,  0 );
      D(2) =  EigenVector3<T>( 0, 0, 1  );

      return D;
  }

  template<typename T>
  inline DirectionTableEigen<T,4>  make4Eigen()
  {

      DirectionTableEigen<T,4> D;

      D(0) = unit( EigenVector3<T>(  1,  1,  1 ) );
      D(1) = unit( EigenVector3<T>(  1,  1, -1 ) );
      D(2) = unit( EigenVector3<T>(  1, -1,  1 ) );
      D(3) = unit( EigenVector3<T>(  1, -1, -1 ) );

      return D;
  }

  template <typename T> inline DirectionTableEigen<T, 6> make6Eigen()
  {
      DirectionTableEigen<T,6> D;

      D(0) = unit( EigenVector3<T>(  1,   1,  0) );
      D(1) = unit( EigenVector3<T>(  1,  -1,  0) );
      D(2) = unit( EigenVector3<T>(  1,   0, 1 ) );
      D(3) = unit( EigenVector3<T>(  0,  1,  1 ) );
      D(4) = unit( EigenVector3<T>( -1,   0, 1 ) );
      D(5) = unit( EigenVector3<T>(  0, -1,  1 ) );

      return D;
  }

  template <typename T> inline DirectionTableEigen<T, 7> make7Eigen()
  {
      return make_union(make3Eigen<T>(), make4Eigen<T>());
  }

  template <typename T> inline DirectionTableEigen<T, 9> make9Eigen()
  {
      return make_union(make3Eigen<T>(), make6Eigen<T>());
  }

  template <typename T> inline DirectionTableEigen<T, 10> make10Eigen()
  {
      return make_union(make4Eigen<T>(), make6Eigen<T>());
  }

  template <typename T> inline DirectionTableEigen<T, 13> make13Eigen()
  {
      return make_union(make7Eigen<T>(), make6Eigen<T>());
  }

  template<typename T, size_t N> class DirectionTableEigenHelper;

  template<typename T>
  class DirectionTableEigenHelper<T,3>
  {
  public:
      static DirectionTableEigen<T, 3> make() { return make3Eigen<T>(); }
  };

  template<typename T>
  class DirectionTableEigenHelper<T,4>
  {
  public:
      static DirectionTableEigen<T,4> make() { return make4Eigen<T>(); }
  };

  template<typename T>
  class DirectionTableEigenHelper<T,6>
  {
  public:
      static DirectionTableEigen<T, 6> make() { return make6Eigen<T>(); }
  };

  template<typename T>
  class  DirectionTableEigenHelper<T,7>
  {
  public:
      static DirectionTableEigen<T, 7> make() { return make7Eigen<T>(); }
  };

  template<typename T>
  class  DirectionTableEigenHelper<T,9>
  {
  public:
      static DirectionTableEigen<T, 9> make() { return make9Eigen<T>(); }
  };

  template<typename T>
  class  DirectionTableEigenHelper<T,10>
  {
  public:
      static DirectionTableEigen<T, 10> make() { return make10Eigen<T>(); }
  };

  template<typename T>
  class  DirectionTableEigenHelper<T,13>
  {
  public:
      static DirectionTableEigen<T, 13> make() { return make13Eigen<T>(); }
  };

  }// namespace geometry

// GEOMETRY_DIRECTION_TABLE_H
#endif
