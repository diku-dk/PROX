#include <big_types.h>
#include <big_generate_PD.h>
#include <big_generate_PSD.h>
#include <big_generate_random.h>

#include <big_min_map_newton.h>
#include <big_projected_gauss_seidel.h>

#include <iostream>

template<typename matrix_type,typename vector_type>
inline void test(matrix_type const & A, vector_type  & x, vector_type const & b, vector_type const & y, size_t * cnt_status)
{
  typedef typename matrix_type::value_type real_type;
  typedef typename matrix_type::size_type  size_type;
  
  size_type max_iter = 100;
  real_type tol_abs  = boost::numeric_cast<real_type>(1e-6);
  real_type tol_rel  = boost::numeric_cast<real_type>(0.000000001);
  real_type tol_stag = boost::numeric_cast<real_type>(0.000000001);
  size_t status      = 0;
  size_type iter     = 0;
  real_type err      = boost::numeric_cast<real_type>(0.0);
  real_type alpha    = boost::numeric_cast<real_type>(0.0001);
  real_type beta     = boost::numeric_cast<real_type>(0.5);

  
  // 2012-09-30 Kenny: Here we warmstart Newton method with PGS solution, try out commenting this
  //                   On a note we could use Mueller's and Nuttapongs MG-LCP solver for warmstarting
  //
  big::projected_gauss_seidel(
                      A, b, x
                      , 10u      // 2012-09-30 Kenny: I just picked some max iter number for PGS... have no idea how it affects the number of required Newton iterations, but even a few PGS seems to be worthwhile
                      , tol_abs
                      , tol_rel
                      , tol_stag
                      , status
                      , iter
                      , err
                      );
  
  big::min_map_newton(
                      A, b, x
                      , max_iter
                      , tol_abs
                      , tol_rel
                      , tol_stag
                      , status
                      , iter
                      , err
                      , alpha
                      , beta
                      );
  
  cnt_status[status]++;
  
  std::cout << "error = "
            << err
            << " found in "
            << iter
            << " Newton iterations "
            << " with final solver status= "
            << big::get_error_message(status)
            << std::endl;
}

  
int main(int argc, char **argv)
{
  typedef ublas::compressed_matrix<double> matrix_type;
  typedef ublas::vector<double>            vector_type;
  typedef vector_type::size_type           size_type;
  
  size_t cnt_status[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  
  size_type N = 100;  // 2012-09-30 Kenny: right now we use SVD as subsolver in Newton method -- higher dimenions than this kind of removes the interactive aspect of this demo... when we replace SVD with a proper iterative method we should be able to scale to much higher dimensions.
  
  matrix_type A;
  vector_type x;
  vector_type b;
  vector_type y;
  
  A.resize(N,N,false);
  x.resize(N,false);
  b.resize(N,false);
  y.resize(N,false);
  
  for(size_type tst=0u;tst<100u;++tst)
  {
    // 2012-09-30 Kenny: This is just some random generated problem -- it is not very fluid like, its dense sym and PD, fluid problems are sparse sym banded and PSD!!!
    big::fast_generate_PD(N,A);
    big::generate_random( N, y);
    
    b.assign(-y);
    big::generate_random( N, y);
    x.clear();
    
    test(A, x, b, y, &cnt_status[0]);
  }
  
  for(size_t i=0;i<9;++i)
  {
    std::cout << cnt_status[i]
              << "\t:\t"
              << big::get_error_message( i )
              << std::endl;
  }
  std::cout << std::endl;

	return 0;
}
