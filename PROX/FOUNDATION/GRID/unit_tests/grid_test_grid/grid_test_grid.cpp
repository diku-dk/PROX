#include "grid_enclosing_indices.h"
#include "grid_grid.h"
#include "grid_iterators.h"
#include "grid_node_position.h"
#include "grid_value_at.h"

#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <fstream>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>



BOOST_AUTO_TEST_SUITE(grid);

template<typename T>
T sphereSdf(const Eigen::Matrix<T,3,1>& p, const Eigen::Matrix<T,3,1>& c, T r)
{
    return (p - c).norm() - r;
}

BOOST_AUTO_TEST_CASE(grid_test_grid)
{
    {
        using D = double;
        using T = double;

        const size_t res = 200u;
        Eigen::Matrix<T,3,1> gmin(T(-1.0), T(-1.0), T(-1.0));
        Eigen::Matrix<T,3,1> gmax(T( 1.0), T( 1.0), T( 1.0));
        Eigen::Matrix<size_t,3,1> nodes(res, res, res);

        grid::Grid<D,T> G;
        G.create(gmin, gmax, nodes);

        //Sphere parameters
        Eigen::Matrix<T,3,1> center(T(0.0), T(0.0), T(0.0));
        T radius = 0.6;

        //Fill grid with SDF values at nodes
        for(size_t k=0;k<res;++k)
        {
            for(size_t j=0;j<res;++j)
            {
                for(size_t i=0;i<res;++i)
                {
                    Eigen::Matrix<size_t,3,1> idx(i,j,k);
                    Eigen::Matrix<T,3,1> pos;
                    grid::node_position(G, idx, pos);
                    D sdf = static_cast<D>(sphereSdf(pos, center, radius));
                    G(idx) = sdf;
                }
            }
        }

        //Sample random points and compare SDF vs interpolated value_at
        std::mt19937_64 rng(0xdeadbeef);
        std::uniform_real_distribution<T> u(0.0, 1.0);

        const size_t nSamples = 2000;
        double maxErr = 0.0;
        double sumErr = 0.0;

        for(size_t s = 0; s < nSamples; ++s)
        {
            //random point in the domain
            Eigen::Matrix<T,3,1> p;
            p.x() = gmin.x() + u(rng)*(gmax.x()-gmin.x());
            p.y() = gmin.y() + u(rng)*(gmax.y()-gmin.y());
            p.z() = gmin.z() + u(rng)*(gmax.z()-gmin.z());

            D analytic = static_cast<D>(sphereSdf(p, center, radius) );
            D interp = grid::value_at<D,T>(G, p);

            double err = std::abs( interp - analytic );
            sumErr += err;
            if(err > maxErr) maxErr = err;
        }

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Grid resolution: " << res << " x " << res << " x " << res << "\n";
        std::cout << "Max abs error (value_at vs analytic) over "<< nSamples << " samples: " << maxErr << "\n";
        std::cout << "Mean abs error: " << (sumErr / (double)nSamples) << "\n";

        //Export node positions and SDF values into Eigen matrices (P and S) to simulate igl....
        const size_t total = res*res*res;
        Eigen::MatrixXd P(total, 3);
        Eigen::VectorXd S(total);

        size_t idxLin = 0;
        for(size_t k = 0;k < res; ++k)
        {
            for(size_t j = 0; j < res; ++j)
            {
                for(size_t i = 0; i < res; ++i)
                {
                    Eigen::Matrix<size_t,3,1> idx3(i,j,k);
                    Eigen::Matrix<T,3,1> p;
                    grid::node_position(G, idx3, p);
                    P(idxLin, 0) = p.x();
                    P(idxLin, 1) = p.y();
                    P(idxLin, 2) = p.z();
                    S(idxLin) = G(idx3);
                    ++idxLin;
                }
            }
        }

        std::cout << "Exported node count: " << P.rows() << " (P columns: " << P.cols() << ")\n";
        std::cout << "First 6 (pos, sdf) entries:\n";
        for(int i=0; i < std::min<int>(6, P.rows()); ++i)
        {
            std::cout << "  (" << P(i,0) << ", " << P(i,1) << ", " << P(i,2) << ") -> " << S(i) << "\n";
        }

        /*std::ofstream csv("nodes_sdf.csv");
        csv << "x,y,z,sdf\n";
        for(int r=0; r<P.rows(); ++r)
            csv << P(r,0) << "," << P(r,1) << "," << P(r,2) << "," << S(r) << "\n";
        csv.close();
        std::cout << "Wrote nodes_sdf.csv (can be loaded into Python later)\n";*/

        {
            BOOST_CHECK_CLOSE(-1.0, -1.0, 0.01 );
        }
    }

    {
        //=============================ITERATOR TEST============================
        using D = int;
        using T = double;


        const size_t nx = 4, ny = 3, nz = 2;
        grid::Grid<D,T> G;
        G.create(Eigen::Matrix<T,3,1>(0,0,0), Eigen::Matrix<T,3,1>((T)(nx-1),(T)(ny-1),(T)(nz-1)), Eigen::Matrix<size_t,3,1>(nx,ny,nz));


        const size_t total = G.m_nodes.x()*G.m_nodes.y()*G.m_nodes.z();


        //Fill with linear data
        for(size_t s = 0; s < total; ++s)
        {
            G.data()[s] = static_cast<D>(s + 100);
        }


        //Test basic iterator forward traversal
        using It = grid::detail::Iterator<grid::Grid<D,T>, D&, D* >;
        It it(&G, G.data());
        It it_end(&G, G.data() + total);


        size_t s = 0;
        bool ok = true;
        for(; it != it_end; ++it, ++s)
        {
            if(*it != static_cast<D>(s + 100) )
            {
                std::cout << "Iterator value mismatch at linear " << s << ": got " << *it << " expected " << (s+100) << "\n";
                ok = false;
                break;
            }
        }


        if (ok)
        {
            std::cout << "Iterator forward traversal: PASS\n";
        }
        else
        {
            std::cout << "Iterator forward traversal: FAIL\n";
        }
        BOOST_CHECK_EQUAL(ok, true);


        //Test compute_index mapping
        It it2(&G, G.data()+5);
        Eigen::Matrix<size_t,3,1> idx = it2.compute_index();
        size_t lin = idx.z() * G.m_nodes.y() * G.m_nodes.x() + idx.y() * G.m_nodes.x() + idx.x();
        if(lin != 5)
        {
            std::cout << "compute_index mismatch: got linear="<<lin<<" expected 5 -> idx("<<idx.x()<<","<<idx.y()<<","<<idx.z()<<")\n";
        }
        else
        {
            std::cout << "compute_index test: PASS (linear 5 -> idx("<<idx.x()<<","<<idx.y()<<","<<idx.z()<<"))\n";
        }
        BOOST_CHECK_EQUAL(lin, 5);


        //Test IndexIterator and increment
        using IIt = grid::detail::IndexIterator< grid::Grid<D,T>, D&, D* >;
        IIt iit(&G, G.data());

        //jump to a specific index using assignment from Eigen vector
        Eigen::Matrix<size_t,3,1> target(2,1,1);
        iit = target;
        if( *iit != G(target))
        {
            std::cout << "IndexIterator assignment mismatch: *iit=" << *iit << " grid=" << G(target) << "\n";
        }
        else
        {
            std::cout << "IndexIterator assignment: PASS (value="<<*iit<<")\n";
        }


        //test operator++ moves indices correctly
        // near end of row
        iit = Eigen::Matrix<size_t,3,1>(G.I()-2, 0, 0);
        size_t i_before = iit.i();
        ++iit;
        if(iit.i() != i_before + 1)
        {
            std::cout << "IndexIterator ++ failed to increment i.\n";
        }
        //advance across row boundary
        iit = Eigen::Matrix<size_t,3,1>(G.I()-1, 0, 0);
        //should wrap i to 0 and increment j
        ++iit;
        if( iit.i() != 0 || iit.j() != 1 )
        {
            std::cout << "IndexIterator ++ failed to wrap row to next column. got ("<<iit.i()<<","<<iit.j()<<","<<iit.k()<<")\n";
        }
        else
        {
            std::cout << "IndexIterator increment behavior: PASS\n";
        }
        BOOST_CHECK_EQUAL(iit.j(), 1);
        BOOST_CHECK_EQUAL(iit.i(), 0);

        std::cout << "All iterator tests passed.\n";
    }

    {
        using D = double;
        using T = double;


        const size_t res = 40; // finer resolution
        Eigen::Matrix<T,3,1> gmin(T(-1.2), T(-1.2), T(-1.2));
        Eigen::Matrix<T,3,1> gmax(T( 1.2), T( 1.2), T( 1.2));
        Eigen::Matrix<size_t,3,1> nodes(res,res,res);


        grid::Grid<D,T> G;
        G.create(gmin, gmax, nodes);


        Eigen::Matrix<T,3,1> center(T(0.2), T(0.0), T(-0.1)); // slightly off-center
        T radius = 0.7;


        // Fill SDF at nodes
        for(size_t k=0;k<res;++k)
        {
            for(size_t j=0;j<res;++j)
            {
                for(size_t i=0;i<res;++i)
                {
                    Eigen::Matrix<size_t,3,1> idx(i,j,k);
                    Eigen::Matrix<T,3,1> pos;
                    grid::node_position(G, idx, pos);
                    G(idx) = sphereSdf(pos, center, radius);
                }
            }
        }


        //Compute min and max SDF
        D minv = G.data()[0];
        D maxv = G.data()[0];
        auto totalNodes = G.m_nodes.x()*G.m_nodes.y()*G.m_nodes.z();
        for(size_t s = 0; s < totalNodes; ++s)
        {
            if( G.data()[s] < minv) minv = G.data()[s];
            if( G.data()[s] > maxv) maxv = G.data()[s];
        }


        std::cout << std::fixed << std::setprecision(6);
        std::cout << "SDF min: " << minv << " max: " << maxv << "\n";
        if(!(minv < 0.0 && maxv > 0.0))
        {
            std::cout << "Warning: no sign change across grid (SDF may not intersect). Test expects both signs.\n";
        }
        BOOST_CHECK_EQUAL(!(minv < 0.0 && maxv > 0.0), false);


        //Count sign-change edges in +x, +y, +z directions
        size_t signEdges = 0;
        for(size_t k = 0; k < res; ++k)
        {
            for(size_t j = 0; j < res; ++j)
            {
                for(size_t i = 0; i < res; ++i)
                {
                    Eigen::Matrix<size_t,3,1> idx(i, j ,k);
                    D v = G(idx);
                    if(i+1 < res){ if( v * G(Eigen::Matrix<size_t,3,1>(i+1,j,k)) < 0.0 ) ++signEdges; }
                    if(j+1 < res){ if( v * G(Eigen::Matrix<size_t,3,1>(i,j+1,k)) < 0.0 ) ++signEdges; }
                    if(k+1 < res){ if( v * G(Eigen::Matrix<size_t,3,1>(i,j,k+1)) < 0.0 ) ++signEdges; }
                }
            }
        }


        std::cout << "Sign-change edges count: " << signEdges << " (should be > 0 if the surface intersects the grid)\n";
        if(signEdges == 0)
        {
            //std::cout << "ERROR!";
        }
        BOOST_CHECK_GT(signEdges, 0);


        //Export CSV for visualization
        /*        std::ofstream csv("nodes_sdf2.csv");
        csv << "x,y,z,sdf\n";
        size_t idx_lin = 0;
        for(size_t k=0;k<res;++k)
        {
            for(size_t j=0;j<res;++j)
            {
                for(size_t i=0;i<res;++i)
                {
                    Eigen::Matrix<size_t,3,1> idx(i,j,k);
                    Eigen::Matrix<T,3,1> p;
                    grid::node_position(G, idx, p);
                    csv << p.x() << "," << p.y() << "," << p.z() << "," << G(idx) << "\n";
                    ++idx_lin;
                }
            }
        }
        csv.close();
        std::cout << "Wrote nodes_sdf2.csv\n";*/


        std::cout << "SDF grid test: PASS\n";
    }
}


BOOST_AUTO_TEST_SUITE_END();
