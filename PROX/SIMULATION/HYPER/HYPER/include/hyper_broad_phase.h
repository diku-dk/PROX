#ifndef HYPER_BROAD_PHASE_H
#define HYPER_BROAD_PHASE_H

#include <hyper_engine.h>

#include <types/geometry_dop.h>
#include <overlap/geometry_overlap_dop_dop.h>

#include <cassert>    // needed for assert

#include <algorithm>  // needed for std::pair and std::sort

namespace hyper
{  

    template<typename MT>
    inline std::pair<Body<MT>*,Body<MT>*>  make_overlap( Body<MT> const * A, Body<MT> const * B )
    {          
      assert( A || !"make_overlap(): A object pointer was null");
      assert( B || !"make_overlap(): B object pointer was null");
      
      Body<MT> * bodyA = const_cast<Body<MT> *> ( A < B ? A : B );
      Body<MT> * bodyB = const_cast<Body<MT> *> ( A < B ? B : A );
      
      return std::make_pair( bodyA, bodyB );
    }

  /**
   *
   * @param efficiency    Upon return this argument gives the ratio of number of found
   *                      overlaps divided by the actual overlap tests done. A ratio of
   *                      close to 1 is optimal whereas a ratio close to zero is very bad.
   */
  template<typename MT>
  inline bool broad_phase(
                          Engine<MT> const & engine
                          , std::vector< std::pair<Body<MT>*,Body<MT>*>   >  & overlaps
                          , typename MT::real_type & efficiency
                          )
  {
    typedef          Engine<MT>                           engine_type;
    typedef typename engine_type::body_type               body_type;
    typedef typename engine_type::const_body_iterator     const_body_iterator;

    // Clean up any potential old left over information
    overlaps.clear();
    
    // Determine box overlaps
    const_body_iterator end   = engine.body_end();
    const_body_iterator iterA = engine.body_begin();
    
    size_t tests = 0u;
    
    for(;iterA!=end;)
    {     
      body_type const & A = (*iterA);

      const_body_iterator iterB = ++iterA;
      for(;iterB!=end;++iterB)
      {
        body_type const & B = (*iterB);

        tests++;

        if( ! geometry::overlap_dop_dop( A.m_tree.m_root, B.m_tree.m_root ) )
          continue;

        // Report that we have found an overlap between the bounding boxes of object A and B.
        overlaps.push_back( make_overlap(&A, &B) );
      }
    }

    // Lexiographic storting of overlaps
    std::sort( overlaps.begin(), overlaps.end() );

    efficiency = 1.0f*overlaps.size() / tests;
    
    // Return a status flag indicating whether we have seen an overlap or not
    return (overlaps.size()>0);
  }    

} //namespace hyper

// HYPER_BROAD_PHASE_H
#endif 
