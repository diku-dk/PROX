#ifndef PROX_RIGID_BODY_H
#define PROX_RIGID_BODY_H

#include <narrow_object.h>   // Narrow phase collision detection interface support
#include <broad_object.h>    // Broad phase collision detection interface support

#include <tiny_matrix.h>
#include <tiny_quaternion_functions.h>        // for tiny::unit

#include <cassert>
#include <vector>

namespace prox
{

template <typename T>
requires(std::is_floating_point_v<T>)
class ForceCallback;

template < typename T >
requires(std::is_floating_point_v<T>)
class RigidBody : public narrow::Object<T>, public broad::Object<T>
{
public:
    using narrow_object = narrow::Object<T>;
    using broad_object = broad::Object<T>;

protected:
    bool m_fixed;         ///< if the body is fixed or not.
    bool         m_scripted;      ///< if the body is scripted or not.
    EigenVector3<T> m_r;             ///< The current position of the node in WCS.
    EigenQuaternion<T> m_Q;             ///< The current orientation of the node in WCS.
    EigenVector3<T> m_V;             ///< Linear velocity of CM in WCS.
    EigenVector3<T> m_W;             ///< Angular velocity around CM in WCS .
    EigenMatrix3<T> m_I_BF;          ///< Inertia tensor wrt CM in BF.
    T            m_mass;          ///< Mass.
    size_t       m_material_idx;  ///< The material index of the rigid body.
    size_t       m_idx;           ///< A body index.
    std::string  m_name;          ///< A name
    T            m_radius;        ///< A bounding radius of the rigid body.

    std::vector< ForceCallback<T>* > m_force_callbacks;

private:
    void copy(RigidBody const& body)
    {
      if(this==&body)
        return;

      this->m_fixed           = body.m_fixed;
      this->m_scripted        = body.m_scripted;
      this->m_r               = body.m_r;
      this->m_Q               = body.m_Q;
      this->m_V               = body.m_V;
      this->m_W               = body.m_W;
      this->m_I_BF            = body.m_I_BF;
      this->m_mass            = body.m_mass;
      this->m_material_idx    = body.m_material_idx;
      this->m_idx             = body.m_idx;
      this->m_name            = body.m_name;
      this->m_radius          = body.m_radius;
      this->m_force_callbacks = body.m_force_callbacks;
    }

  public:

    RigidBody()
    {
      this->clear();
    }

    virtual ~RigidBody(){	}

    RigidBody (RigidBody const & body)
    : narrow_object(body)
    , broad_object(body)
    {

      this->copy(body);
    }

    RigidBody & operator= (RigidBody const & body)
    {
      this->copy(body);
      return *this;
    }

  public:

    void clear()
    {
      this->m_fixed = false;
      this->m_scripted = false;
      this->m_r = {0, 0, 0};
      this->m_Q = EigenQuaternion<T>::Identity();
      this->m_I_BF = EigenMatrix3<T>({
          {1, 0, 0},
          {0, 1, 0},
          {0, 0, 1}
      });
      this->m_mass  = 1;
      this->m_V = {0, 0, 0};
      this->m_W = {0, 0, 0};
      this->m_material_idx = 0u;
      this->m_idx = 0u;
      this->m_name = "";
      this->m_radius = 1;
      this->m_force_callbacks.clear();
    }

    void set_idx( size_t const & idx ) { this->m_idx = idx; }
    size_t const & get_idx() const { return this->m_idx; }

    void set_material_idx( size_t const & idx ) { this->m_material_idx = idx; }
    size_t const & get_material_idx() const { return this->m_material_idx; }

    void set_fixed(bool const & fixed) { this->m_fixed = fixed; }
    bool is_fixed() const { return this->m_fixed; }

    void set_scripted(bool const & scripted) { this->m_scripted = scripted; }
    bool is_scripted() const { return this->m_scripted; }

    void set_position(EigenVector3<T> const& r) { this->m_r = r; }

    auto const& get_position() const { return this->m_r; }

    auto const& get_orientation() const { return this->m_Q; }

    void set_orientation(EigenQuaternion<T> const& Q) { this->m_Q = Q.normalized(); }

    void set_velocity(EigenVector3<T> const& V) { this->m_V = V; }

    auto const& get_velocity() const { return this->m_V; }

    void set_spin(EigenVector3<T> const& W) { this->m_W = W; }

    auto const& get_spin() const { return this->m_W; }

    void set_mass(T const& mass) { this->m_mass = mass; }

    T const& get_mass() const { return this->m_mass; }

    void set_inertia_bf(EigenMatrix3<T> const& I_BF) { this->m_I_BF = I_BF; }

    auto const& get_inertia_bf() const { return this->m_I_BF; }

    void set_name(std::string const & name) { this->m_name = name; }
    std::string const & get_name() const { return this->m_name; }

    void set_radius( T const & value )
    {
      assert( value > 0 || !"set_radius(): must be a positive value");
      this->m_radius = value;
    }

    /**
     * Broad phase interface implementation, must get information from narrow phase object interface about the geometry.
     */
    void get_box(T & mx,T & my,T & mz,T & Mx,T & My,T & Mz) const
    {
      assert( this->m_radius > 0 || !"get_box(): radius must be positive"                   );

      mx = this->m_r(0) - this->m_radius;
      my = this->m_r(1) - this->m_radius;
      mz = this->m_r(2) - this->m_radius;
      Mx = this->m_r(0) + this->m_radius;
      My = this->m_r(1) + this->m_radius;
      Mz = this->m_r(2) + this->m_radius;

      assert( mx < Mx || !"get_box(): min x must be less than max x");
      assert( my < My || !"get_box(): min y must be less than max y");
      assert( mz < Mz || !"get_box(): min z must be less than max z");
    }

    const auto& get_force_callbacks() const { return this->m_force_callbacks; }

    auto& get_force_callbacks() { return this->m_force_callbacks; }
};

} // namespace prox

// PROX_RIGID_BODY_H
#endif
