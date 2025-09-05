#include <simulators_prox.h>
#include <simulators_prox_data.h>

#include <cassert>

namespace simulators
{

size_t ProxEngine::get_number_of_contacts()
{
    assert(m_data || !"ProxEngine::get_number_of_contacts(): Internal error, data null pointer");

    return m_data->m_contacts.size();
}

void ProxEngine::get_contact_position(size_t const& contact_number, float& x, float& y, float& z)
{
    assert(m_data || !"ProxEngine::get_contact_position(): Internal error, data null pointer");

    assert(contact_number < m_data->m_contacts.size()
           || !"ProxEngine::get_contact_position(): internal error: contact index out of bounds");

    auto [px, py, pz] = m_data->m_contacts[contact_number].position;
    x = px;
    y = py;
    z = pz;
}

void ProxEngine::get_contact_normal(size_t const& contact_number, float& x, float& y, float& z)
{
    assert(m_data || !"ProxEngine::get_contact_normal(): Internal error, data null pointer");

    assert(contact_number < m_data->m_contacts.size()
           || !"ProxEngine::get_contact_normal(): internal error: contact index out of bounds");

    auto [px, py, pz] = m_data->m_contacts[contact_number].normal;
    x = px;
    y = py;
    z = pz;
}

void ProxEngine::get_contact_depth(size_t const& contact_number, float& depth)
{
    assert(m_data || !"ProxEngine::get_contact_depth(): Internal error, data null pointer");

    assert(contact_number < m_data->m_contacts.size()
           || !"ProxEngine::get_contact_depth(): internal error: contact index out of bounds");

    depth = m_data->m_contacts[contact_number].depth;
}

}// namespace simulators
