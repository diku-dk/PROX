#include <geometry.h>
#include <eigenhelperall.h>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <vector>

using T = float;
class ContactInfo
{
public:
    EigenVector3<T> m_point;
    EigenVector3<T> m_normal;
    T m_distance;
};

class MyCallback : public geometry::ContactsCallback<T>
{
public:
    std::vector<ContactInfo> m_contacts;

public:
    void tempParenthesisOperatorImpl(const EigenVector3<T>& point,
                                     const EigenVector3<T>& normal,
                                     const T& distance)
    {
        ContactInfo info;

        info.m_point = point;
        info.m_normal = normal;
        info.m_distance = distance;

        m_contacts.push_back(info);
    }
};

BOOST_AUTO_TEST_SUITE(geometry);

BOOST_AUTO_TEST_CASE(contacts_obb_cylinder_test)
{
    {
        EigenVector3<T> const center = EigenVector3<T>(0.0, 0.0, 0.0);
        EigenVector3<T> const half_ext = EigenVector3<T>(1.0, 1.0, 1.0);
        const EigenQuaternion<T> q = EigenQuaternion<T>::Identity();

        EigenVector3<T> const center2 = EigenVector3<T>(0, 0, 0);
        EigenVector3<T> const axis = EigenVector3<T>(1, 0, 0);
        T const height = 1.0;
        T const radius = 1.0;

        geometry::OBBEigen<T> const& obb
            = geometry::make_obb<T>(center, q, half_ext);
        geometry::CylinderEigen<T> const& cyl
            = geometry::make_cylinder(radius, height, axis, center2);

        MyCallback mycallback;

        bool const test = geometry::contacts_obb_cylinder<T>(obb, cyl, 0.0,
                                                             mycallback, false);

        BOOST_CHECK(!test);
    }
}

BOOST_AUTO_TEST_SUITE_END();
