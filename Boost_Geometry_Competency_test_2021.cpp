
#include <iostream>
#include<list>
#include<utility>

#include<boost/geometry.hpp>

#define epsilon 1e-1

using namespace boost::geometry;

template<typename Point>
struct edge;

template
<
    typename Point
>
struct facet
{
    std::list<Point> m_facet;
    edge<Point>* m_edge;
    inline void add_point(Point const& p)
    {
        m_facet.push_back(p);
    }
    template<typename Iterator>
    inline void insert_point(Point const& p, Iterator const* it)
    {
        m_facet.insert(it, p);
    }
};

template
<
    typename Point
>
struct vertex
{
    Point m_vertex;
    edge<Point>* m_edge;
};

template
<
    typename Point
>
struct edge
{
    facet<Point>* m_left_facet, * m_right_facet;
    vertex<Point>* m_v1, * m_v2;
    edge<Point>* m_upper_left, * m_upper_right, * m_lower_left, * m_lower_right;
};

template<typename Point>
struct unprocessed_point
{
    Point m_point;
};

template
<
    typename Point
>
struct polyhedron
{
    std::vector<facet<Point>> m_facet;
    std::vector<edge<Point>> m_edge;
    std::vector<vertex<Point>> m_vertex;
};

template
<
    typename Point
>
struct conflict_graph
{
    std::list<std::pair<facet<Point>*,std::vector<unprocessed_point<Point>*>>> m_facet_list;
    std::list<std::pair<unprocessed_point<Point>*,std::vector<facet<Point>*>>> m_point_list;
};

enum location
{
    above = 0,
    below = 1,
    on = 2
};

inline void revert_enum(enum location & loc)
{
    if (loc == above)
        loc = below;
    else if (loc == below)
        loc = above;
}

template
<
    typename Point,
    typename volume_result = long double,
    bool ClockWise = false
>
inline location is_visible(Point const& P1, Point const& P2, Point const& P3,Point const &check)
{
    model::d3::point_xyz<volume_result> p1, p2, p3;
    p1.set<0>(boost::numeric_cast<volume_result>(get<0>(P1)) - boost::numeric_cast<volume_result>(get<0>(check)));
    p1.set<1>(boost::numeric_cast<volume_result>(get<1>(P1)) - boost::numeric_cast<volume_result>(get<1>(check)));
    p1.set<2>(boost::numeric_cast<volume_result>(get<2>(P1)) - boost::numeric_cast<volume_result>(get<2>(check)));

    p2.set<0>(boost::numeric_cast<volume_result>(get<0>(P2)) - boost::numeric_cast<volume_result>(get<0>(check)));
    p2.set<1>(boost::numeric_cast<volume_result>(get<1>(P2)) - boost::numeric_cast<volume_result>(get<1>(check)));
    p2.set<2>(boost::numeric_cast<volume_result>(get<2>(P2)) - boost::numeric_cast<volume_result>(get<2>(check)));

    p3.set<0>(boost::numeric_cast<volume_result>(get<0>(P3)) - boost::numeric_cast<volume_result>(get<0>(check)));
    p3.set<1>(boost::numeric_cast<volume_result>(get<1>(P3)) - boost::numeric_cast<volume_result>(get<1>(check)));
    p3.set<2>(boost::numeric_cast<volume_result>(get<2>(P3)) - boost::numeric_cast<volume_result>(get<2>(check)));

    volume_result result = volume_result(0);
    result = get<0>(p1) *
        (get<1>(p2) * get<2>(p3) -
            get<2>(p2) * get<1>(p3)) -
        get<0>(p2) *
        (get<1>(p1) * get<2>(p3) -
            get<2>(p1) * get<1>(p3)) +
        get<0>(p3) *
        (get<1>(p1) * get<2>(p2) -
            get<2>(p1) * get<1>(p2));

    location res;
    if (result < -epsilon)
        res = above;
    else if (result > epsilon)
        res = below;
    else
        res = on;
    if (ClockWise)
        revert_enum(res);
    return res;
}

template
<
    typename Point,
    typename collinearity_result = long double
>
inline bool is_collinear(Point const& P1, Point const& P2, Point const& P3)
{
    model::d2::point_xy<collinearity_result> p1, p2, p3;
    
    p1.set<0>(boost::numeric_cast<collinearity_result>(get<1>(P1)) - boost::numeric_cast<collinearity_result>(get<0>(P1)));
    p1.set<1>(boost::numeric_cast<collinearity_result>(get<2>(P1)) - boost::numeric_cast<collinearity_result>(get<0>(P1)));

    p2.set<0>(boost::numeric_cast<collinearity_result>(get<1>(P2)) - boost::numeric_cast<collinearity_result>(get<0>(P2)));
    p2.set<1>(boost::numeric_cast<collinearity_result>(get<2>(P2)) - boost::numeric_cast<collinearity_result>(get<0>(P2)));
    
    p3.set<0>(boost::numeric_cast<collinearity_result>(get<1>(P3)) - boost::numeric_cast<collinearity_result>(get<0>(P3)));
    p3.set<1>(boost::numeric_cast<collinearity_result>(get<2>(P3)) - boost::numeric_cast<collinearity_result>(get<0>(P3)));

    model::d3::point_xyz<collinearity_result> result;

    result.set<0>((boost::numeric_cast<collinearity_result>(get<0>(p2)) *
        boost::numeric_cast<collinearity_result>(get<1>(p3))) -
        boost::numeric_cast<collinearity_result>(get<1>(p2)) *
        boost::numeric_cast<collinearity_result>(get<0>(p3)));

    result.set<1>((boost::numeric_cast<collinearity_result>(get<1>(p1)) *
        boost::numeric_cast<collinearity_result>(get<0>(p3))) -
        boost::numeric_cast<collinearity_result>(get<0>(p1)) *
        boost::numeric_cast<collinearity_result>(get<1>(p3)));

    result.set<2>((boost::numeric_cast<collinearity_result>(get<0>(p1)) *
        boost::numeric_cast<collinearity_result>(get<1>(p2))) -
        boost::numeric_cast<collinearity_result>(get<1>(p1)) *
        boost::numeric_cast<collinearity_result>(get<0>(p2)));

    if (std::abs(get<0>(result)) < epsilon && std::abs(get<1>(result)) < epsilon && std::abs(get<2>(result) < epsilon))
        return true;
    else
        return false;
}

template
<
    typename Geometry,
    typename Point
>
inline bool find_point2D(Point const& p1, Point const& p2, Geometry const& geometry,Point &result)
{
    for (auto it = boost::begin(geometry) + 2; it != boost::end(geometry); it++)
    {
        if (!is_collinear(p1, p2, *it))
        {
            result = *it;
            return true;
        }
    }
    return false;
}

template
<
    typename Point
>
class convex_hull_3D
{
public:
    template<typename Geometry>
    inline void initialize_hull(Geometry const& geometry)
    {
        BOOST_CONCEPT_ASSERT((concepts::MultiPoint<Geometry>));
        BOOST_STATIC_ASSERT(geometry.size() > 3);
        typedef typename boost::range_value<Geometry>::type point_type;
        std::vector<point_type> initial_points;
        initial_points.push_back(*(boost::begin(geometry)));
        initial_points.push_back(*(boost::begin(geometry) + 1));

    }
private:
    class convex_hull_container
    {
        friend class convex_hull_3D;
        polyhedron<Point> m_polyhedron;
        conflict_graph<Point> m_conflict_graph;
    };
};

int main()
{
    std::cout << "Hello World!\n";
    typedef model::d3::point_xyz<double> point3d;
    typedef model::ring<point3d> rng;
    std::cout << is_visible(point3d(0, 0, 1), point3d(1, 0, 0), point3d(0, 1, 0), point3d(0.33, 0.33, 0.34)) << "\n";
}
