
#include <iostream>
#include<list>
#include<utility>
#include<initializer_list>
#include<random>
#include<chrono> 

#include<boost/geometry.hpp>

#define epsilon 1e-4

using namespace boost::geometry;

template
<
    typename Point
>
struct facet
{
    std::list<Point> m_facet;
    inline void add_point(Point const& p)
    {
        m_facet.push_back(p);
    }
    template<typename Iterator>
    inline void insert_point(Point const& p, Iterator const* it)
    {
        m_facet.insert(it, p);
    }
    inline facet(std::initializer_list<Point> l)
    {
        for (auto it = l.begin(); it != l.end(); it++)
        {
            m_facet.push_back(*it);
        }
    }
    inline facet()
    {}
    inline void print_facet()   // function for testing purposes
    {
        for (auto it = boost::begin(m_facet); it != boost::end(m_facet); it++)
        {
            std::cout << wkt(*it) << " ";
        }
        std::cout << "\n";
    }
    template<std::size_t dim>
    inline Point get()
    {
        BOOST_ASSERT((dim < m_facet.size()));
        auto it = boost::begin(m_facet);
        std::size_t count = dim;
        while (count > 0)
        {
            count--;
            it++;
        }
        return *it;
    }
};

template
<
    typename Point
>
struct vertex
{
    Point m_vertex;
    vertex(Point const& p)
    {
        m_vertex.set<0>(get<0>(p));
        m_vertex.set<1>(get<1>(p));
        m_vertex.set<2>(get<2>(p));
    }
    inline void get_pointer()
    {
        return &m_vertex;
    }
};

template
<
    typename Point
>
struct edge
{
    facet<Point>* m_facet1, * m_facet2;
    vertex<Point>* m_v1, * m_v2;
    inline edge(facet<Point>* f1, facet<Point>* f2, vertex<Point>* v1, vertex<Point>* v2)
    {
        m_facet1 = f1;
        m_facet2 = f2;
        m_v1 = v1;
        m_v2 = v2;
    }
};

template<typename Point>
struct unprocessed_point
{
    Point m_point;
    inline unprocessed_point(Point const& p)
    {
        m_point.set<0>(get<0>(p));
        m_point.set<1>(get<1>(p));
        m_point.set<2>(get<2>(p));
    }
};

template
<
    typename Point
>
struct polyhedron
{
    std::list<facet<Point>> m_face;
    std::list<edge<Point>> m_edge;
    std::list<vertex<Point>> m_vertex;
    
    template<std::size_t dim>
    inline facet<Point>* get_face()
    {
        BOOST_ASSERT((dim < m_face.size()));
        auto it = boost::begin(m_face);
        std::size_t count = dim;
        while (count > 0)
        {
            count--;
            it++;
        }
        facet<Point> f = *it;
        return &f;
    }

    template<std::size_t dim>
    inline vertex<Point>* get_vertex()
    {
        BOOST_ASSERT((dim < m_vertex.size()));
        auto it = boost::begin(m_vertex);
        std::size_t count = dim;
        while (count > 0)
        {
            count--;
            it++;
        }
        vertex<Point> v = *it;
        return &v;
    }
    inline void print_poly_facet() // function for testing purposes
    {
        for (auto it = boost::begin(m_face); it != boost::end(m_face); it++)
            it->print_facet();
   }
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
    
    p1.set<0>(boost::numeric_cast<collinearity_result>(get<0>(P2)) - boost::numeric_cast<collinearity_result>(get<0>(P1)));
    p1.set<1>(boost::numeric_cast<collinearity_result>(get<0>(P3)) - boost::numeric_cast<collinearity_result>(get<0>(P1)));

    p2.set<0>(boost::numeric_cast<collinearity_result>(get<1>(P2)) - boost::numeric_cast<collinearity_result>(get<1>(P1)));
    p2.set<1>(boost::numeric_cast<collinearity_result>(get<1>(P3)) - boost::numeric_cast<collinearity_result>(get<1>(P1)));
    
    p3.set<0>(boost::numeric_cast<collinearity_result>(get<2>(P2)) - boost::numeric_cast<collinearity_result>(get<2>(P1)));
    p3.set<1>(boost::numeric_cast<collinearity_result>(get<2>(P3)) - boost::numeric_cast<collinearity_result>(get<2>(P1)));

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
    typename Point,
    std::size_t dim_count = dimension<Point>::value,
    std::size_t dim = 0
>
struct is_equal
{
    static inline bool apply(Point const&p1,Point const &p2)
    {
        return ((get<dim>(p1) == get<dim>(p2)) && is_equal<Point, dim_count, dim + 1>::apply(p1, p2));
    }
};
template
<
    typename Point,
    std::size_t dim_count
>
struct is_equal<Point, dim_count, dim_count>
{
    static inline bool apply(Point const& p1, Point const& p2)
    {
        return true;
    }
};

template
<
    typename Geometry,
    typename Point
>
inline bool find_point3D(Point const& p1, Point const& p2, Point const& p3, Geometry const& geometry, Point & result)
{
    for (auto it = boost::begin(geometry) + 2; it != boost::end(geometry); it++)
    {
        if (! is_equal<Point>::apply(*it,p3) && is_visible(p1, p2, p3, *it) != on)
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
    template
        <
         typename Geometry
        >
    inline void initialize_hull(Geometry const& geometry)
    {
        BOOST_CONCEPT_ASSERT((concepts::MultiPoint<Geometry>));
        BOOST_ASSERT((geometry.size() > 3));
        typedef typename boost::range_value<Geometry>::type point_type;
        std::vector<point_type> initial_points;
        initial_points.push_back(*(boost::begin(geometry)));
        initial_points.push_back(*(boost::begin(geometry) + 1));
        point_type result;
        bool res = find_point2D(initial_points[0], initial_points[1], geometry, result);
        BOOST_ASSERT(res);
        initial_points.push_back(result);
        res = find_point3D(initial_points[0], initial_points[1], initial_points[2], geometry, result);
        BOOST_ASSERT(res);
        initial_points.push_back(result);
        construct_initial_polyhedron(initial_points);
        std::vector<unprocessed_point<Point>> shuffled_unprocessed_points;
        for (auto it = boost::begin(geometry) + 2; it != boost:end(geometry); it++)
        {
            if (!is_equal<Point>::apply(initial_points[2], *it) && !is_equal<Point>::apply(initial_points[3], *it))
            {
                shuffled_unprocessed_points.push_back(unprocessed_point<Point>(*it));
            }
        }
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle(boost::begin(shuffled_unprocessed_points), boost::end(shuffled_unprocessed_points), std::default_random_engine(seed));
    }
   template
        <
        typename Point
        >
    inline void construct_initial_polyhedron(std::vector<Point> const& initials)
    {
        std::vector<Point> ccw_order = { initials[0],initials[1],initials[2] };
        if (is_visible(ccw_order[0], ccw_order[1], ccw_order[2], initials[3]) == above)
        {
            std::reverse(boost::begin(ccw_order), boost::end(ccw_order));
        }
        facet<Point> face;
        
        // create four faces of the initial hull

        face = { ccw_order[0],ccw_order[1],ccw_order[2],ccw_order[0] }; // face 1   (1 2 3 1)
        m_polyhedron.m_face.push_back(face);
        face = { ccw_order[0],initials[3],ccw_order[1],ccw_order[0] };  // face 2   (1 4 2 1)
        m_polyhedron.m_face.push_back(face); 
        face = { ccw_order[1],initials[3],ccw_order[2],ccw_order[1] };  // face 3   (2 4 3 2)
        m_polyhedron.m_face.push_back(face);
        face = { ccw_order[0],ccw_order[2],initials[3],ccw_order[0] };  // face 4   (1 3 4 1)
        m_polyhedron.m_face.push_back(face);

        // store initial 4 vertices of the hull
        
        m_polyhedron.m_vertex.push_back(vertex<Point>(ccw_order[0]));
        m_polyhedron.m_vertex.push_back(vertex<Point>(ccw_order[1]));
        m_polyhedron.m_vertex.push_back(vertex<Point>(ccw_order[2]));
        m_polyhedron.m_vertex.push_back(vertex<Point>(initials[3]));

        // store initial 6 edges of the hull
        
        m_polyhedron.m_edge.push_back(edge<Point>(m_polyhedron.get_face<1>(), m_polyhedron.get_face<0>(), m_polyhedron.get_vertex<0>(), m_polyhedron.get_vertex<1>()));
        m_polyhedron.m_edge.push_back(edge<Point>(m_polyhedron.get_face<3>(), m_polyhedron.get_face<0>(), m_polyhedron.get_vertex<0>(), m_polyhedron.get_vertex<2>()));
        m_polyhedron.m_edge.push_back(edge<Point>(m_polyhedron.get_face<1>(), m_polyhedron.get_face<3>(), m_polyhedron.get_vertex<0>(), m_polyhedron.get_vertex<3>()));
        m_polyhedron.m_edge.push_back(edge<Point>(m_polyhedron.get_face<0>(), m_polyhedron.get_face<2>(), m_polyhedron.get_vertex<1>(), m_polyhedron.get_vertex<2>()));
        m_polyhedron.m_edge.push_back(edge<Point>(m_polyhedron.get_face<1>(), m_polyhedron.get_face<2>(), m_polyhedron.get_vertex<1>(), m_polyhedron.get_vertex<3>()));
        m_polyhedron.m_edge.push_back(edge<Point>(m_polyhedron.get_face<3>(), m_polyhedron.get_face<2>(), m_polyhedron.get_vertex<2>(), m_polyhedron.get_vertex<3>()));
        //m_polyhedron.print_poly_facet();
    }
    
   template
       <
       typename Point
       >
   inline void construct_initial_conflict_graph(std::vector<unprocessed_point<Point>> const& points)
   {

   }

   // should be declared as private, but public in this case for ease of testing
        polyhedron<Point> m_polyhedron;
        conflict_graph<Point> m_conflict_graph;
};

int main()
{
    std::cout << "Hello World!\n";
    typedef model::d3::point_xyz<double> point3d;
    typedef model::multi_point<point3d> mulpoly;
    typedef model::ring<point3d> rng;
    //std::cout << is_visible(point3d(0, 0, 1), point3d(1, 0, 0), point3d(0, 1, 0), point3d(0.33, 0.33, 0.34)) << "\n";
    mulpoly mul;
    read_wkt("MULTIPOINT(0 0 0, 1 0 0,0 0 1,0 1 0,2 2 2)", mul);
    convex_hull_3D<point3d> pt;
    pt.initialize_hull(mul);
    for (auto it = boost::begin(pt.m_polyhedron.m_face); it != boost::end(pt.m_polyhedron.m_face); it++)
    {
        it->print_facet();
        if (is_visible(it->get<0>(), it->get<1>(), it->get<2>(), point3d(-1,-1,0)) == above)
            std::cout << "Visible\n";
        else
            std::cout << "Not visible\n";
     }
}
