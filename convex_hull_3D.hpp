
#include <iostream>
#include<list>
#include<vector>
#include<utility>
#include<initializer_list>
#include<random>
#include<chrono> 
#include<set>
#include<map>
#include<unordered_set>

#include<boost/geometry.hpp>

#define epsilon 1e-4

using namespace boost::geometry;
using namespace std::chrono;

template<typename Point,std::size_t dim_count,std::size_t dim>
struct is_equal;

// wrapper for facet of a polyhedron

template
<
    typename Point
>
struct facet
{
    std::vector<Point*> m_facet_ptr;
    inline void add_point(Point const& p) 
    {
        Point* ptr = new Point(p);
        m_facet_ptr.push_back(ptr);
    }

    inline void remove_point()    
    {
        m_facet_ptr.pop_back();
    }

    void insert_point(Point const& p, Point const &p1,Point const & p2)
    {
        std::size_t index = 0;
        for (auto it = boost::begin(m_facet_ptr); it != boost::end(m_facet_ptr) - 1; it++)
        {
            if ((is_equal<Point>::apply((**it), p1) && is_equal<Point>::apply(**(it+1), p2)) || (is_equal<Point>::apply((**it), p2) && is_equal<Point>::apply(**(it + 1), p1)))
            {
                Point* ptr = new Point(p);
                m_facet_ptr.insert(boost::begin(m_facet_ptr) + index + 1, ptr);
                return;
            }
            index++;
        }
    }
    inline facet(std::initializer_list<Point> l)
    {
        for (auto it = l.begin(); it != l.end(); it++)
        {
            Point* ptr = new Point(*it);
            m_facet_ptr.push_back(ptr);
        }
    }
    inline facet()
    {}

    void determine_point_order(Point const& P1, Point const& P2, Point& p1, Point& p2)
    {
        for (auto it = boost::begin(m_facet_ptr); it != boost::end(m_facet_ptr); it++)
        {
            if (is_equal<Point>::apply((**it), P1) && is_equal<Point>::apply(**(it+1),P2))
            {
                p1 = P2;
                p2 = P1;
                return;
            }
            else if (is_equal<Point>::apply((**it), P2) && is_equal<Point>::apply(**(it + 1), P1))
            {
                p1 = P1;
                p2 = P2;
                return;
            }
        }
    }

    inline void print_facet()   // function for testing purposes
    {
        for (auto it = boost::begin(m_facet_ptr); it != boost::end(m_facet_ptr); it++)
        {
            std::cout << wkt(**it) << " ";
        }
        std::cout << "\n";
    }
    template<std::size_t dim>
    inline Point get()
    {
        BOOST_ASSERT((dim < m_facet_ptr.size()));
        return *m_facet_ptr[dim];
    }
};

// wrapper for vertex of a polyhedron

template
<
    typename Point
>
struct vertex
{
    Point m_vertex;
    inline vertex(Point const& p)
    {
        m_vertex.set<0>(get<0>(p));
        m_vertex.set<1>(get<1>(p));
        m_vertex.set<2>(get<2>(p));
    }
    inline vertex(vertex<Point> const& v)
    {
        m_vertex = v.m_vertex;
    }
    inline void get_pointer()
    {
        return &m_vertex;
    }
    inline void print()
    {
        std::cout << wkt(m_vertex) << "\n";
    }

};

// wrapper that store edge information of a polyhedron, edge's end points,incidence relation,etc

template
<
    typename Point
>
struct edge
{
    facet<Point>* m_facet1, * m_facet2;
    vertex<Point>* m_v1, * m_v2;
   
    inline edge()
    {}
    
    inline edge(facet<Point>* f1, facet<Point>* f2, vertex<Point>* v1, vertex<Point>* v2)
    {
        m_facet1 = f1;
        m_facet2 = f2;
        m_v1 = v1;
        m_v2 = v2;
    }
    void print_edge()
    {
        std::cout << "Adjacent faces are :\n";
        m_facet1->print_facet();
        m_facet2->print_facet();
        std::cout << "Adjacent vertices are:\n";
        std::cout << wkt(m_v1->m_vertex) << "\n";
        std::cout << wkt(m_v2->m_vertex) << "\n";
    }
};

// struct to store unprocessed points(points yet to be added to the hull)

template<typename Point>
struct unprocessed_point
{
    Point m_point;
    
    inline unprocessed_point()
    {}
    
    inline unprocessed_point(std::initializer_list<Point> l)
    {
        for (auto it = l.begin(); it != l.end(); it++)
        {
            m_point.push_back(*it);
        }
    }
    inline unprocessed_point(Point const& p)
    {
        m_point.set<0>(get<0>(p));
        m_point.set<1>(get<1>(p));
        m_point.set<2>(get<2>(p));
    }
    inline unprocessed_point(unprocessed_point<Point> const & u_point)
    {
        m_point = u_point.m_point;
    }
    // for testing purposes
    inline void print()
    {
        std::cout << wkt(m_point) << "\n";
    }
};

// struct to represent a polyhedron, contains vertices,faces, edges and incidence relation

template
<
    typename Point
>
struct polyhedron
{
    std::vector<facet<Point>*> m_face_ptr;
    std::vector<edge<Point>*> m_edge_ptr;
    std::vector<vertex<Point>*> m_vertex_ptr;

    template<std::size_t dim>
    inline facet<Point>* get_face()
    {
        BOOST_ASSERT((dim < m_face_ptr.size()));
        return m_face_ptr[dim];
    }

    template<std::size_t dim>
    inline vertex<Point>* get_vertex()
    {
        BOOST_ASSERT((dim < m_vertex_ptr.size()));
        return m_vertex_ptr[dim];
    }
    inline void add_face(facet<Point> const& face)
    {
        facet<Point>* ptr = new facet<Point>(face);
        m_face_ptr.push_back(ptr);
    }
    inline void add_edge(edge<Point> const& edges)
    {
        edge<Point>* ptr = new edge<Point>(edges);
        m_edge_ptr.push_back(ptr);
    }
    inline void add_vertex(vertex<Point> const& vertices)
    {
        vertex<Point>* ptr = new vertex<Point>(vertices);
        m_vertex_ptr.push_back(ptr);
    }
    inline void remove_face()
    {
        m_face_ptr.pop_back();
    }
    inline void remove_edge()
    {
        m_edge_ptr.pop_back();
    }
    inline void remove_vertex()
    {
        m_vertex_ptr.pop_back();
    }
    
    // functions for testing purposes
    inline void print_poly_facet()
    {
        for (auto it = boost::begin(m_face_ptr); it != boost::end(m_face_ptr); it++)
        {
            (*it)->print_facet();
        }
   }
    void print_edges()
    {
        for (auto it = boost::begin(m_edge_ptr); it != boost::end(m_edge_ptr); it++)
        {
            (*it)->print_edge();
        }
    }
    void print_vertex()
    {
        for (auto it = boost::begin(m_vertex_ptr); it != boost::end(m_vertex_ptr); it++)
        {
            (*it)->print();
        }
    }
};

// An important feature of randomized incremental algorithm, conflict graph. Helps in achieving speedup

template
<
    typename Point
>
struct conflict_graph
{
    std::vector<std::pair<facet<Point>*,std::vector<unprocessed_point<Point>*>>> m_facet_list;
    std::vector<std::pair<unprocessed_point<Point>*,std::vector<facet<Point>*>>> m_point_list;
    
    inline conflict_graph()
    {}

    // for testing purposes
    inline void print_graph()
    {
        std::cout << "Printing facet list:\n";
        for (auto it = boost::begin(m_facet_list); it != boost::end(m_facet_list); it++)
        {
            (*(*it).first).print_facet();
            std::cout << " : \n";
            for (auto itr = boost::begin((*it).second); itr != boost::end((*it).second); itr++)
            {
                (*(*itr)).print();
            }
            
        }
        std::cout << "Printing point list:\n";
        for (auto it = boost::begin(m_point_list); it != boost::end(m_point_list); it++)
        {
            (*(*it).first).print();
            std::cout << " : \n";
            for (auto itr = boost::begin((*it).second); itr != boost::end((*it).second); itr++)
            {
                (*(*itr)).print_facet();
            }
        }
    }
};

// comparator definition for vertex structure, using lamda expressions instead of old operator() standard was efficient, but using operator() for simplicity

template
<
    typename Point
>
struct comp
{
    bool operator()(vertex<Point> const& v1, vertex<Point> const& v2) const
    {
        if (get<0>(v1.m_vertex) != get<0>(v2.m_vertex))
            return get<0>(v1.m_vertex) < get<0>(v2.m_vertex);
        else if (get<1>(v1.m_vertex) != get<1>(v2.m_vertex))
            return get<1>(v1.m_vertex) < get<1>(v2.m_vertex);
        return get<2>(v1.m_vertex) < get<2>(v2.m_vertex);
    }
};

// enum for denoting orientation of a point with respect to a planar face

enum location
{
    above = 0,
    below = 1,
    on = 2
};

//for enum reversal

inline void revert_enum(enum location& loc)
{
    if (loc == above)
        loc = below;
    else if (loc == below)
        loc = above;
}

// for determining equality of two Points
template
<
    typename Point,
    std::size_t dim_count = dimension<Point>::value,
    std::size_t dim = 0
>
struct is_equal
{
    static inline bool apply(Point const& p1, Point const& p2)
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

// function to obtain orientation of point w.r.t a plane

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

// determine collinearity of three points

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

// find 3rd non collinear point for initial hull 

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

// find 4th non-coplanar point for initial hull

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



// for point removal

template
<
    typename T
>
void erase_points(std::vector<T>& point_list, T points)
{
    auto del = boost::begin(point_list);
    for (auto it = boost::begin(point_list); it != boost::end(point_list); it++)
    {
        if ((*it) == points)
        {
            del = it;
            break;
        }
    }
    point_list.erase(del);
}

template
<
    typename Point
>
class convex_hull_3D
{
public:
    typedef std::vector<std::pair<facet<Point>*, std::pair<facet<Point>*, facet<Point>*>>> face_list;
    
    inline void add_unprocessed_point(unprocessed_point<Point> u_point)
    {
        m_unprocessed_points.push_back(u_point);
        unprocessed_point<Point>* ptr = new unprocessed_point<Point>(u_point);
        m_unprocessed_points_ptr.push_back(ptr);
    }

    inline void remove_unprocessed_point()
    {
        m_unprocessed_points.pop_back();
        m_unprocessed_points_ptr.pop_back();
    }

    // construct initial hull(tetradedron), throw exception if hull cannot be constructed
    
    template
        <
         typename Geometry
        >
    void initialize_hull(Geometry const& geometry)
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
        for (auto it = boost::begin(geometry) + 2; it != boost::end(geometry); it++)
        {
            if (!is_equal<Point>::apply(initial_points[2], *it) && !is_equal<Point>::apply(initial_points[3], *it))
            {
                add_unprocessed_point(unprocessed_point<Point>(*it));
            }
        }
        std::vector<std::pair<unprocessed_point<Point>, unprocessed_point<Point>*>> point_shuffle;
        for (auto it = boost::begin(m_unprocessed_points); it != boost::end(m_unprocessed_points); it++)
        {
            std::size_t index = boost::numeric_cast<std::size_t>(it - boost::begin(m_unprocessed_points));
            point_shuffle.push_back(std::make_pair(*it, m_unprocessed_points_ptr[index]));
        }
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle(boost::begin(point_shuffle), boost::end(point_shuffle), std::default_random_engine(seed));
        m_unprocessed_points.clear();
        m_unprocessed_points_ptr.clear();
        for (auto it = boost::begin(point_shuffle); it != boost::end(point_shuffle); it++)
        {
            m_unprocessed_points.push_back(it->first);
            m_unprocessed_points_ptr.push_back(it->second);
        }
        construct_initial_conflict_graph();
    }

    // add vertices, edges and faces to initial hull

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
        m_polyhedron.add_face(face);
        face = { ccw_order[0],initials[3],ccw_order[1],ccw_order[0] };  // face 2   (1 4 2 1)
        m_polyhedron.add_face(face);
        face = { ccw_order[1],initials[3],ccw_order[2],ccw_order[1] };  // face 3   (2 4 3 2)
        m_polyhedron.add_face(face);
        face = { ccw_order[0],ccw_order[2],initials[3],ccw_order[0] };  // face 4   (1 3 4 1)
        m_polyhedron.add_face(face);

        // store initial 4 vertices of the hull
        
        m_polyhedron.add_vertex(vertex<Point>(ccw_order[0]));
        m_polyhedron.add_vertex(vertex<Point>(ccw_order[1]));
        m_polyhedron.add_vertex(vertex<Point>(ccw_order[2]));
        m_polyhedron.add_vertex(vertex<Point>(initials[3]));

        // store initial 6 edges of the hull
        
        m_polyhedron.add_edge(edge<Point>(m_polyhedron.get_face<1>(), m_polyhedron.get_face<0>(), m_polyhedron.get_vertex<0>(), m_polyhedron.get_vertex<1>()));
        m_polyhedron.add_edge(edge<Point>(m_polyhedron.get_face<3>(), m_polyhedron.get_face<0>(), m_polyhedron.get_vertex<0>(), m_polyhedron.get_vertex<2>()));
        m_polyhedron.add_edge(edge<Point>(m_polyhedron.get_face<1>(), m_polyhedron.get_face<3>(), m_polyhedron.get_vertex<0>(), m_polyhedron.get_vertex<3>()));
        m_polyhedron.add_edge(edge<Point>(m_polyhedron.get_face<0>(), m_polyhedron.get_face<2>(), m_polyhedron.get_vertex<1>(), m_polyhedron.get_vertex<2>()));
        m_polyhedron.add_edge(edge<Point>(m_polyhedron.get_face<1>(), m_polyhedron.get_face<2>(), m_polyhedron.get_vertex<1>(), m_polyhedron.get_vertex<3>()));
        m_polyhedron.add_edge(edge<Point>(m_polyhedron.get_face<3>(), m_polyhedron.get_face<2>(), m_polyhedron.get_vertex<2>(), m_polyhedron.get_vertex<3>()));
    }
    
   // initialization and construction of conflict graph

   void construct_initial_conflict_graph()
   {
       initialize_conflict_graph();
       auto point_it = boost::begin(m_conflict_graph.m_point_list);
       for (auto it = boost::begin(m_unprocessed_points); it != boost::end(m_unprocessed_points); it++)
       {
           auto face_it = boost::begin(m_conflict_graph.m_facet_list);
           for (auto f_it = boost::begin(m_polyhedron.m_face_ptr); f_it != boost::end(m_polyhedron.m_face_ptr); f_it++)
           {
               std::size_t index1 = boost::numeric_cast<std::size_t>(f_it - boost::begin(m_polyhedron.m_face_ptr));
               std::size_t index2 = boost::numeric_cast<std::size_t>(it - boost::begin(m_unprocessed_points));
               if (is_visible<Point>((*f_it)->get<0>(), (*f_it)->get<1>(), (*f_it)->get<2>(), it->m_point) == above)
               {
                   point_it->second.push_back(*f_it);
                   face_it->second.push_back(m_unprocessed_points_ptr[index2]);
               }
               face_it++;
           }
           point_it++;
       }
   }

   // initialize conflict graph

   void initialize_conflict_graph()
   {
       for (auto it = boost::begin(m_polyhedron.m_face_ptr); it != boost::end(m_polyhedron.m_face_ptr); it++)
       {
           std::vector<unprocessed_point<Point>*> poin;
           m_conflict_graph.m_facet_list.push_back(std::make_pair(*it, poin));
       }
       m_conflict_graph.m_point_list.clear();
       for (auto it = boost::begin(m_unprocessed_points); it != boost::end(m_unprocessed_points); it++)
       {
           std::size_t index = boost::numeric_cast<std::size_t>(it - boost::begin(m_unprocessed_points));
           std::vector<facet<Point>*> face;
           m_conflict_graph.m_point_list.push_back(std::make_pair(m_unprocessed_points_ptr[index], face));
       }
   }

   // get hidden face given edge

   facet<Point>* get_invisible_face(edge<Point> * edge,std::set<facet<Point>*> const &face_set)
   {
       facet<Point>* invisible;
       if (face_set.find((*edge).m_facet1) == face_set.end())
       {
           invisible = (*edge).m_facet1;
       }
       else
       {
           invisible = (*edge).m_facet2;
       }
       return invisible;
   }
   
   // reorder edge list to maintain continuity of points of same face

   void order_edge_list(std::vector<edge<Point>*> & edge_list,std::set<facet<Point>*> const & face_set)
   {
       auto it1 = boost::begin(edge_list);
       auto it3 = boost::begin(edge_list);
       std::size_t index1=0, index2=edge_list.size()-1;
       while (get_invisible_face(*it1, face_set) == get_invisible_face(*it3, face_set))
       {
           it1++;
           index1++;
           if (it1 == boost::end(edge_list))
               return;
       }
       std::vector<edge<Point>*> edge_temp;
       for (auto it = boost::begin(edge_list)+index1; it != boost::end(edge_list)-1; it++)
       {
           edge_temp.push_back(*it);
       }
       for (auto it = boost::begin(edge_list); it != boost::begin(edge_list)+index1; it++)
       {
           edge_temp.push_back(*it);
       }
       edge_temp.push_back(*(boost::begin(edge_list)+index1));
       edge_list = edge_temp;
   }

   // processed each unprocessed point and update the hull and conflict graph

   void construct_hull()
   {
       for (auto it = boost::rbegin(m_unprocessed_points); it != boost::rend(m_unprocessed_points); it++)   
       {
           std::vector<edge<Point>*> edge_list;
           face_list new_facet_list;
           std::set<facet<Point>*> face_set;  // we can other use other efficient spatial index instead of set, but using set for simplicity 
           std::set<edge<Point>*> edge_set;
           std::set<vertex<Point>*> vertex_set;
           create_horizon_edge_list(edge_list,face_set,edge_set,vertex_set);
           if (face_set.size() == 0)
           {
               update_conflict_graph(new_facet_list);
               continue;
           }
           order_edge_list(edge_list, face_set);
           create_vertex_set(edge_set, vertex_set);
           vertex<Point> v1(it->m_point);
           m_polyhedron.add_vertex(v1);
           update_hull(edge_list, face_set, m_polyhedron.m_vertex_ptr.back(), new_facet_list,vertex_set,edge_set);
           update_conflict_graph(new_facet_list);
           clean_up(face_set, edge_set, vertex_set);
       }
   }

   // final clean up and removal of faces,edges and vertices after point addition

   void clean_up(std::set<facet<Point>*> const& face_set, std::set<edge<Point>*> const& edge_set, std::set<vertex<Point>*> const& vertex_set)
   {
       std::set<facet<Point>*> face_set_new;  
       std::set<edge<Point>*> edge_set_new;
       std::set<vertex<Point>*> vertex_set_new;
       for (auto it = boost::begin(m_polyhedron.m_face_ptr); it != boost::end(m_polyhedron.m_face_ptr); it++)
       {
           if (face_set.find(*it) == face_set.end())
           {
               face_set_new.insert(*it);
           }
       }
       for (auto it = boost::begin(m_polyhedron.m_edge_ptr); it != boost::end(m_polyhedron.m_edge_ptr); it++)
       {
           if (edge_set.find(*it) == edge_set.end())
           {
               edge_set_new.insert(*it);
           }
       }

       for (auto it = boost::begin(m_polyhedron.m_vertex_ptr); it != boost::end(m_polyhedron.m_vertex_ptr); it++)
       {
           if (vertex_set.find(*it) == vertex_set.end())
           {
               vertex_set_new.insert(*it);
           }
       }
       m_polyhedron.m_face_ptr.clear();
       m_polyhedron.m_edge_ptr.clear();
       m_polyhedron.m_vertex_ptr.clear();

       for (auto it = boost::begin(face_set_new); it != boost::end(face_set_new); it++)
       {
           m_polyhedron.m_face_ptr.push_back(*it);
       }

       for (auto it = boost::begin(edge_set_new); it != boost::end(edge_set_new); it++)
       {
           m_polyhedron.m_edge_ptr.push_back(*it);
       }

       for (auto it = boost::begin(vertex_set_new); it != boost::end(vertex_set_new); it++)
       {
           m_polyhedron.m_vertex_ptr.push_back(*it);
       }
   }

   // update hull

   void update_hull(std::vector<edge<Point>*> const& edge_list,std::set<facet<Point>*> const & face_set,vertex<Point>* const p,face_list & new_facet_list,std::set<vertex<Point>*> &vertex_set,std::set<edge<Point>*> & edge_set)
   {
       facet<Point>* temp,*not_visible_temp,*temp1;
       temp = temp1 = not_visible_temp = nullptr;
       for (auto it = boost::begin(edge_list); it != boost::end(edge_list) - 1; it++)
       {
           facet<Point> *not_visible,*visible;
          if (face_set.find((**it).m_facet1) == face_set.end())
           {
               not_visible = (**it).m_facet1;
               visible = (**it).m_facet2;
           }
           else
           {
               not_visible = (**it).m_facet2;
               visible = (**it).m_facet1;
           }
           if (is_visible(not_visible->get<0>(), not_visible->get<1>(), not_visible->get<2>(), p->m_vertex) == on)
           {
               std::set<vertex<Point>*> vertex_temp;
               std::set<vertex<Point>,comp<Point>>  vertex_del;
               auto itr = it;
               for (itr = it; itr != boost::end(edge_list) - 1; itr++)
               {
                   if (face_set.find((**itr).m_facet1) == face_set.end())
                   {
                       not_visible_temp = (**itr).m_facet1;
                   }
                   else
                   {
                       not_visible_temp = (**itr).m_facet2;
                   }
                   if (not_visible_temp != not_visible)
                   {
                       break;
                   }
                   else
                   {
                       edge_set.insert(*itr);
                       if (vertex_temp.find((**itr).m_v1) == vertex_temp.end())
                       {
                           vertex_temp.insert((**itr).m_v1);
                       }
                       else
                       {
                           vertex_temp.erase((**itr).m_v1);
                           vertex_del.insert(*((**itr).m_v1));
                           vertex_set.insert((**itr).m_v1);
                       }
                       if (vertex_temp.find((**itr).m_v2) == vertex_temp.end())
                       {
                           vertex_temp.insert((**itr).m_v2);
                       }
                       else
                       {
                           vertex_temp.erase((**itr).m_v2);
                           vertex_del.insert(*((**itr).m_v2));
                           vertex_set.insert((**itr).m_v2);
                       }
                   }
               }
               itr--;
               vertex<Point>* v1, * v2;
               v1 = v2 = NULL;
               for (auto i : vertex_temp)
               {
                   if (v1 == NULL)
                       v1 = i;
                   else if (v2 == NULL)
                       v2 = i;
                   else
                       break;
               }
               preprocess(vertex_del, not_visible);
               not_visible->insert_point(p->m_vertex, v1->m_vertex, v2->m_vertex);
               if (it != boost::begin(edge_list))
               {
                   m_polyhedron.add_edge(edge<Point>(temp, not_visible, p, (**it).m_v1));
               }
               else
               {
                   temp1 = not_visible;
               }
               temp = not_visible;
               it = itr;

           }
           else
           {
               Point p1, p2;
               not_visible->determine_point_order((**it).m_v1->m_vertex, (**it).m_v2->m_vertex, p1, p2);
               facet<Point> face = {p1,p2,p->m_vertex,p1};
               m_polyhedron.add_face(face);
               update_edge_facets(*it, visible, m_polyhedron.m_face_ptr.back());
               if (it != boost::begin(edge_list))
               {
                   m_polyhedron.add_edge(edge<Point>(temp, m_polyhedron.m_face_ptr.back(), p, (**it).m_v1));
               }
               else
               {
                   temp1 = m_polyhedron.m_face_ptr.back();
               }
               temp = m_polyhedron.m_face_ptr.back();
               new_facet_list.push_back(std::make_pair(m_polyhedron.m_face_ptr.back(), std::make_pair(visible, not_visible)));
           }
       }
       auto it = boost::begin(edge_list);
       m_polyhedron.add_edge(edge<Point>(temp, temp1, p, (**it).m_v1));
   }
 
   // get alternate face

   inline void update_edge_facets(edge<Point>* e, facet<Point>* old, facet<Point>* new_facet)
   {
       if (e->m_facet1 == old)
       {
           e->m_facet1 = new_facet;
       }
       else
       {
           e->m_facet2 = new_facet;
       }
   }

   // update conflict graph

   void update_conflict_graph(face_list const & new_facet_list)
   {
       if (m_conflict_graph.m_point_list.size() < 2)
           return;
       std::map<unprocessed_point<Point>*, std::size_t> point_map;
       std::map<facet<Point>*, std::size_t> face_map;
       for (auto it = boost::begin(m_conflict_graph.m_point_list); it != boost::end(m_conflict_graph.m_point_list); it++)
       {
           std::size_t index = boost::numeric_cast<std::size_t>(it - boost::begin(m_conflict_graph.m_point_list));
           point_map[(*it).first] = index;
       }
       for (auto it = boost::begin(m_conflict_graph.m_facet_list); it != boost::end(m_conflict_graph.m_facet_list); it++)
       {
           std::size_t index = boost::numeric_cast<std::size_t>(it - boost::begin(m_conflict_graph.m_facet_list));
           face_map[(*it).first] = index;
       }
       std::set<std::size_t> face_removed;
       for (auto it = boost::begin(m_conflict_graph.m_point_list.back().second); it != boost::end(m_conflict_graph.m_point_list.back().second); it++)
       {
           face_removed.insert(face_map[*it]);
           erase_points(m_conflict_graph.m_facet_list[face_map[*it]].second, m_conflict_graph.m_point_list.back().first);
       }
       for (auto it = boost::begin(new_facet_list); it != boost::end(new_facet_list); it++)
       {
           std::vector<unprocessed_point<Point>*> union_result;
           facet<Point>* face = (*it).first;
           create_union(m_conflict_graph.m_facet_list[face_map[(*it).second.first]].second, m_conflict_graph.m_facet_list[face_map[(*it).second.second]].second, union_result,face);
           m_conflict_graph.m_facet_list.push_back(std::make_pair((*it).first, union_result));
           for (auto itr = boost::begin(union_result); itr != boost::end(union_result); itr++)
           {
               m_conflict_graph.m_point_list[point_map[*itr]].second.push_back(face);
           }
       }
       
       for (auto it = boost::rbegin(face_removed); it != boost::rend(face_removed); it++)
       {
           std::set<unprocessed_point<Point>*> processed_points;
           for (auto itr = boost::begin(m_conflict_graph.m_facet_list[*it].second); itr != boost::end(m_conflict_graph.m_facet_list[*it].second); itr++)
           {
               if (processed_points.find(*itr) == processed_points.end())
               {
                   processed_points.insert(*itr);
                   erase_points(m_conflict_graph.m_point_list[point_map[*itr]].second, m_conflict_graph.m_facet_list[*it].first);
               }
           }
           m_conflict_graph.m_facet_list.erase(boost::begin(m_conflict_graph.m_facet_list) + (*it));
       }
       m_conflict_graph.m_point_list.pop_back();
   }

   // preprocess face and remove vertices which are no longer part of hull

   void preprocess(std::set<vertex<Point>,comp<Point>> const& vertex_del, facet<Point>* face)
   {
       facet<Point> face1;
       for (auto it = boost::begin(face->m_facet_ptr); it != boost::end(face->m_facet_ptr); it++)
       {
           if (vertex_del.find(vertex<Point>(**it)) == vertex_del.end())
           {
               face1.add_point(**it);
           }
       }
       if (!is_equal<Point>::apply(**(boost::begin(face1.m_facet_ptr)),**(boost::rbegin(face1.m_facet_ptr))))
       {
           Point p = **(boost::begin(face1.m_facet_ptr));
           face1.add_point(p);
       }
       (*face) = face1;
   }

   // perform union and return new conflict graph values

   void create_union(std::vector<unprocessed_point<Point>*> const& v1, std::vector<unprocessed_point<Point>*> const& v2, std::vector<unprocessed_point<Point>*>& result,facet<Point>* face)
   {
       std::set<unprocessed_point<Point>*> u_set;
       for (auto it = boost::begin(v1); it != boost::end(v1); it++)
       {
           u_set.insert(*it);
       }
       for (auto it = boost::begin(v2); it != boost::end(v2); it++)
       {
           u_set.insert(*it);
       }
       result.clear();
       for (auto it = boost::begin(u_set); it != boost::end(u_set); it++)
       {
           if (is_visible(face->get<0>(), face->get<1>(), face->get<2>(), (*it)->m_point) == above)
           {
               result.push_back(*it);
           }
       }
   }

   // determine vertices to be deleted 
  
   void create_vertex_set(std::set<edge<Point>*> const & edge_set,std::set<vertex<Point>*> & vertex_set)
   {
       for (auto it = boost::begin(edge_set); it != boost::end(edge_set); it++)
       {
           if (vertex_set.find((**it).m_v1) != vertex_set.end())
           {
               vertex_set.erase((**it).m_v1);
           }
           if (vertex_set.find((**it).m_v2) != vertex_set.end())
           {
               vertex_set.erase((**it).m_v2);
           }
       }
   }

   // create horizon list

   void create_horizon_edge_list(std::vector<edge<Point>*> & edge_list, std::set<facet<Point>*>& face_set,std::set<edge<Point>*> & edge_set,std::set<vertex<Point>*> & vertex_set)
   {  
       std::vector<facet<Point>*> visible_faces = m_conflict_graph.m_point_list.back().second;
       if (visible_faces.size() == 0)
           return;
       for (auto it = boost::begin(visible_faces); it != boost::end(visible_faces); it++)
       {
           face_set.insert(*it);
       }
       for (auto it = boost::begin(m_polyhedron.m_edge_ptr); it != boost::end(m_polyhedron.m_edge_ptr); it++)
       {
           std::size_t count = 0;
           if (face_set.find((**it).m_facet1) != face_set.end())
           {
               count++;
           }
           if (face_set.find((**it).m_facet2) != face_set.end())
           {
               count++;
           }
           std::size_t index = boost::numeric_cast<size_t>(it - boost::begin(m_polyhedron.m_edge_ptr));
           if (count == 1)
           {
               edge_list.push_back(m_polyhedron.m_edge_ptr[index]);
           }
           else if (count == 2)
           {
               edge_set.insert(m_polyhedron.m_edge_ptr[index]);
               vertex_set.insert((**it).m_v1);
               vertex_set.insert((**it).m_v2);
           }
       }
       std::map<vertex<Point>*, std::vector<std::size_t>> vertex_map1,vertex_map2;
       for (auto it = boost::begin(edge_list); it != boost::end(edge_list); it++)
       {
           std::size_t index = boost::numeric_cast<std::size_t>(it - boost::begin(edge_list));
           vertex_map1[(**it).m_v1].push_back(index); // we can use more efficient spatial index, but using map for simplicity
           vertex_map2[(**it).m_v2].push_back(index);
       }
       edge<Point> *edge1 = edge_list.front();
       std::vector<edge<Point>*> edge_temp;
       edge_temp.push_back(edge1);
       edge<Point> *temp,*temp1,*temp3;
       temp1 = temp = temp3 = NULL;
       vertex<Point>* temp2;
       temp2 = edge1->m_v2;
       do
       {
           temp = edge_temp.back();
           std::set<edge<Point>*> vec;
           if (vertex_map1[temp2].size())
           {
               for (auto it : vertex_map1[temp2])
                   vec.insert(edge_list[it]);
           }
           if (vertex_map2[temp2].size())
           {
               for (auto it : vertex_map2[temp2])
                   vec.insert(edge_list[it]);
           }
           for (auto it = boost::begin(vec); it != boost::end(vec); it++)
           {
               if ((*it) != temp)
               {
                   temp1 = *it;
                   break;
               }
           }
           temp2 = return_alternate(*temp1, temp2);
           if (temp->m_v2 != temp1->m_v1)
               exchange_vertex(temp1);
           edge_temp.push_back(temp1);
       } while (temp1!=edge1);
       edge_list.clear();
       for (auto it = boost::begin(edge_temp); it != boost::end(edge_temp); it++)
       {
           edge_list.push_back(*it);
       }
   }

   // utility function

   inline void exchange_vertex(edge<Point>* e)
   {
       vertex<Point>* temp;
       temp = e->m_v1;
       e->m_v1 = e->m_v2;
       e->m_v2 = temp;
   }

   // utility function

   inline vertex<Point>* return_alternate(edge<Point> const& edge, vertex<Point>* const& vertex)
   {
       if (edge.m_v1 == vertex)
           return edge.m_v2;
       return edge.m_v1;
   }

   // should be declared as private, but public in this case for ease of testing
        polyhedron<Point> m_polyhedron;
        conflict_graph<Point> m_conflict_graph;
        std::vector<unprocessed_point<Point>> m_unprocessed_points;
        std::vector<unprocessed_point<Point>*> m_unprocessed_points_ptr;
};

// free function for hull construction

template
<
    typename Geometry1
>
polyhedron<model::d3::point_xyz<double>> convex_hull3D(Geometry1 const& input)
{
    convex_hull_3D<model::d3::point_xyz<double>> hull;
    hull.initialize_hull(input);
    hull.construct_hull();
    return hull.m_polyhedron;
}
