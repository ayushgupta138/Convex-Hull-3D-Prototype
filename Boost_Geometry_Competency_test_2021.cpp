
#include <iostream>
#include<list>
#include<vector>
#include<utility>
#include<initializer_list>
#include<random>
#include<chrono> 
#include<set>
#include<map>

#include<boost/geometry.hpp>

#define epsilon 1e-4

using namespace boost::geometry;

template<typename Point,std::size_t dim_count,std::size_t dim>
struct is_equal;

template
<
    typename Point
>
struct facet
{
    std::vector<Point> m_facet;
    inline void add_point(Point const& p)
    {
        m_facet.push_back(p);
    }
    template<typename Iterator>
    void insert_point(Point const& p, Point const &p1,Point const & p2)
    {
        for (auto it = boost::begin(m_facet); it != boost::end(m_facet) - 1; it++)
        {
            if ((is_equal<Point>::apply((*it), p1) && is_equal<Point>::apply(*(it+1), p2) || (is_equal<Point>::apply((*it), p2) && is_equal<Point>::apply(*(it + 1), p1))
            {
                m_facet.insert(it, p);
                return;
            }
        }
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

    void determine_point_order(Point const& P1, Point const& P2, Point& p1, Point& p2)
    {
        for (auto it = boost::begin(m_facet); it != boost::end(m_facet); it++)
        {
            if (is_equal<Point>::apply((*it), P1))
            {
                p1 = P2;
                p2 = P1;
                return;
            }
            else if (is_equal<Point>::apply((*it), P2))
            {
                p1 = P1;
                p2 = P2;
                return;
            }
        }
    }

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
        return m_facet[dim];
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
    // for testing purposes
    inline void print()
    {
        std::cout << wkt(m_point) << "\n";
    }
};

template
<
    typename Point
>
struct polyhedron
{
    std::vector<facet<Point>> m_face;
    std::vector<edge<Point>> m_edge;
    std::vector<vertex<Point>> m_vertex;
    
    template<std::size_t dim>
    inline facet<Point>* get_face()
    {
        BOOST_ASSERT((dim < m_face.size()));
        return &m_face[dim];
    }

    template<std::size_t dim>
    inline vertex<Point>* get_vertex()
    {
        BOOST_ASSERT((dim < m_vertex.size()));
        return &m_vertex[dim];
    }
    inline void print_poly_facet() // function for testing purposes
    {
        for (auto it = boost::begin(m_face); it != boost::end(m_face); it++)
        {
            it->print_facet();
        }
   }
    void print_edges()
    {
        for (auto it = boost::begin(m_edge); it != boost::end(m_edge); it++)
        {
            it->print_edge();
        }
    }
};

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
                m_unprocessed_points.push_back(unprocessed_point<Point>(*it));
            }
        }
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle(boost::begin(m_unprocessed_points), boost::end(m_unprocessed_points), std::default_random_engine(seed));
        construct_initial_conflict_graph();
        construct_hull();
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
    
   void construct_initial_conflict_graph()
   {
       initialize_conflict_graph();
       auto point_it = boost::begin(m_conflict_graph.m_point_list);
       for (auto it = boost::begin(m_unprocessed_points); it != boost::end(m_unprocessed_points); it++)
       {
           auto face_it = boost::begin(m_conflict_graph.m_facet_list);
           for (auto f_it = boost::begin(m_polyhedron.m_face); f_it != boost::end(m_polyhedron.m_face); f_it++)
           {
               std::size_t index1 = boost::numeric_cast<std::size_t>(f_it - boost::begin(m_polyhedron.m_face));
               std::size_t index2 = boost::numeric_cast<std::size_t>(it - boost::begin(m_unprocessed_points));
               if (is_visible<Point>(f_it->get<0>(), f_it->get<1>(), f_it->get<2>(), it->m_point) == above)
               {
                   point_it->second.push_back(&m_polyhedron.m_face[index1]);
                   face_it->second.push_back(&m_unprocessed_points[index2]);
               }
               face_it++;
           }
           point_it++;
       }
       //m_conflict_graph.print_graph();
   }

 
   void initialize_conflict_graph()
   {
       for (auto it = boost::begin(m_polyhedron.m_face); it != boost::end(m_polyhedron.m_face); it++)
       {
           std::size_t index = boost::numeric_cast<std::size_t>(it - boost::begin(m_polyhedron.m_face));
           m_conflict_graph.m_facet_list.push_back(std::make_pair(&m_polyhedron.m_face[index], std::vector<unprocessed_point<Point>*>()));
       }
       for (auto it = boost::begin(m_unprocessed_points); it != boost::end(m_unprocessed_points); it++)
       {
           std::size_t index = boost::numeric_cast<std::size_t>(it - boost::begin(m_unprocessed_points));
           m_conflict_graph.m_point_list.push_back(std::make_pair(&m_unprocessed_points[index], std::vector<facet<Point>*>()));
       }
   }

   void construct_hull()
   {
       for (auto it = boost::end(m_unprocessed_points)-1; it != boost::begin(m_unprocessed_points)+1; it--)  // there is need for forward iterators 
       {
           std::size_t index = boost::numeric_cast<std::size_t>(it - boost::begin(m_unprocessed_points));
           std::vector<edge<Point>*> edge_list,new_edge_list;
           std::vector<facet<Point>*> new_facet_list; 
           std::set<facet<Point>*> face_set;  // we can other use other efficient spatial index instead of set, but using set for simplicity 
           std::set<edge<Point>*> edge_set;
           std::set<vertex<Point>*> vertex_set;
           create_horizon_edge_list(edge_list,face_set,edge_set,vertex_set);
           create_vertex_set(edge_set, vertex_set);
       }
   }

   void update_hull(std::vector<edge<Point>*> const& edge_list,std::set<facet<Point>*> const & face_set,Point* const & p,std::size_t const & index)
   {
       
       m_polyhedron.m_vertex.push_back(vertex<Point>(p));
       facet<Point>* temp,temp1;
       for (auto it = boost::begin(edge_list); it != boost::end(edge_list) - 1; it++)
       {
           facet<Point>* not_visible;
           if (face_set.find((**it).m_facet1) == face_set.end())
               not_visible = (**it).m_facet1;
           else
               not_visible = (**it).m_facet2;
           if (is_visible(not_visible->get<0>(), not_visible->get<1>(), not_visible->get<2>(), *p) == on)
           {
               not_visible->insert_point(*p, *((**it).m_v1), *((**it).m_v2));
               if (it != boost::begin(edge_list))
               {
                   m_polyhedron.m_edge.push_back(edge<Point>(temp, not_visible, p, (**it).m_v1);
               }
               else
               {
                   temp1 = not_visible;
               }
               temp = not_visible;
           }
           else
           {
               Point p1, p2;
               not_visible->determine_point_order((**it).m_v1.m_point, (**it).m_v2.m_point, p1, p2);
               facet<Point> face = {p1,p2,*p,p1};
               m_polyhedron.m_face.push_back(face);
               if (it != boost::begin(edge_list))
               {
                   m_polyhedron.m_edge.push_back(edge<Point>(temp, &face, p, (**it).m_v1));
               }
               else
               {
                   temp1 = &face;
               }
               temp = &face;
           }
       }
       auto it = boost::begin(edge_list);
       m_polyhedron.m_edge.push_back(edge<Point>(temp, temp1, p, (**it).m_v1);
   }
   
   void update_conflict_graph()
   {

   }

   void create_vertex_set(std::set<edge<Point>*> const & edge_set,std::set<vertex<Point>*> & vertex_set)
   {
       std::set<vertex<Point>*> vertex_temp,temp;
       for (auto it = boost::begin(edge_set); it != boost::end(edge_set); it++)
       {
           vertex_temp.insert((**it).m_v1);
       }
       for (auto it = boost::begin(vertex_set); it != boost::end(vertex_set); it++)
       {
           if (vertex_temp.find((*it)) == vertex_temp.end())
               temp.insert((*it));
       }
       vertex_set = temp;
   }

   
   void create_horizon_edge_list(std::vector<edge<Point>*> & edge_list, std::set<facet<Point>*>& face_set,std::set<edge<Point>*> & edge_set,std::set<vertex<Point>*> & vertex_set)
   {
       std::vector<facet<Point>*> visible_faces = m_conflict_graph.m_point_list.back().second;
    
       for (auto it = boost::begin(visible_faces); it != boost::end(visible_faces); it++)
       {
           face_set.insert(*it);
       }
       for (auto it = boost::begin(m_polyhedron.m_edge); it != boost::end(m_polyhedron.m_edge); it++)
       {
           std::size_t count = 0;
           if (face_set.find((*it).m_facet1) != face_set.end())
           {
               count++;
           }
           if (face_set.find((*it).m_facet2) != face_set.end())
           {
               count++;
           }
           std::size_t index = boost::numeric_cast<size_t>(it - boost::begin(m_polyhedron.m_edge));
           if (count == 1)
           {
               edge_list.push_back(&m_polyhedron.m_edge[index]);
           }
           else if (count == 0)
           {
               edge_set.insert(&m_polyhedron.m_edge[index]);
               vertex_set.insert((*it).m_v1);
               vertex_set.insert((*it).m_v2);
           }
       }
       std::map<vertex<Point>*, std::size_t> vertex_map;
       for (auto it = boost::begin(edge_list); it != boost::end(edge_list); it++)
       {
           std::size_t index = boost::numeric_cast<std::size_t>(it - boost::begin(edge_list));
           vertex_map[(**it).m_v1] = index; // we can use more efficient spatial index, but using map for simplicity
       }
       edge<Point> *edge1 = edge_list.front();
       std::vector<edge<Point>*> edge_temp;
       edge_temp.push_back(edge1);
       edge<Point> *temp,*temp1;
       do
       {
           temp = edge_temp.back();
           temp1 = edge_list[vertex_map[(*temp).m_v2]];
           edge_temp.push_back(temp1);
       } while (temp1!=edge1);
       edge_list.clear();
       for (auto it = boost::begin(edge_temp); it != boost::end(edge_temp); it++)
       {
           edge_list.push_back(*it);
       }
       edge_list.push_back(edge1);
   }

   // should be declared as private, but public in this case for ease of testing
        polyhedron<Point> m_polyhedron;
        conflict_graph<Point> m_conflict_graph;
        std::vector<unprocessed_point<Point>> m_unprocessed_points;
};

int main()
{
    std::cout << "Hello World!\n";
    typedef model::d3::point_xyz<double> point3d;
    typedef model::multi_point<point3d> mulpoly;
    typedef model::ring<point3d> rng;
    //std::cout << is_visible(point3d(0, 0, 1), point3d(1, 0, 0), point3d(0, 1, 0), point3d(0.33, 0.33, 0.34)) << "\n";
    mulpoly mul;
    read_wkt("MULTIPOINT(0 0 0, 1 0 0 , 0 0 1 , 0 1 0 , -1 -1 0 )", mul);
   
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
