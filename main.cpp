#include "convex_hull_3D.hpp"
#include<iostream>
#include<chrono>
#include<random>

using namespace boost::geometry;
using namespace std::chrono;

int main()
{
    typedef model::d3::point_xyz<double> point3d;
    typedef model::multi_point<point3d> mulpoly;
    typedef model::ring<point3d> rng;
    mulpoly mul;
    // for 100 random points input
    int x, y, z;
    for (int i = 0; i < 100; i++)
    {
        std::cin >> x >> y >> z;
        append(mul, point3d(x, y, z));
    }
    std::cout << "\n\n";
    //read_wkt("MULTIPOINT(0 0 0,1 1 1,1 2 3,0 0 2,0 2 0,0 2 2,1 1 4)", mul);
    polyhedron<point3d> result;
    auto start = high_resolution_clock::now();
    result = convex_hull3D(mul);
    auto stop = high_resolution_clock::now();
    result.print_poly_facet();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Time taken by 3D convex hull function : " << duration.count() << " ms\n";
}
