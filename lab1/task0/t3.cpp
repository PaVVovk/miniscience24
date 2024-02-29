#include <set>
#include <gmsh.h>
#include <vector>


int main(int argc, char **argv)
{
    gmsh::initialize();

    gmsh::model::add("t3");

    double lc = 0.01;

    const double xc = 0;
    const double yc = 0;
    const double zc = 0;

    gmsh::model::geo::addPoint(xc, yc, zc, lc, 1);

    const double r = 0.1;

    gmsh::model::geo::addPoint(xc - r, yc, zc, lc, 2);
    gmsh::model::geo::addPoint(xc, yc + r, zc, lc, 3);
    gmsh::model::geo::addPoint(xc + r, yc, zc, lc, 4);
    gmsh::model::geo::addPoint(xc, yc - r, zc, lc, 5);

    const double h = 0.25;


    for (int i = 1; i < 4; i++){
        gmsh::model::geo::addCircleArc(i+1, 1, i+2, i);
    }
    gmsh::model::geo::addCircleArc(5, 1, 2, 4);

    gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);

    std::vector<std::pair<int, int> > ov2;
    gmsh::model::geo::extrude({{2, 1}}, 0, 0, h, ov2);

    gmsh::model::geo::synchronize();

    gmsh::model::mesh::generate(3);

    gmsh::write("t3.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();

    return 0;
}
