#include <set>
#include <cmath>
#include <gmsh.h>
#include <vector>


int main(int argc, char **argv)
{
    gmsh::initialize();

    gmsh::model::add("tor");

    const double xc = 0;
    const double yc = 0;
    const double zc = 0;

    const double r_in = 0.1;
    const double r_out = 0.15;
    const double R = 0.45;

    gmsh::model::occ::addTorus(xc, yc, zc, R, r_out, 1);
    gmsh::model::occ::addTorus(xc, yc, zc, R, r_in, 2);

    using Pair = std::vector<std::pair<int, int>>;
    Pair out;
    std::vector<Pair> outmap;

    gmsh::model::occ::cut({{3, 1}}, {{3, 2}}, out, outmap, 3);

    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeFactor", 0.25);

    gmsh::model::mesh::generate(3);

    gmsh::write("tor.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();

    return 0;
}
