#include <set>
#include <cmath>
#include <gmsh.h>
#include <vector>


int main(int argc, char **argv)
{
    gmsh::initialize();

    gmsh::model::add("cylinder");

    const double xc = 0;
    const double yc = 0;
    const double zc = 0;

    const double R = 7e-3;
    const double r = 5e-4;
    const double h = 0.4;

    const double lc = 5 * 1e-4;
    const double lc_ratio = sqrt(r / (2 * R));

    gmsh::model::occ::addCylinder(xc, yc, zc, 0, 0, h, R / 2, 1);

    gmsh::model::occ::synchronize();

    gmsh::model::mesh::field::add("Cylinder", 1);
    gmsh::model::mesh::field::setNumber(1, "VIn", lc * lc_ratio);
    gmsh::model::mesh::field::setNumber(1, "VOut", lc);
    gmsh::model::mesh::field::setNumber(1, "XCenter", xc);
    gmsh::model::mesh::field::setNumber(1, "YCenter", yc);
    gmsh::model::mesh::field::setNumber(1, "ZCenter", zc);
    gmsh::model::mesh::field::setNumber(1, "ZAxis", h);
    gmsh::model::mesh::field::setNumber(1, "Radius", r / 2);

    gmsh::model::mesh::field::setAsBackgroundMesh(1);

    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

    gmsh::option::setNumber("Mesh.MshFileVersion", 2.2); //для конвертации в .xml с помощью dolfin-а
    gmsh::option::setNumber("Mesh.Algorithm", 5);

    gmsh::model::mesh::generate(3);

    gmsh::write("cylinder.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();

    return 0;
}
