#include <set>
#include <cmath>
#include <gmsh.h>
#include <vector>


int main(int argc, char **argv)
{
    gmsh::initialize();

    gmsh::model::add("rect");

    const double xc = 0;
    const double yc = 0;
    const double zc = 0;

    const double R = 7e-3;
    const double r = 5e-4;
    const double h = 0.4;

    const double lc = 2.5 * 1e-4;
    const double lc_ratio = sqrt(r / (2 * R));

    gmsh::model::occ::addRectangle(xc, yc - R / 2, zc, h, R / 2, 1);
    gmsh::model::occ::addRectangle(xc, yc, zc, h, R / 2, 2);

    gmsh::model::occ::synchronize();

    gmsh::model::mesh::field::add("Box", 1);
    gmsh::model::mesh::field::setNumber(1, "VIn", lc * lc_ratio);
    gmsh::model::mesh::field::setNumber(1, "VOut", lc);
    gmsh::model::mesh::field::setNumber(1, "XMin", xc);
    gmsh::model::mesh::field::setNumber(1, "XMax", xc + h);
    gmsh::model::mesh::field::setNumber(1, "YMin", yc - r / 2);
    gmsh::model::mesh::field::setNumber(1, "YMax", yc + r / 2);

    gmsh::model::mesh::field::setAsBackgroundMesh(1);

    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

    gmsh::option::setNumber("Mesh.MshFileVersion", 2.2); //для конвертации в .xml с помощью dolfin-а
    gmsh::option::setNumber("Mesh.Algorithm", 5);

    gmsh::model::mesh::generate(2);

    gmsh::write("rect.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();

    return 0;
}
