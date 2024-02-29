#include <set>
#include <cmath>
#include <gmsh.h>
#include <vector>


int main(int argc, char **argv)
{
    gmsh::initialize();

    gmsh::model::add("cat");

    try
    {
        gmsh::merge("cat.stl");
    }
    catch(...)
    {
        gmsh::logger::write("Could not load STL mesh: bye!");
        gmsh::finalize();
        return 0;
    }

    // Angle between two triangles above which an edge is considered as sharp:
    double angle = 40;

    // For complex geometries, patches can be too complex, too elongated or too
    // large to be parametrized; setting the following option will force the
    // creation of patches that are amenable to reparametrization:
    bool forceParametrizablePatches = false;

    // For open surfaces include the boundary edges in the classification process:
    bool includeBoundary = true;

    // Force curves to be split on given angle:
    double curveAngle = 180;

    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
                                        forceParametrizablePatches,
                                        curveAngle * M_PI / 180.);

    // Create a geometry for all the discrete curves and surfaces in the mesh, by
    // computing a parametrization for each one
    gmsh::model::mesh::createGeometry();

    // Create a volume from all the surfaces
    std::vector<std::pair<int, int>> s;
    gmsh::model::getEntities(s, 2);
    std::vector<int> sl;
    for(auto surf : s) sl.push_back(surf.second);
    int l = gmsh::model::geo::addSurfaceLoop(sl);
    gmsh::model::geo::addVolume({l});

    gmsh::model::geo::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeFactor", 0.1);

    gmsh::model::mesh::generate(3);

    gmsh::write("cat.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();

    return 0;
}
