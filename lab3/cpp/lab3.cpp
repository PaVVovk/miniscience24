#include <dolfin.h>
#include "heat.h"

using namespace dolfin;

// Source term (f description)
class Source : public Expression {
public:
    double t = 0.0;
private:
    double tau = 1.0;
    double ampl = 1e3;
    double sigma = 3.0;

    void eval(Array<double> &values, const Array<double> &x) const
    {
        double dx = x[0];
        double dy = x[1] - 11;
        double dz = x[2] - 20;
        double gauss = exp(-(dx * dx + dy * dy + dz * dz) / (2 * sigma * sigma)) / pow(2 * M_PI * sigma * sigma, 3/2);
        values[0] = ampl * gauss * exp(-t / tau);
    }
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
    bool inside(const Array<double> &x, bool on_boundary) const
    {
        return on_boundary;
    }
};

int main()
{
    // Create mesh and function space
    auto mesh = std::make_shared<Mesh>("cat.xml");
    auto V = std::make_shared<heat::FunctionSpace>(mesh);

    // Define boundary condition
    auto u_d = std::make_shared<Constant>(36.6);
    auto boundary = std::make_shared<DirichletBoundary>();
    DirichletBC bc(V, u_d, boundary);

    double T = 15.0;
    double steps = 450;
    double dt = T / steps;
    double cond = 10;

    auto u0 = std::make_shared<Function>(V);
    u0->interpolate(*u_d);

    double t = 0.0;

    // Define variational forms
    heat::BilinearForm a(V, V);
    heat::LinearForm L(V);
    auto f = std::make_shared<Source>();
    
    a.c = std::make_shared<Constant>(cond);
    a.dt = std::make_shared<Constant>(dt);

    L.dt = std::make_shared<Constant>(dt);
    L.f = f;
    L.u0 = u0;
    
    File file("heat_cpp/heat.pvd");

    for (int i = 0; i < steps; i++){
        t += dt;
        f->t = t;

        // Compute solution
        auto u1 = std::make_shared<Function>(V);
        solve(a == L, *u1, bc);

        // Save solution in VTK format
        
        file << *u1;

        assign(u0, u1);
    }
    return 0;
}
