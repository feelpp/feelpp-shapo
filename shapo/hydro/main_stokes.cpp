
#include <shapo/hydro/stokes.hpp>

int main(int argc, char** argv)
{
    using namespace Feel;
    using namespace Feel::shapo::hydro;

    Feel::po::options_description shape_optimization_options;
    shape_optimization_options.add_options()
        ("l0", Feel::po::value<double>()->default_value( 1.0 ),"Initial value of the Lagrange multiplier")
        ("b0", Feel::po::value<double>()->default_value( 1.0 ),"Penalization constant for the volume constraint")
        ("alpha", Feel::po::value<double>()->default_value( 1.5 ),"Multiplicative constant increasing the penalization constant")
        ("itmax", Feel::po::value<int>()->default_value( 10 ),"Max number of shape optimization iterations")
        ("tolerance", Feel::po::value<double>()->default_value( 1.0 ),"Tolerance for the descent direction")
        ("btarget", Feel::po::value<double>()->default_value( 10 ),"Maximal value of the augmented Lagrangian b parameter")
        ("mu", Feel::po::value<double>()->default_value( 1.0 ),"Viscosity in Stokes equations")
        ("stepdescent", Feel::po::value<double>()->default_value( 1.0 ),"Descent step for shape deformation")       
        ("Usolid", Feel::po::value<std::string>()->default_value( "{1,0,0}" ),"Boundary condition primal problem on Shape")
        ("Vsolid", Feel::po::value<std::string>()->default_value( "{1,0,0}" ),"Boundary condition dual problem on Shape")
        ("file-result-name",Feel::po::value<std::string>()->default_value( "results_shape_opti.csv" ),"File where saving cost, constraint, volume")
        ("minimize-objective-function",Feel::po::value<bool>()->default_value(true),"Minimize (true) or maximize (false) the cost function by changing the sign")
        ("all-exporters",Feel::po::value<bool>()->default_value(false),"Save only the primal exporter")
        ("linesearch.c",Feel::po::value<double>()->default_value(0.1),"c parameter for linesearch algo for the descent step")
        ("linesearch.tau",Feel::po::value<double>()->default_value(0.1),"tau parameter for linesearch algo for the descent step")
        ("linesearch.alpha",Feel::po::value<double>()->default_value(0.1),"alpha parameter for linesearch algo for the descent step")
        ("niter_augment_b",Feel::po::value<int>()->default_value( 5 ),"Max number of shape optimization iterations")
    ;
    shape_optimization_options.add( Feel::feel_options() );
    shape_optimization_options.add( Feel::backend_options("SP"));
    shape_optimization_options.add( Feel::backend_options("SA"));
    shape_optimization_options.add( Feel::backend_options("SD"));
    shape_optimization_options.add( Feel::backend_options("L2proj"));
    Environment env( _argc=argc, _argv=argv,
                     _desc=shape_optimization_options,
                     _worldcomm= Environment::worldComm());

    stokesShapeOpt<Simplex<3>,
    Lagrange<2, Vectorial>,
    Lagrange<2, Vectorial>,
    Lagrange<1, Vectorial>> shape_opt;

#if 1
    shape_opt.solveShapeOptimization();

#else
    shape_opt.init();
    shape_opt.solvePrimalProblem(0);
    shape_opt.solveAdjointProblem(0);
    auto diff_solution = integrate(_range=markedelements(shape_opt.mesh(),"Fluid"),_expr=idv(shape_opt.primalVariable)-idv(shape_opt.adjointVariable)).evaluate();
    auto diff_solution_2 = integrate(_range=markedelements(shape_opt.mesh(),"Fluid"),_expr=idv(shape_opt.primalVariable)).evaluate();
    std::cout << "Integral primal solution " << diff_solution_2 << std::endl;
    std::cout << "Difference between solutions " << diff_solution << std::endl;

    auto u = expr<3,1>("{4*sin(pi*x)*cos(2*pi*y),10*sin(3*pi*x)*cos(4*pi*y),0}:x:y");
    auto u_imposed = shape_opt.primalSpace->element();
    u_imposed.on(_range=elements(shape_opt.mesh()),_expr=u);
    auto shape_space = Pch<1>(shape_opt.mesh());
    auto shape_elem = shape_space->element();
    shape_elem.on(_range=elements(shape_opt.mesh()),_expr=inner(sym(gradv(u_imposed)),sym(gradv(u_imposed))));
    shape_opt.shapeGradient->on(_range=markedfaces(shape_opt.mesh(),"Shape"),_expr=inner(sym(gradv(u_imposed)),sym(gradv(u_imposed)))*N());

    auto ex = exporter(_mesh=shape_opt.mesh(),_name="exporter_test",_geo="change");
    ex->step(0)->setMesh(shape_opt.mesh());    
    ex->step(0)->add("u",u_imposed);
    ex->step(0)->add("shape_elem",idv(shape_elem),markedfaces(shape_opt.mesh(),"Shape"),"nodal");
    ex->step(0)->add("shapeGradient",idv(shape_opt.shapeGradient),markedfaces(shape_opt.mesh(),"Shape"),"nodal");
    ex->step(0)->add("expr_shape_grad",expr("pow(M_PI, 2)*(450.0*pow(0.26666666666666666*sin(M_PI*x)*sin(2*M_PI*y) - cos(3*M_PI*x)*cos(4*M_PI*y), 2) + 1600.0*pow(sin(3*M_PI*x), 2)*pow(sin(4*M_PI*y), 2) + 16.0*pow(cos(M_PI*x), 2)*pow(cos(2*M_PI*y), 2)):x:y"),markedfaces(shape_opt.mesh(),"Shape"),"nodal");
    ex->save();


#endif
    
    return 0;
}