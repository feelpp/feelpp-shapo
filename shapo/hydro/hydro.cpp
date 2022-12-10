#define SHAPO_NOEXTERN 1

#include <shapo/hydro/hydro.hpp>
#include <feel/feelfilters/loadmesh.hpp>

namespace Feel::shapo::hydro {

    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    ShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>::ShapeOpt()
    {
        this->init();
    }
    // Routine initialising the variables
    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    void ShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>::init()
    {
        LOG(INFO) << fmt::format("[ShapeOpt::init] starts") << std::endl;
        // Define the mesh for the fluid domain
        this->mesh_initial = loadMesh(_mesh=new mesh_type,_filename=soption("gmsh.filename"));
        M_mesh = createSubmesh(_mesh=this->mesh_initial,_range=markedelements(this->mesh_initial,"Fluid"));
        LOG(INFO) << fmt::format("[ShapeOpt::init] mesh loaded and extracted Fluid submesh") << std::endl;

        primalSpace = primal_space_type::New(_mesh=this->mesh());
        adjointSpace = adjoint_space_type::New(_mesh=this->mesh());
        descentSpace = descent_space_type::New(_mesh=this->mesh());

        primalVariable = primalSpace->elementPtr("primalVariable");
        adjointVariable = adjointSpace->elementPtr("adjointVariable");
        descentVariable = descentSpace->elementPtr("descentVariable");
        shapeGradient = descentSpace->elementPtr("shapeGradient");
        LOG(INFO) << fmt::format("[ShapeOpt::init] spaces and elements created") << std::endl;

        exporterPrimal = exporter(_mesh=this->mesh(),_name="exporterPrimal",_geo="change");
        exporterAdjoint = exporter(_mesh=this->mesh(),_name="exporterAdjoint",_geo="change");
        exporterDescent = exporter(_mesh=this->mesh(),_name="exporterDescent",_geo="change");
        exporterOpti= exporter(_mesh=this->mesh(),_name="exporterOpti",_geo="change");
        LOG(INFO) << fmt::format("[ShapeOpt::init] exporters created") << std::endl;
        
        // Read the values of these variables from the options
        lnp1 = doption("l0");
        bnp1 = doption("b0");
        itmax=ioption("itmax");
        tolerance=doption("tolerance");
        alpha=doption("alpha");
        stepdescent = doption("stepdescent");
        btarget=doption("btarget");
        mu=doption("mu");    
        allExporters = boption("all-exporters");
        if(boption("minimize-objective-function"))
            minOrMax=1;
        else 
            minOrMax=-1;
        niter_augment_b=ioption("niter_augment_b");

        // Linesearch parameters
        double c_ls=doption("linesearch.c"); // 0<c<1
        double tau_ls=doption("linesearch.tau"); // 0<tau<1
        double alpha_ls = doption("linesearch.alpha");// alpha>0

        // Compute the volume of the whole domain
        this->volCompDomain=integrate(_range=elements(this->mesh_initial),_expr=cst(1.)).evaluate()(0,0);
        // Compute the volume of the shape
        volumeConstraint =this->volCompDomain;
        volumeConstraint -= integrate(_range=markedelements(this->mesh(),"Fluid"),_expr=cst(1.)).evaluate()(0,0);
        volShape = volumeConstraint;
        LOG(INFO) << fmt::format("[ShapeOpt::init] volume of the shape computed: {} constraint: {}", this->volCompDomain, volumeConstraint) << std::endl;

        cost_function = 1e5;     
        
        std::cout << "Volume of the initial shape " << this->volShape <<std::endl;
        if(!this->compute_initial_volume)
        {
            this->compute_initial_volume=true;
            this->initial_constraint_value=this->volShape;
        }   

        // Initialise the csv file for results storage
        M_fileresults_name=soption("file-result-name");
        std::ofstream M_fileresults( M_fileresults_name.c_str(), std::ios::trunc );
        M_fileresults<<"Cost function"<<","<<"Volume of the shape"<<","<<"L2 norm of the descent step"<<"\n";
        M_fileresults.close();
        LOG(INFO) << fmt::format("[ShapeOpt::init] done") << std::endl;
    }


    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    void ShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>::solveShapeOptimization()
    {
        this->init();
        double L2normDescent = 10;
        int i=0;
        while(L2normDescent> this->tolerance && i< this->itmax)
        {
            this->solvePrimalProblem(i);
            this->solveAdjointProblem(i);
            this->computeShapeGradient(i);                                
            this->solveDescentDirection(i);                    
            this->findDescentStep(i);
            L2normDescent = integrate(_range = markedfaces(this->mesh(),"Shape"),_expr=inner(idv(this->descentVariable),idv(this->descentVariable))).evaluate()(0,0);
            L2normDescent = math::sqrt(L2normDescent);
            this->updateMesh();
            this->updateLagrangianParameters(i);
            std::cout << "L2 norm descent" << L2normDescent <<std::endl;
            i++;
            this->updateHistory(L2normDescent);
            std::cout << "i= " << i << std::endl;
        }
    }


    // Routine updating the Lagrange multiplier and penalty parameter for the augmented Lagrangian
    // algorithm
    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    void ShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>::updateLagrangianParameters(int iter)
    {
        // Compute the current volume of the shape
        this->volShape = this->volCompDomain;
        this->volShape -= integrate(_range=markedelements(this->mesh(),"Fluid"),_expr=cst(1.)).evaluate()(0,0);
        
        std::cout << "Volume constraint " << this->volShape - this->volumeConstraint << std::endl; 
        std::cout << "Relative volume constraint " << ( this->volShape - this->volumeConstraint )/this->initial_constraint_value << std::endl; 
        
        // Update the Lagrange multiplier (l = l + b*(|Body|-|Body|_0)/normalization_factor)
        this->lnp1 += this->bnp1*( this->volShape - this->volumeConstraint )/this->initial_constraint_value;

        // Update the penalty parameter if it's less than btarget and at the appropriate iteration steps
        if (this->bnp1<this->btarget && iter%this->niter_augment_b==0)
        {
            this->bnp1 *=alpha;
        }
        std::cout << "lnp1 " << this->lnp1 << "bnp1 " << this->bnp1 << std::endl;                
    }

    // Routine saving the cost function, constraints and L2 norm of the descent function
    // in a csv file
    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    void ShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>::updateHistory(double l2norm)
    {
        // this->history_cost_function.push_back(this->cost_function);
        // this->history_constraints.push_back(this->volShape);
        // this->history_l2_norm_disp.push_back(l2norm);

        M_fileresults.open( M_fileresults_name.c_str(), std::ios::app );
        M_fileresults<<this->cost_function_save<<","<<this->volShape<<","<<l2norm<<"\n";
        M_fileresults.close();
    }

    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    bool ShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>::ArmijoCondition(double current_cost, double initial_cost,double mConstant)
    {
        return (- current_cost + initial_cost< this->alpha_ls*this->c_ls*mConstant);
    }

    
    template class ShapeOpt<Simplex<3>, Lagrange<2, Vectorial>, Lagrange<2, Vectorial>, Lagrange<1, Vectorial>>;
}