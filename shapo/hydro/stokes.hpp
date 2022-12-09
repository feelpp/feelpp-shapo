
#include<feel/feelfilters/loadmesh.hpp>
#include<feel/feelvf/on.hpp>
#include <feel/feelmesh/meshmover.hpp>
#include<feel/feeldiscr/createsubmesh.hpp>
#include <feel/feeldiscr/quality.hpp>
#include <feel/feelmesh/remesher.hpp>

#include <shapo/hydro/hydro.hpp>

// NB: All normals are inversed N() -> -N() because the normal in the paper is pointing outward with respect to the shape; 
// here N() points inside the the empty shape
namespace Feel::shapo::hydro
{
// ------------- STOKES SHAPEOPT CLASS -----------
    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    class stokesShapeOpt : public ShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>
    {
        public:
            typedef ShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType> super_type;
            void solvePrimalProblem(int i); 
            void solveAdjointProblem(int i);
            void computeShapeGradient(int i);
            void solveDescentDirection(int i);        
            void updateMesh();
            void findDescentStep(int i);
            double computeCostFunction(int i, int j);
    };

    // Routine initialising the variables
    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    void ShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>::init()
    {
        // Define the mesh for the fluid domain
        this->mesh_initial = loadMesh(_mesh=new mesh_type,_filename=soption("gmsh.filename"));
        M_mesh = createSubmesh(_mesh=this->mesh_initial,_range=markedelements(this->mesh_initial,"Fluid"));

        primalSpace = primal_space_type::New(_mesh=this->mesh());
        adjointSpace = adjoint_space_type::New(_mesh=this->mesh());
        descentSpace = descent_space_type::New(_mesh=this->mesh());

        primalVariable = primalSpace->elementPtr("primalVariable");
        adjointVariable = adjointSpace->elementPtr("adjointVariable");
        descentVariable = descentSpace->elementPtr("descentVariable");
        shapeGradient = descentSpace->elementPtr("shapeGradient");

        exporterPrimal = exporter(_mesh=this->mesh(),_name="exporterPrimal",_geo="change");
        exporterAdjoint = exporter(_mesh=this->mesh(),_name="exporterAdjoint",_geo="change");
        exporterDescent = exporter(_mesh=this->mesh(),_name="exporterDescent",_geo="change");
        exporterOpti= exporter(_mesh=this->mesh(),_name="exporterOpti",_geo="change");

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
    }

    // Routine for the solution of the primal Stokes problem
    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    void stokesShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>::solvePrimalProblem(int i)
    {
        auto Ph = Pch<PrimalBasisType::nOrder-1>(this->mesh());
        auto p = Ph->elementPtr();
        auto q = Ph->element();
        auto v = this->primalSpace->element();
        
        // Define the matrices and vectors
        BlocksBaseGraphCSR myblockGraph(2,2);
        BlocksBaseVector<double> myblockVec(2);
        BlocksBaseVector<double> myblockVecSol(2);

        myblockGraph(0,0) = stencil(_test=this->primalSpace,_trial=this->primalSpace, _diag_is_nonzero=false, _close=false)->graph();
        myblockGraph(0,1) = stencil(_test=this->primalSpace,_trial=Ph, _diag_is_nonzero=false, _close=false)->graph();
        myblockGraph(1,0) = stencil(_test=Ph,_trial=this->primalSpace, _diag_is_nonzero=false, _close=false)->graph();        
        auto A = backend()->newBlockMatrix(_block=myblockGraph);
        myblockVec(0,0) = backend()->newVector( this->primalSpace );
        myblockVec(1,0) = backend()->newVector( Ph );
        auto F = backend()->newBlockVector(_block=myblockVec, _copy_values=false);

        myblockVecSol(0,0) = this->primalVariable;
        myblockVecSol(1,0) = p;
        auto U = backend()->newBlockVector(_block=myblockVecSol, _copy_values=false);
        
        // grad grad
        form2(_trial=this->primalSpace,_test=this->primalSpace,_matrix=A,
         _rowstart=0, _colstart=0) = integrate(_range=markedelements(this->mesh(),"Fluid"),
                                        _expr = cst(-2*this->mu)*inner( sym( gradt(this->primalVariable) ) ,sym( grad(v))));
        
        // div terms                                
        form2(_trial=this->primalSpace,_test=Ph,_matrix=A,
        _rowstart=1, _colstart=0) += integrate(_range=markedelements(this->mesh(),"Fluid"),
                                            _expr = inner(divt(this->primalVariable),id(q)));
        form2(_trial=Ph,_test=this->primalSpace,_matrix=A,
        _rowstart=0, _colstart=1) += integrate(_range=markedelements(this->mesh(),"Fluid"),
                                            _expr = inner(div(v),idt(p)) );
        
        // boundary conditions
        form2( _trial=this->primalSpace, _test=this->primalSpace ,_matrix=A )
                +=on(_range=markedfaces(this->mesh(),"External"), _rhs=F, _element=*this->primalVariable,
                _expr=oneX()*cst(0.0));
        form2( _trial=this->primalSpace, _test=this->primalSpace ,_matrix=A )
                +=on(_range=markedfaces(this->mesh(),"Shape"), _rhs=F, _element=*this->primalVariable,
                _expr=expr<3,1>(soption("Usolid")));                        

        // Solve the PDE
        backend(_name="SP",_rebuild=true)->solve( _matrix=A, _rhs=F, _solution=U );
        myblockVecSol.localize(U);

        // Post-processing
        // Compute the normal stresses 
        auto normal_stresses = this->primalSpace->element();
        auto boundary_speed = this->primalSpace->elementPtr();
        normal_stresses.on(_range=markedfaces(this->mesh(),"Shape"), _expr =(-idv(p)*Id<FEELPP_DIM>() + cst(2*this->mu)*sym(gradv(this->primalVariable)))*(-N()) );
        boundary_speed->on(_range=markedfaces(this->mesh(),"Shape"), _expr = expr<3,1>(soption("Vsolid")));
        BlocksBaseVector<double> myblockStresses(2);

        myblockStresses(0,0) = boundary_speed;
        myblockStresses(1,0) = backend()->newVector( Ph );
        auto stresses_larger = backend()->newBlockVector(_block=myblockStresses, _copy_values=true);
        myblockStresses.localize(stresses_larger);

        // Compute the cost function
        this->cost_function = integrate(_range=markedfaces(this->mesh(),"Shape"),_expr=inner(-idv(normal_stresses),idv(boundary_speed))).evaluate()(0,0);        
    
        std::cout << "cost function primal solution: " << this->cost_function << std::endl;
        if(!this->compute_initial_cost_function)
        {
            this->compute_initial_cost_function=true;
            this->initial_cost_function=this->cost_function;
        }

        // Normalize the cost function
        this->cost_function = this->cost_function/this->initial_cost_function;

        this->cost_function_save = this->cost_function;

        // Exporters
        if(i<100)
        {
            auto normalGradientprod = Ph->element();
            normalGradientprod.on(_range=elements(this->mesh()),_expr=inner(sym(gradv(this->primalVariable)),sym(gradv(this->primalVariable)))/this->initial_cost_function+cst(this->lnp1/this->initial_constraint_value)+cst(this->bnp1/(this->initial_constraint_value*this->initial_constraint_value))*(cst(this->volShape)-cst(this->volumeConstraint)));            
            auto normalGradientprodOnly= Ph->element();
            normalGradientprodOnly.on(_range=elements(this->mesh()),_expr=inner(sym(gradv(this->primalVariable)),sym(gradv(this->primalVariable)))/this->initial_cost_function);            
            auto constraintOnly= Ph->element();
            constraintOnly.on(_range=elements(this->mesh()),_expr=inner(cst(this->lnp1/this->initial_constraint_value)+cst(this->bnp1/(this->initial_constraint_value*this->initial_constraint_value))*(cst(this->volShape)-cst(this->volumeConstraint))));            
            this->exporterPrimal->step(i)->setMesh(this->mesh());
            this->exporterPrimal->step(i)->add("velocityPrimal",this->primalVariable);
            this->exporterPrimal->step(i)->add("pressurePrimal",p);
            this->exporterPrimal->step(i)->add("normalGradientprod",normalGradientprod);
            this->exporterPrimal->step(i)->add("normalGradientprodOnly",normalGradientprodOnly);
            this->exporterPrimal->step(i)->add("constraintOnly",constraintOnly);
            this->exporterPrimal->save();
        }
        if(i>100 && this->allExporters)
        {
            this->exporterOpti->step(i)->setMesh(this->mesh());
            this->exporterOpti->step(i)->add("velocityPrimal",this->primalVariable);
            this->exporterOpti->step(i)->add("pressurePrimal",p);
            this->exporterOpti->step(i)->add("stressN",idv(normal_stresses),markedfaces(this->mesh(),"Shape"),"nodal");
            this->exporterOpti->step(i)->add("boundarySpeed",idv(boundary_speed),markedfaces(this->mesh(),"Shape"),"nodal");
            this->exporterOpti->save();
        }

    }

    // Routine for the solution of the adjoint Stokes problem
    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    void stokesShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>::solveAdjointProblem(int i)
    {
        auto Ph = Pch<PrimalBasisType::nOrder-1>(this->mesh());
        auto p = Ph->elementPtr();
        auto q = Ph->element();
        auto v = this->adjointSpace->element();        

        // Create matrices and vectors
        BlocksBaseGraphCSR myblockGraph(2,2);
        BlocksBaseVector<double> myblockVec(2);
        BlocksBaseVector<double> myblockVecSol(2);
                
        myblockGraph(0,0) = stencil(_test=this->adjointSpace,_trial=this->adjointSpace, _diag_is_nonzero=false, _close=false)->graph();
        myblockGraph(0,1) = stencil(_test=this->adjointSpace,_trial=Ph, _diag_is_nonzero=false, _close=false)->graph();
        myblockGraph(1,0) = stencil(_test=Ph,_trial=this->adjointSpace, _diag_is_nonzero=false, _close=false)->graph();        
        auto A = backend()->newBlockMatrix(_block=myblockGraph);

        myblockVec(0,0) = backend()->newVector( this->adjointSpace );
        myblockVec(1,0) = backend()->newVector( Ph );
        auto F = backend()->newBlockVector(_block=myblockVec, _copy_values=false);

        myblockVecSol(0,0) = this->adjointVariable;
        myblockVecSol(1,0) = p;
        auto U = backend()->newBlockVector(_block=myblockVecSol, _copy_values=false);
        
        // grad grad
        form2(_trial=this->adjointSpace,_test=this->adjointSpace,_matrix=A,
         _rowstart=0, _colstart=0) = integrate(_range=markedelements(this->mesh(),"Fluid"),
                                        _expr = cst(-2*this->mu)*inner( sym( gradt(this->adjointVariable) ) ,sym( grad(v))));
        // div terms                                
        form2(_trial=this->adjointSpace,_test=Ph,_matrix=A,
        _rowstart=1, _colstart=0) += integrate(_range=markedelements(this->mesh(),"Fluid"),
                                            _expr = inner(divt(this->adjointVariable),id(q)));
        form2(_trial=Ph,_test=this->adjointSpace,_matrix=A,
        _rowstart=0, _colstart=1) += integrate(_range=markedelements(this->mesh(),"Fluid"),
                                            _expr = inner(div(v),idt(p)) );
        
        // Boundary conditions
        form2( _trial=this->adjointSpace, _test=this->adjointSpace ,_matrix=A )
                +=on(_range=markedfaces(this->mesh(),"External"), _rhs=F, _element=*this->adjointVariable,
                _expr=oneX()*cst(0.0));
        form2( _trial=this->adjointSpace, _test=this->adjointSpace ,_matrix=A )
                +=on(_range=markedfaces(this->mesh(),"Shape"), _rhs=F, _element=*this->adjointVariable,
                _expr=expr<3,1>(soption("Vsolid")));
        
        // Solve the PDE
        backend(_name="SA",_rebuild=true)->solve( _matrix=A, _rhs=F, _solution=U );
        myblockVecSol.localize(U);

        // Export the results
        if(this->allExporters)
        {
            this->exporterAdjoint->step(i)->setMesh(this->mesh());
            this->exporterAdjoint->step(i)->add("velocityAdjoint",this->adjointVariable);
            this->exporterAdjoint->step(i)->add("pressureAdjoint",p);
            this->exporterAdjoint->save();
        }
    }

    // Routine for the computation of the descent direction via harmonic extension
    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    void stokesShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>::solveDescentDirection(int i)
    {
        auto B = form2(_trial=this->descentSpace,_test=this->descentSpace);
        auto g = form1(_test=this->descentSpace);
        auto v = this->descentSpace->elementPtr();        

        // grad grad
        B = integrate(_range=markedelements(this->mesh(),"Fluid"),_expr=inner(gradt(this->descentVariable),grad(v)));
        
        // right-hand side: normalized shape derivative and constraints 
        g = integrate(_range=markedfaces(this->mesh(),"Shape"),_expr=cst(-1/this->initial_cost_function)*inner(idv(this->shapeGradient),id(v)));
        g +=integrate(_range=markedfaces(this->mesh(),"Shape"),_expr=inner((-cst(this->lnp1/this->initial_constraint_value)-cst(this->bnp1/(this->initial_constraint_value*this->initial_constraint_value))*(cst(this->volShape)-cst(this->volumeConstraint)))*(-N()),id(v)));
        
        // Boundary conditions
        B+=on(_range=markedfaces(this->mesh(),"External"), _rhs=g, _element=*this->descentVariable,
                _expr=oneX()*cst(0.0));
        
        // Solve the PDE
        B.solve(_solution=this->descentVariable,_rhs=g,_name="SD",_rebuild=true);
        
        // Exporter
        if(this->allExporters)
        {
            this->exporterDescent->step(i)->setMesh(this->mesh());
            this->exporterDescent->step(i)->add("velocityDescent",this->descentVariable);
            this->exporterDescent->step(i)->add("shapeGradient",idv(this->shapeGradient),markedfaces(this->mesh(),"Shape"),"nodal");
            this->exporterDescent->step(i)->add("constraint",(-cst(this->lnp1/this->initial_constraint_value)-cst(this->bnp1/(this->initial_constraint_value*this->initial_constraint_value))*(cst(this->volShape)-cst(this->volumeConstraint)))*(-N()),markedfaces(this->mesh(),"Shape"),"nodal");
            this->exporterDescent->save();
        }
    }

    // Routine computing the shape gradient
    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    void stokesShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>::computeShapeGradient(int i)
    {       
        // L2 projection of the shape gradient onto P2 piecewise continuous elements
        auto mass_mat = form2(_test=this->primalSpace,_trial=this->primalSpace);
        auto right_hs = form1(_test=this->primalSpace);
        auto p2_shapeGrad = this->primalSpace->element();
        mass_mat = integrate(_range=elements(this->mesh()),_expr=inner(idt(p2_shapeGrad),id(p2_shapeGrad)));

        this->shapeGradient->on(_range=markedfaces(this->mesh(),"Shape"),_expr=cst(2*this->mu)*(inner(sym(gradv(this->primalVariable)),sym(gradv(this->adjointVariable))))*(-N()));
        right_hs = integrate(_range=elements(this->mesh()),_expr=inner(idv(this->shapeGradient),id(p2_shapeGrad)));

        mass_mat.solve(_solution=p2_shapeGrad,_rhs=right_hs,_name="L2proj",_rebuild=true);

        // Choose minimization or maximization
        this->shapeGradient->on(_range=markedfaces(this->mesh(),"Shape"),_expr=cst(this->minOrMax)*idv(p2_shapeGrad));
        
        // if(!this->computeInitialShapeGrad)
        // {
        //     this->computeInitialShapeGrad=true;
        //     this->initialShapeGrad =integrate(_range=markedfaces(this->mesh(),"Shape"),_expr=inner(idv(this->shapeGradient),idv(this->shapeGradient))).evaluate()(0,0);
        //     this->initialShapeGrad = math::sqrt(this->initialShapeGrad);
        // }
        // Exporter
        if(this->allExporters)
        {
            this->exporterDescent->step(i)->setMesh(this->mesh());
            this->exporterDescent->step(i)->add("velocityDescent",this->descentVariable);
            this->exporterDescent->step(i)->add("shapeGradient",idv(this->shapeGradient),markedfaces(this->mesh(),"Shape"),"nodal");
            this->exporterDescent->save();
        }
    }


    // Routine computing the Augmented Lagrangian
    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    double stokesShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>::computeCostFunction(int i, int j)
    {
        this->solvePrimalProblem((j+1)*1000+i); // Update the cost function
        auto cost_function = this->cost_function;
        std::cout << "cost function" << this->cost_function <<std::endl;
        this->volShape = this->volCompDomain;
        this->volShape -= integrate(_range=markedelements(this->mesh(),"Fluid"),_expr=cst(1.)).evaluate()(0,0);
        std::cout << "Volume of the shape " << this->volShape <<std::endl; 

        // add the Lagrange multiplier and penalty parameter
        cost_function += this->lnp1*(this->volShape-this->volumeConstraint)/this->initial_constraint_value;
        std::cout << "lnp1 part of the cost function " << this->lnp1*(this->volShape-this->volumeConstraint)/this->initial_constraint_value <<std::endl;
        cost_function +=this->bnp1*0.5*(this->volShape-this->volumeConstraint)*(this->volShape-this->volumeConstraint)/(this->initial_constraint_value*this->initial_constraint_value);
        std::cout << "bnp1 part of the cost function " << this->bnp1*0.5*(this->volShape-this->volumeConstraint)*(this->volShape-this->volumeConstraint)/(this->initial_constraint_value*this->initial_constraint_value) <<std::endl;

        return cost_function;
    }

    // Routine choosing the descent step along the descent direction
    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    void stokesShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>::findDescentStep(int j)
    {
        int i =1;
        auto initial_cost = this->computeCostFunction(0,j);
        std::cout << "initial cost " << initial_cost << std::endl;
        auto current_cost = initial_cost +1;
        
        // Compute the descent step according to Backtracking Linesearch
        double mConstant = integrate(_range=markedfaces(this->mesh(),"Shape"),_expr=inner(idv(this->shapeGradient),idv(this->descentVariable))).evaluate()(0,0);
        mConstant += integrate(_range=markedfaces(this->mesh(),"Shape"),_expr=inner((cst(this->lnp1/this->initial_constraint_value)+cst(this->bnp1/(this->initial_constraint_value*this->initial_constraint_value))*(cst(this->volShape)-cst(this->volumeConstraint)))*(-N()),idv(this->descentVariable))).evaluate()(0,0);

        this->stepdescent = integrate(_range = markedfaces(this->mesh(),"Shape"),_expr=inner(idv(this->descentVariable),idv(this->descentVariable))).evaluate()(0,0);
        this->stepdescent = this->alpha_ls/math::sqrt(this->stepdescent);
        
        std::cout << " c times m " << this->c_ls*mConstant << std::endl;

        double current_volume = this->volCompDomain;
        
        while(this->ArmijoCondition(current_cost, initial_cost, mConstant) && this->stepdescent>1e-6)
        {                
            auto disp = this->descentSpace->elementPtr(cst(this->stepdescent)*idv(this->descentVariable));
            meshMove(this->M_mesh,*disp );
            std::cout << "Currentstepsize" << this->stepdescent <<std::endl;            
            current_volume = this->volCompDomain;
            current_volume -= integrate(_range=markedelements(this->mesh(),"Fluid"),_expr=cst(1.)).evaluate()(0,0);
            current_cost = this->computeCostFunction(i,j);            
            double mConstant = integrate(_range=markedfaces(this->mesh(),"Shape"),_expr=inner(idv(this->shapeGradient),idv(this->descentVariable))).evaluate()(0,0);
            mConstant += integrate(_range=markedfaces(this->mesh(),"Shape"),_expr=inner((cst(this->lnp1/this->initial_constraint_value)+cst(this->bnp1/(this->initial_constraint_value*this->initial_constraint_value))*(cst(this->volShape)-cst(this->volumeConstraint)))*(-N()),idv(this->descentVariable))).evaluate()(0,0);
            std::cout <<  this->stepdescent << std::endl;
            disp->scale(-1.);
            meshMove(this->M_mesh,*disp );            
            disp->scale(-1.);    
            std::cout << "current cost: " << current_cost <<std::endl;   
            i++;
            std::cout << "why is it so small???" << math::abs(( current_volume - this->volumeConstraint )/this->initial_constraint_value) << std::endl;
            //if(current_cost>= initial_cost || math::abs(( current_volume - this->volumeConstraint )/this->initial_constraint_value)>5e-2)
            if (this->ArmijoCondition(current_cost, initial_cost, mConstant) )
                this->stepdescent *= this->tau_ls;                     
        }        
        std::cout << "loop left with current cost "<< current_cost<<" and initial cost "<< initial_cost <<std::endl;
        std::cout << "loop left with step size "<< this->stepdescent<<std::endl;
    }

    // Routine updating the mesh via mesh adaptation and remeshing, if necessary
    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    void stokesShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>::updateMesh()
    {
        // Mesh motion
        auto disp = this->descentSpace->elementPtr(cst(this->stepdescent)*idv(this->descentVariable));
        meshMove(this->M_mesh,*disp );        
        auto etaqmin = etaQ( this->mesh() )->min();
        std::cout << "etaqmin " << etaqmin << std::endl;

        // Remesh if quality is low, and reconstruct all the spaces
        if (etaqmin < 0.4)
        {                        
            auto r = remesher(this->M_mesh,std::vector<int>{},std::vector<int>{});
            std::cout << this->mesh()->markerNames() << std::endl;
            auto new_mesh_remeshed = r.execute();
            std::cout << new_mesh_remeshed->markerNames() << std::endl;
            new_mesh_remeshed->updateForUse();
            new_mesh_remeshed->setMarkerNames(this->mesh()->markerNames());
            this->M_mesh = new_mesh_remeshed;
            this->primalSpace = super_type::primal_space_type::New(_mesh=this->mesh());
            this->adjointSpace = super_type::adjoint_space_type::New(_mesh=this->mesh());
            this->descentSpace = super_type::descent_space_type::New(_mesh=this->mesh());
            this->primalVariable = this->primalSpace->elementPtr("primalVariable");
            this->adjointVariable = this->adjointSpace->elementPtr("adjointVariable");
            this->descentVariable = this->descentSpace->elementPtr("descentVariable");
            this->shapeGradient = this->descentSpace->elementPtr("shapeGradient");
        }
            
    }
#if !defined(SHAPO_NOEXTERN)
extern template class stokesShapeOpt<Simplex<3>, Lagrange<2, Vectorial>, Lagrange<2, Vectorial>, Lagrange<1, Vectorial>>;
extern template class ShapeOpt<Simplex<3>, Lagrange<2, Vectorial>, Lagrange<2, Vectorial>, Lagrange<1, Vectorial>>;
#endif

} // end namespace Feel
