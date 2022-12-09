namespace Feel::shapo::hydro 
{
    template< typename ConvexType, typename PrimalBasisType, typename AdjointBasisType, typename DescentBasisType>
    class ShapeOpt
    {
    public:     

        // Mesh of the computational domain - types
        typedef ConvexType convex_type;
        typedef Mesh<convex_type> mesh_type;
        typedef std::shared_ptr<mesh_type> mesh_ptrtype;
        typedef typename mesh_type::trace_mesh_type trace_mesh_type;
        typedef typename mesh_type::trace_mesh_ptrtype trace_mesh_ptrtype;

        // Finite element spaces and fields - types
        // Primal problem
        typedef FunctionSpace<mesh_type, bases<PrimalBasisType> > primal_space_type;
        typedef std::shared_ptr<primal_space_type> primal_space_ptrtype;
        typedef typename primal_space_type::element_type element_primal_type;
        typedef std::shared_ptr<element_primal_type> element_primal_ptrtype;
        // Adjoint problem
        typedef FunctionSpace<mesh_type, bases<AdjointBasisType> > adjoint_space_type;
        typedef std::shared_ptr<adjoint_space_type> adjoint_space_ptrtype;
        typedef typename adjoint_space_type::element_type element_adjoint_type;
        typedef std::shared_ptr<element_adjoint_type> element_adjoint_ptrtype;

        mesh_ptrtype M_mesh,mesh_initial;
        primal_space_ptrtype primalSpace;
        adjoint_space_ptrtype adjointSpace;
        element_primal_ptrtype primalVariable;
        element_adjoint_ptrtype adjointVariable;

        // Parameters of the augmented Lagrangian algorithm
        double lnp1; // Lagrange multiplier
        double bnp1; // penalty parameter
        int itmax=10; // max number of iterations in for the shape optimization
        double tolerance = 1.; // L2 norm tolerance - stopping criterion
        double alpha; // multiplicative factor for penalty parameter in Augmented Lagrangian
        double volumeConstraint = 1; // volume at time 0 of the shape
        double volShape = 1.1; // current volume of the shape
        double volCompDomain = 1; // volume of the whole computational domain
        double btarget = 0; // maximal value of the penalty parameter
        double stepdescent=1; // value of the descent step
        double cost_function; // value of the Augmented Lagrangian
        double cost_function_save; // saved value of the cost function
        double initial_cost_function; // initial value of the cost function 
        double initial_constraint_value; // initial value of the constraints
        double initialShapeGrad; // initial value of shape gradient
        bool compute_initial_cost_function=false; // compute initial cost function for normalization of the cost
        bool compute_initial_volume=false; // compute initial value of constraints for normalization of the constraints
        bool computeInitialShapeGrad=false; // compute initial value of shape Grad

        double mu=1.; // fluid viscosity 


        double minOrMax; // 1 to minimize the cost and -1 to maximize the cost
        bool allExporters; // true to export primal, adjoint, descent results; else export only primal
        int niter_augment_b=1; // tells how often the penalty parameter should be updated

        // Linesearch parameters
        double c_ls=0.1; // 0<c<1
        double tau_ls=0.1; // 0<tau<1
        double alpha_ls = 0.1; // alpha>0

        // // Collect cost function, constraints and descent L2 norm
        // std::vector<double> history_cost_function={};
        // std::vector<double> history_constraints={};
        // std::vector<double> history_l2_norm_disp={};

        // Descent direction problem - types
        typedef FunctionSpace<mesh_type, bases<DescentBasisType> > descent_space_type;
        typedef std::shared_ptr<descent_space_type> descent_space_ptrtype;
        typedef typename descent_space_type::element_type element_descent_type;
        typedef std::shared_ptr<element_descent_type> element_descent_ptrtype;

        descent_space_ptrtype descentSpace;
        element_descent_ptrtype descentVariable;
        element_descent_ptrtype shapeGradient;
    
        // Exporters
        typedef Exporter<mesh_type,1> export_type;
        typedef std::shared_ptr<export_type> export_ptrtype;        
        typedef Exporter<trace_mesh_type,1> export_trace_type;
        typedef std::shared_ptr<export_trace_type> export_trace_ptrtype;

        export_ptrtype exporterPrimal, exporterAdjoint, exporterDescent,exporterOpti;

        // File to save cost function, constraint value and L2 norm
        std::ofstream M_fileresults;
        std::string M_fileresults_name;

        void init();
        virtual void solvePrimalProblem(int i){};
        virtual void solveAdjointProblem(int i){};
        virtual void solveDescentDirection(int i){};
        virtual void computeShapeGradient(int i){};
        virtual void findDescentStep(int i){};        
        virtual void updateMesh(){};
        virtual double computeCostFunction(int i,int j)=0;
        void updateLagrangianParameters(int i);
        void updateHistory(double l2norm);
        bool ArmijoCondition(double current_cost, double initial_cost, double mConstant);
        mesh_ptrtype mesh() {return this->M_mesh;}
        
        // Routine solving the shape optimization problem
        void solveShapeOptimization()
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

    };


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


} // end namespace Feel
