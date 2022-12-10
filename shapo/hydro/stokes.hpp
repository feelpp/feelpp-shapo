
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
        using super = ShapeOpt<ConvexType,PrimalBasisType,AdjointBasisType,DescentBasisType>;
        public:
            
            stokesShapeOpt()
                : super()
            {
            }
            void solvePrimalProblem(int i); 
            void solveAdjointProblem(int i);
            void computeShapeGradient(int i);
            void solveDescentDirection(int i);        
            void updateMesh();
            void findDescentStep(int i);
            double computeCostFunction(int i, int j);
    };

    
#if !defined(SHAPO_NOEXTERN)
extern template class stokesShapeOpt<Simplex<3>, Lagrange<2, Vectorial>, Lagrange<2, Vectorial>, Lagrange<1, Vectorial>>;
#endif

} // end namespace Feel
