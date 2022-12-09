#define SHAPO_NOEXTERN 1
#include<shapo/hydro/stokes.hpp>

namespace Feel::shapo::hydro
{
    //using namespace Feel;

    template class stokesShapeOpt<Simplex<3>, Lagrange<2, Vectorial>, Lagrange<2, Vectorial>, Lagrange<1, Vectorial>>;
    //template class stokesShapeOpt<Simplex<3>, Lagrange<3, Vectorial>, Lagrange<3, Vectorial>, Lagrange<2, Vectorial>>;
} // namespace Feel::shapo::hydro