#define SHAPO_NOEXTERN 1
#include<shapo/hydro/hydro.hpp>

namespace shapo::hydro
{
    template class stokesShapeOpt<Simplex<3>, Lagrange<2, Vectorial>, Lagrange<2, Vectorial>, Lagrange<1, Vectorial>>;
    //template class stokesShapeOpt<Simplex<3>, Lagrange<3, Vectorial>, Lagrange<3, Vectorial>, Lagrange<2, Vectorial>>;
} // namespace shapo::hydro