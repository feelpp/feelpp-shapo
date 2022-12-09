#define SHAPO_NOEXTERN 1
#include<shapo/hydro/hydro.hpp>

namespace shapo
{
    namespace hydro
    {
        template class stokesShapeOpt<Simplex<3>, Lagrange<2, Vectorial>, Lagrange<2, Vectorial>, Lagrange<1, Vectorial>>;
        //template class stokesShapeOpt<Simplex<3>, Lagrange<3, Vectorial>, Lagrange<3, Vectorial>, Lagrange<2, Vectorial>>;

    } // namespace hydro

} // namespace shapo