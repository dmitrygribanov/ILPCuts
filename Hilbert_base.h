#pragma once

#include "base_module.h"
#include "Eigen_solving_utilities.h"
#include "ntl_module.h"

namespace ILPCuts
{
    namespace Utility
    {
        template <typename int_type>
        SparseMatrix<int_type> compute_Hilbert_basis(const SparseMatrix<int_type> &A);
    } // namespace Utility

} // namespace ILPCuts
