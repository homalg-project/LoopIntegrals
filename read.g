# SPDX-License-Identifier: GPL-2.0-or-later
# LoopIntegrals: Compute master integrals using commutative and noncommutative methods from computational algebraic geometry
#
# Reading the implementation part of the package.
#

ReadPackage( "LoopIntegrals", "gap/LoopIntegrals.gi");

if IsPackageMarkedForLoading( "JuliaInterface", ">= 0.3" ) then
    ReadPackage( "LoopIntegrals", "gap/Julia.gi" );
fi;
