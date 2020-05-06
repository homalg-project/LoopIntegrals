#
# LoopIntegrals
#
# Reading the implementation part of the package.
#

ReadPackage( "LoopIntegrals", "gap/LoopIntegrals.gi");

if IsPackageMarkedForLoading( "JuliaInterface", ">= 0.3" ) then
    ReadPackage( "LoopIntegrals", "gap/Julia.gi" );
fi;
