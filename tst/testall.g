#
# LoopIntegrals
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#

LoadPackage( "FinSetsForCAP", ">= 2019.12.15" );

LoadPackage( "IO_ForHomalg" );
HOMALG_IO.show_banners := false;

options := rec(
    exitGAP := true,
    testOptions := rec(
        compareFunction := "uptowhitespace"
    ),
);

LoadPackage( "LoopIntegrals" );

TestDirectory( DirectoriesPackageLibrary( "LoopIntegrals", "tst" ), options );

FORCE_QUIT_GAP( 1 ); # if we ever get here, there was an error
