gap> LoadPackage( "LoopIntegrals", false );
true
gap> package_loading_info_level := InfoLevel( InfoPackageLoading );;
gap> SetInfoLevel( InfoPackageLoading, PACKAGE_INFO );;
gap> LoadPackage( "LoopIntegrals" );
true
gap> SetInfoLevel( InfoPackageLoading, package_loading_info_level );;
