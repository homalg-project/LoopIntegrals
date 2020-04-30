#! @Chunk 2LoopCuspAnomalousDimensions

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "k1, k2", "v1, v2", 2 );
#! <A loop diagram with loop momenta [ k1, k2 ] & external momenta [ v1, v2 ]>
cos:= v1*v2;;
SetAbbreviation( cos, "cos" );
rel := [ v1^2 - 1, v2^2 - 1 ];;
SetRelationsOfMomenta( LD, rel );
SetIndependentLorentzInvariants( LD, [ k1^2, k1*k2, k1*v1, k1*v2, k2^2, k2*v1, k2*v2, cos ] );
SetPropagators( LD, -[ 2*k2*v1 - 1, 2*k2*v2 - 1, (k1 - k2)^2, 2*k1*v1 - 1, 2*k1*v2 - 1, k1^2 ] );
SetNumerators( LD, -[ k2^2 ] );
SetExtraLorentzInvariants( LD, [ cos ] );
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <An unevaluated non-zero 8 x 6 matrix over an external ring>,
#!   <An unevaluated non-zero 6 x 6 matrix over an external ring> ]
S := SyzygiesOfRows( E12[1], E12[2] );
#! <A non-zero 68 x 8 matrix over an external ring>
EntriesOfHomalgMatrix( S[1] );
#! [ 2*D6*cos+2*D6, 0, -D6, -D6, 0, 2*D6*cos+2*D6, -D6, -D6 ]
EntriesOfHomalgMatrix( S[2] );
#! [ 0, 0, 0, 0, 0, 2*D3*cos+2*D3, -D3, -D3 ]
#! @EndExample
