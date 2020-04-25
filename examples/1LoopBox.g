#! @Chunk 1LoopBox

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1", "k1..2,k4", 2 );
#! <A loop diagram with loop momenta [ l1 ] & external momenta [ k1, k2, k4 ]>
s12 := 2*k1*k2;;
SetAbbreviation( s12, "s12" );
s14 := 2*k1*k4;;
SetAbbreviation( s14, "s14" );
SetIndependentLorentzInvariants( LD,
        [ l1^2, l1*k1, l1*k2, l1*k4, s12, s14 ] );
rel1 := List( ExternalMomenta( LD ), k -> k^2 );;
rel2 := [ (k1+k2+k4)^2 ];;
SetRelationsOfMomenta( LD, Concatenation( rel1, rel2 ) );
SetPropagators( LD, [ l1^2, (l1-k1)^2, (l1-k1-k2)^2, (l1+k4)^2 ] );
SetNumerators( LD, [ ] );
SetExtraLorentzInvariants( LD, [ s12, s14 ] );
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <An unevaluated non-zero 4 x 4 matrix over an external ring>,
#!   <An unevaluated non-zero 4 x 4 matrix over an external ring> ]
Display( E12[1] );
#! 2*D1,     D1+D2,     D1+D3-s12, D1+D4,
#! D1-D2,    D1-D2,     D1-D2-s12, D1-D2+s14,
#! D2-D3+s12,D2-D3,     D2-D3,     D2-D3-s14,
#! -D1+D4,   -D1+D4-s14,-D1+D4+s12,-D1+D4
S := SyzygiesOfRows( E12[1], E12[2] );
#! <A non-zero 8 x 4 matrix over an external ring>
Sred := ReducedBasisOfRowModule( S );
#! <A non-zero 6 x 4 matrix over an external ring>
Display( EntriesOfHomalgMatrixAsListList( Sred{[ 1 .. 3 ]} ) );
#! [ [ D2-D4, D4, 0, D2 ],
#!   [ D1-D3, -D1, -D1, 0 ],
#!   [ 2*D3*D4-2*D4^2-D4*s12, D2*D4-D3*D4+2*D4^2+D4*s12, D1*D4+D2*D4, 2*D2*D4 ] ]
#! @EndExample
