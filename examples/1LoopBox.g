#! @Chunk 1LoopBox

LoadPackage( "LoopIntegrals" );

#! @Example
LOOP_INTEGRALS.Dimension := 2;
#! 2
LD := LoopDiagram( "l1", "k1..2,k4" );
#! <A loop diagram with loop momenta [ l1 ] & external momenta [ k1, k2, k4 ]>
SetPropagators( LD, [ l1^2, (l1-k1)^2, (l1-k1-k2)^2, (l1+k4)^2 ] );
SetNumerators( LD, [ ] );
rel1 := List( ExternalMomenta( LD ), k -> k^2 );;
rel2 := [ (k1+k2+k4)^2 ];;
SetRelationsOfMomenta( LD, Concatenation( rel1, rel2 ) );
SetIndependetLorentzInvariants( LD,
        [ l1^2, l1*k1, l1*k2, l1*k4, 2*k1*k2, 2*k1*k4 ] );
SetExtraLorentzInvariants( LD, [ 2*k1*k2, 2*k1*k4 ] );
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <A non-zero 4 x 4 matrix over an external ring>,
#!   <A non-zero 4 x 4 matrix over an external ring> ]
Display( E12[1] );
#! 2*D1,    D1+D2,    D1+D3-x5, D1+D4,
#! D1-D2,   D1-D2,    D1-D2-x5, D1-D2+x6,
#! D2-D3+x5,D2-D3,    D2-D3,    D2-D3-x6,
#! -D1+D4,  -D1+D4-x6,-D1+D4+x5,-D1+D4
S := SyzygiesOfRows( E12[1], E12[2] );
#! <A non-zero 8 x 4 matrix over an external ring>
Sred := ReducedBasisOfRowModule( S );
#! <A non-zero 6 x 4 matrix over an external ring>
Display( EntriesOfHomalgMatrixAsListList( Sred{[ 1 .. 3 ]} ) );
#! [ [ D2-D4, D4, 0, D2 ],
#!   [ D1-D3, -D1, -D1, 0 ],
#!   [ 2*D3*D4-2*D4^2-D4*x5, D2*D4-D3*D4+2*D4^2+D4*x5, D1*D4+D2*D4, 2*D2*D4 ] ]
#! @EndExample
