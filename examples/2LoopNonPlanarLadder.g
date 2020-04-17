#! @Chunk 2LoopNonPlanarLadder

LoadPackage( "LoopIntegrals" );

#! @Example
LOOP_INTEGRALS.Dimension := 2;
#! 2
LD := LoopDiagram( "l1..2", "p1..2" );
#! <A loop diagram with loop momenta [ l1, l2 ] & external momenta [ p1, p2 ]>
SetPropagators( LD,
        [ (l1-p1)^2, (l1+p2)^2, l2^2, (l2-p1)^2,
          (l1-l2)^2, (l1-l2+p2)^2 ] );
SetNumerators( LD, [ l1^2 ] );
rel := List( ExternalMomenta( LD ), p -> p^2 );;
SetRelationsOfMomenta( LD, rel );
SetIndependetLorentzInvariants( LD,
        [ l1^2, l2^2, l1*l2,
          l1*p1, l2*p1, l1*p2, l2*p2,
          p1*p2 ] );
SetExtraLorentzInvariants( LD, [ p1*p2 ] );
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <A non-zero 8 x 6 matrix over an external ring>,
#!   <A non-zero 6 x 6 matrix over an external ring> ]
Display( E12[1] );
#! D1+N7,     D2+N7,      0,          0,       -D3+D5+N7,   D2-D3+D5,
#! D4-D5+N7,  D2+D3-D6,   0,          0,       -D3-D5+N7,   D2-D3-D6,
#! -D1+N7,    -D1+N7+2*x8,0,          0,       -D1-D3+D4+N7,_[3,6],
#! D2-N7-2*x8,D2-N7,      0,          0,       -D5+D6,      -D5+D6,
#! 0,         0,          D3-D5+N7,   D1+D3-D5,D3-D5-N7,    -D2+D3-D5,
#! 0,         0,          2*D3,       D3+D4,   D3+D5-N7,    -D2+D3+D6,
#! 0,         0,          D3-D4,      D3-D4,   D1+D3-D4-N7, _[7,6],
#! 0,         0,          D2+D5-D6-N7,_[8,4],  D5-D6,       D5-D6
S := SyzygiesOfRows( E12[1], E12[2] );
#! <A non-zero 52 x 8 matrix over an external ring>
#! @EndExample
