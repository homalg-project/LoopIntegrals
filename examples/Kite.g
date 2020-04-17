#! @Chunk Kite

LoadPackage( "LoopIntegrals" );

#! @Example
LOOP_INTEGRALS.Dimension := 4;
#! 4
LD := LoopDiagram( "l1..2", "p" );
#! <A loop diagram with loop momenta [ l1, l2 ] & external momenta [ p ]>
SetPropagators( LD, [ l1^2, (l1+p)^2, l2^2, (l2-p)^2, (l1+l2)^2 ] );
SetNumerators( LD, [ ] );
SetRelationsOfMomenta( LD, [ ] );
SetIndependetLorentzInvariants( LD, [ l1^2, l2^2, l1*l2, l1*p, l2*p, p^2 ] );
SetExtraLorentzInvariants( LD, [ p^2 ] );
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <A non-zero 6 x 5 matrix over an external ring>,
#!   <A non-zero 5 x 5 matrix over an external ring> ]
Display( E12[1] );
#! 2*D1,     D1+D2-x6,    0,        0,           D1-D3+D5,
#! -D1-D3+D5,-D1-D4+D5+x6,0,        0,           -D1+D3+D5,
#! -D1+D2-x6,-D1+D2+x6,   0,        0,           -D1+D2+D3-D4,
#! 0,        0,           -D1-D3+D5,-D2-D3+D5+x6,D1-D3+D5,
#! 0,        0,           2*D3,     D3+D4-x6,    -D1+D3+D5,
#! 0,        0,           D3-D4+x6, D3-D4-x6,    -D1+D2+D3-D4
Display( E12[2] );
#! D1,0, 0, 0, 0,
#! 0, D2,0, 0, 0,
#! 0, 0, D3,0, 0,
#! 0, 0, 0, D4,0,
#! 0, 0, 0, 0, D5
S := SyzygiesOfRows( E12[1], E12[2] );
#! <A non-zero 19 x 6 matrix over an external ring>
Display( EntriesOfHomalgMatrix( CertainRows( S, [ 1 ] ) ) );
#! [ 2*D2, 0, 0, D3+D4-x6, D1+D2+2*D4-2*D5-x6, -D1+D3+D5 ]
Display( EntriesOfHomalgMatrix( CertainRows( S, [ 2 ] ) ) );
#! [ 2*D2+D3+D4-2*D5-x6, D1+D2-x6, -D1+D3-D5, 0, 2*D4, 0 ]
Display( EntriesOfHomalgMatrix( CertainRows( S, [ 3 ] ) ) );
#! [ D1-D2, 0, D1, 0, D3-D4, -D3 ]
Display( EntriesOfHomalgMatrix( CertainRows( S, [ 4 ] ) ) );
#! [ 0, 0, 0, 0, D4*D5, 0 ]
#! @EndExample
