#! @Chunk 2LoopNonPlanarLadder

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1..2", "p1..2", 2 );
#! <A loop diagram with loop momenta [ l1, l2 ] & external momenta [ p1, p2 ]>
s := 2*p1*p2;;
SetAbbreviation( s, "s" );
rel := List( ExternalMomenta( LD ), p -> p^2 );;
SetRelationsOfMomenta( LD, rel );
SetIndependentLorentzInvariants( LD,
        [ l1^2, l2^2, l1*l2,
          l1*p1, l2*p1, l1*p2, l2*p2,
          s ] );
SetPropagators( LD,
        -[ (l1-p1)^2, (l1+p2)^2, l2^2, (l2-p1)^2,
           (l1-l2)^2, (l1-l2+p2)^2 ] );
SetNumerators( LD, -[ l1^2 ] );
SetExtraLorentzInvariants( LD, [ s ] );
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <A 8 x 6 matrix over an external ring>,
#!   <A 6 x 6 matrix over an external ring> ]
Display( E12[1] );
#! D1+N7,   D2+N7,   0,          0,            -D3+D5+N7,   D2-D3+D5,
#! D4-D5+N7,D2+D3-D6,0,          0,            -D3-D5+N7,   D2-D3-D6,
#! -D1+N7,  -s-D1+N7,0,          0,            -D1-D3+D4+N7,-s-D1-D3+D4+N7,
#! s+D2-N7, D2-N7,   0,          0,            -D5+D6,      -D5+D6,
#! 0,       0,       D3-D5+N7,   D1+D3-D5,     D3-D5-N7,    -D2+D3-D5,
#! 0,       0,       2*D3,       D3+D4,        D3+D5-N7,    -D2+D3+D6,
#! 0,       0,       D3-D4,      D3-D4,        D1+D3-D4-N7, s+D1+D3-D4-N7,
#! 0,       0,       D2+D5-D6-N7,s+D2+D5-D6-N7,D5-D6,       D5-D6
S := SyzygiesOfRows( E12 );
#! <A non-zero 58 x 8 matrix over an external ring>
Display( EntriesOfHomalgMatrix( S[1] ) );
#! [ D1-D2, 0, D2, D1,
#!  -D3+D4, D1-D2+D3-D4-D5+D6, D2-D6, D4 ]
Display( EntriesOfHomalgMatrix( S[2] ) );
#! [ s+2*D2-2*N7, 0, -D2+N7, -D1-N7,
#!   0, s+2*D2+2*D5-2*D6-2*N7, -D2-D5+D6+N7, -D3-D4 ]
Display( EntriesOfHomalgMatrix( S[3] ) );
#! [ 0, 0, 0, 0,
#!   -D3*D6+D4*D6, D1*D6+2*D3*D6-2*D4*D6-D6*N7, -D3*D6-D5*D6+D6*N7, 0 ]
#! @EndExample
