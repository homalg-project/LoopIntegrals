#! @Chunk 2LoopNonPlanarLadder

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1..2", "p1..2" );
#! <A loop diagram with loop momenta [ l1, l2 ] & external momenta [ p1, p2 ]>
s := 2*p1*p2;;
SetAbbreviation( s, "s" );
rel := List( ExternalMomenta( LD ), p -> p^2 );;
SetRelationsOfExternalMomenta( LD, rel );
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
#! [ <A 6 x 8 matrix over an external ring>,
#!   <A 6 x 6 matrix over an external ring> ]
Display( E12[1] );
#! D1+N7,    D4-D5+N7, -D1+N7,  s+D2-N7,0,        0,        0,     0,
#! D2+N7,    D2+D3-D6, -s-D1+N7,D2-N7,  0,        0,        0,     0,
#! 0,        0,        0,       0,      D3-D5+N7, 2*D3,     D3-D4, _[3,8],
#! 0,        0,        0,       0,      D1+D3-D5, D3+D4,    D3-D4, _[4,8],
#! -D3+D5+N7,-D3-D5+N7,_[5,3],  -D5+D6, D3-D5-N7, D3+D5-N7, _[5,7],D5-D6,
#! D2-D3+D5, D2-D3-D6, _[6,3],  -D5+D6, -D2+D3-D5,-D2+D3+D6,_[6,7],D5-D6
S := SyzygiesOfColumns( E12 );
#! <A non-zero 8 x 48 matrix over an external ring>
Display( EntriesOfHomalgMatrix( CertainColumns( S, [ 1 ] ) ) );
#! [ D1-D2, 0, D2, D1,
#!  -D3+D4, D1-D2+D3-D4-D5+D6, D2-D6, D4 ]
Display( EntriesOfHomalgMatrix( CertainColumns( S, [ 2 ] ) ) );
#! [ s+2*D2-2*N7, 0, -D2+N7, -D1-N7,
#!   0, s+2*D2+2*D5-2*D6-2*N7, -D2-D5+D6+N7, -D3-D4 ]
Display( EntriesOfHomalgMatrix( CertainColumns( S, [ 3 ] ) ) );
#! [ 0, 0, 0, 0,
#!   -D3*D6+D4*D6, D1*D6+2*D3*D6-2*D4*D6-D6*N7, -D3*D6-D5*D6+D6*N7, 0 ]
#! @EndExample
