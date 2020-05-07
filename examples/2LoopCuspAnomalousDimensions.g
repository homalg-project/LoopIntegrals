#! @Chunk 2LoopCuspAnomalousDimensions

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "k1, k2", "v1, v2", 2 );
#! <A loop diagram with loop momenta [ k1, k2 ] & external momenta [ v1, v2 ]>
cos:= v1*v2;;
SetAbbreviation( cos, "cos" );
rel := [ v1^2 - 1, v2^2 - 1 ];;
SetRelationsOfMomenta( LD, rel );
SetIndependentLorentzInvariants( LD,
        [ k1^2, k1*k2, k1*v1, k1*v2, k2^2, k2*v1, k2*v2, cos ] );
SetPropagators( LD,
        -[ 2*k2*v1 - 1, 2*k2*v2 - 1, (k1 - k2)^2,
           2*k1*v1 - 1, 2*k1*v2 - 1, k1^2 ] );
SetNumerators( LD, -[ k2^2 ] );
SetExtraLorentzInvariants( LD, [ cos ] );
e12 := PairOfMatricesOfLoopDiagramInPropagators( LD : rational := true );
#! [ <A 8 x 6 matrix over an external ring>,
#!   <A 6 x 6 matrix over an external ring> ]
Display( e12[1] );
#! 0,       0,       D3+D6-N7, D4-1,    D5-1,    2*D6,
#! 0,       0,       -D3+D6-N7,D1-1,    D2-1,    -D3+D6+N7,
#! 0,       0,       -D1+D4,   -2,      (-2*cos),D4-1,
#! 0,       0,       -D2+D5,   (-2*cos),-2,      D5-1,
#! D4-1,    D5-1,    -D3-D6+N7,0,       0,       0,
#! D1-1,    D2-1,    D3-D6+N7, 0,       0,       0,
#! -2,      (-2*cos),D1-D4,    0,       0,       0,
#! (-2*cos),-2,      D2-D5,    0,       0,       0
s := SyzygiesOfRows( e12 );
#! <A non-zero 46 x 8 matrix over an external ring>
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <A 8 x 6 matrix over an external ring>,
#!   <A 6 x 6 matrix over an external ring> ]
Display( E12[1] );
#! 0,     0,     D3+D6-N7, D4-1,  D5-1,  2*D6,
#! 0,     0,     -D3+D6-N7,D1-1,  D2-1,  -D3+D6+N7,
#! 0,     0,     -D1+D4,   -2,    -2*cos,D4-1,
#! 0,     0,     -D2+D5,   -2*cos,-2,    D5-1,
#! D4-1,  D5-1,  -D3-D6+N7,0,     0,     0,
#! D1-1,  D2-1,  D3-D6+N7, 0,     0,     0,
#! -2,    -2*cos,D1-D4,    0,     0,     0,
#! -2*cos,-2,    D2-D5,    0,     0,     0
S := SyzygiesOfRows( E12 );
#! <A non-zero 79 x 8 matrix over an external ring>
EntriesOfHomalgMatrix( S[1] );
#! [ -2*D4*D6-2*D5*D6, 0, D5*D6, D4*D6,
#!   2*cos*D6-2*D6, -2*cos*D6-2*D4*D6-2*D5*D6+2*D6, D5*D6, D4*D6 ]
EntriesOfHomalgMatrix( S[2] );
#! [ -2*D1*D6-2*D2*D6+4*D6, 2*cos*D6-2*D6, D2*D6-D6, D1*D6-D6,
#!   0, 2*cos*D6-2*D1*D6-2*D2*D6+2*D6, D2*D6-D6, D1*D6-D6 ]
EntriesOfHomalgMatrix( S[3] );
#! [ 2*cos*D6+2*D6, 0, -D6, -D6,
#!   0, 2*cos*D6+2*D6, -D6, -D6 ]
#! @EndExample
