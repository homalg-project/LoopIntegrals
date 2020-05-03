#! @Chunk Kite

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1..2", "p", 4 );
#! <A loop diagram with loop momenta [ l1, l2 ] & external momenta [ p ]>
s := p^2;;
SetAbbreviation( s, "s" );
SetRelationsOfMomenta( LD, [ ] );
SetIndependentLorentzInvariants( LD, [ l1^2, l2^2, l1*l2, l1*p, l2*p, s ] );
SetPropagators( LD, -[ l1^2, (l1+p)^2, l2^2, (l2-p)^2, (l1+l2)^2 ] );
SetNumerators( LD, -[ ] );
SetExtraLorentzInvariants( LD, [ s ] );
e12 := PairOfMatricesOfLoopDiagramInIndependentLorentzInvariants( LD );
#! [ <An unevaluated non-zero 6 x 5 matrix over an external ring>,
#!   <An unevaluated non-zero 5 x 5 matrix over an external ring> ]
Display( e12[1] );
#! -2*l1l1,-2*l1l1-2*l1p,0,      0,            -2*l1l1-2*l1l2,
#! -2*l1l2,-2*l1l2-2*l2p,0,      0,            -2*l2l2-2*l1l2,
#! -2*l1p, -2*l1p-2*s,   0,      0,            -2*l1p-2*l2p,
#! 0,      0,            -2*l1l2,-2*l1l2+2*l1p,-2*l1l1-2*l1l2,
#! 0,      0,            -2*l2l2,-2*l2l2+2*l2p,-2*l2l2-2*l1l2,
#! 0,      0,            -2*l2p, -2*l2p+2*s,   -2*l1p-2*l2p
Display( e12[2] );
#! -l1l1,0,            0,    0,            0,
#! 0,    -l1l1-2*l1p-s,0,    0,            0,
#! 0,    0,            -l2l2,0,            0,
#! 0,    0,            0,    -l2l2+2*l2p-s,0,
#! 0,    0,            0,    0,            -l1l1-l2l2-2*l1l2
e12 := PairOfMatricesOfLoopDiagramInIndependentLorentzInvariants( LD :
               abbreviation := false );
#! [ <An unevaluated non-zero 6 x 5 matrix over an external ring>,
#!   <An unevaluated non-zero 5 x 5 matrix over an external ring> ]
Display( e12[1] );
#! -2*x1,-2*x1-2*x4,0,    0,         -2*x1-2*x3,
#! -2*x3,-2*x3-2*x5,0,    0,         -2*x2-2*x3,
#! -2*x4,-2*x4-2*x6,0,    0,         -2*x4-2*x5,
#! 0,    0,         -2*x3,-2*x3+2*x4,-2*x1-2*x3,
#! 0,    0,         -2*x2,-2*x2+2*x5,-2*x2-2*x3,
#! 0,    0,         -2*x5,-2*x5+2*x6,-2*x4-2*x5
Display( e12[2] );
#! -x1,0,          0,  0,          0,
#! 0,  -x1-2*x4-x6,0,  0,          0,
#! 0,  0,          -x2,0,          0,
#! 0,  0,          0,  -x2+2*x5-x6,0,
#! 0,  0,          0,  0,          -x1-x2-2*x3
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <A 6 x 5 matrix over an external ring>,
#!   <A 5 x 5 matrix over an external ring> ]
Display( E12[1] );
#! 2*D1,     D1+D2+s,    0,        0,          D1-D3+D5,
#! -D1-D3+D5,-D1-D4+D5-s,0,        0,          -D1+D3+D5,
#! -D1+D2+s, -D1+D2-s,   0,        0,          -D1+D2+D3-D4,
#! 0,        0,          -D1-D3+D5,-D2-D3+D5-s,D1-D3+D5,
#! 0,        0,          2*D3,     D3+D4+s,    -D1+D3+D5,
#! 0,        0,          D3-D4-s,  D3-D4+s,    -D1+D2+D3-D4
Display( E12[2] );
#! D1,0, 0, 0, 0,
#! 0, D2,0, 0, 0,
#! 0, 0, D3,0, 0,
#! 0, 0, 0, D4,0,
#! 0, 0, 0, 0, D5
S := SyzygiesOfRows( E12 );
#! <A non-zero 19 x 6 matrix over an external ring>
Display( EntriesOfHomalgMatrix( S[1] ) );
#! [ 2*D2, 0, 0, D3+D4+s, D1+D2+2*D4-2*D5+s, -D1+D3+D5 ]
Display( EntriesOfHomalgMatrix( S[2] ) );
#! [ 2*D2+D3+D4-2*D5+s, D1+D2+s, -D1+D3-D5, 0, 2*D4, 0 ]
Display( EntriesOfHomalgMatrix( S[3] ) );
#! [ D1-D2, 0, D1, 0, D3-D4, -D3 ]
Display( EntriesOfHomalgMatrix( S[4] ) );
#! [ 0, 0, 0, 0, D4*D5, 0 ]
ibp1 := ShiftOperator( S[1], LD );
#! |[ -2*a2*s-2*a4*s-2*a5*s-2*a2*D1-2*a4*D1-2*a5*D1
#!    -4*a1*D2-2*a2*D2-2*a3*D2-4*a5*D2-2*a3*D4-2*a4*D4-2*a5*D4
#!    +2*a3*D5+2*a4*D5+2*a5*D5
#!    +(D+2)*s+(D+2)*D1+(3*D+2)*D2+(2*D)*D4+(-2*D)*D5 ]|
ViewList( DecomposeInMonomials( ibp1 ) );
#! [ [ |[ -2*a2-2*a4-2*a5+(D+2) ]|, |[ D1 ]| ],
#!   [ |[ -4*a1-2*a2-2*a3-4*a5+(3*D+2) ]|, |[ D2 ]| ],
#!   [ |[ -2*a3-2*a4-2*a5+(2*D) ]|, |[ D4 ]| ],
#!   [ |[ 2*a3+2*a4+2*a5+(-2*D) ]|, |[ D5 ]| ],
#!   [ |[ -2*a2*s-2*a4*s-2*a5*s+(D+2)*s ]|, |[ 1 ]| ] ]
IBPRelation( S[1], LD, [ 1, 1, 1, 1, 1 ] );
#! |[ (D-4)*s+(D-4)*D1+(3*D-10)*D2+(2*D-6)*D4+(-2*D+6)*D5 ]|
ibp2 := ShiftOperator( S[2], LD );
#! |[ -2*a2*s-2*a4*s-2*a5*s-2*a1*D2-2*a2*D2-2*a5*D2
#!    -2*a2*D3-2*a4*D3-2*a5*D3-2*a1*D4-4*a3*D4-2*a4*D4-4*a5*D4
#!    +2*a1*D5+2*a2*D5+2*a5*D5
#!    +(D+2)*s+(2*D)*D2+(D+2)*D3+(3*D+2)*D4+(-2*D)*D5 ]|
ViewList( DecomposeInMonomials( ibp2 ) );
#! [ [ |[ -2*a1-2*a2-2*a5+(2*D) ]|, |[ D2 ]| ],
#!   [ |[ -2*a2-2*a4-2*a5+(D+2) ]|, |[ D3 ]| ],
#!   [ |[ -2*a1-4*a3-2*a4-4*a5+(3*D+2) ]|, |[ D4 ]| ],
#!   [ |[ 2*a1+2*a2+2*a5+(-2*D) ]|, |[ D5 ]| ],
#!   [ |[ -2*a2*s-2*a4*s-2*a5*s+(D+2)*s ]|, |[ 1 ]| ] ]
IBPRelation( S[2], LD, [ 1, 1, 1, 1, 1 ] );
#! |[ (D-4)*s+(2*D-6)*D2+(D-4)*D3+(3*D-10)*D4+(-2*D+6)*D5 ]|
ibp3 := ShiftOperator( S[3], LD );
#! |[ -a1*s+a2*s-a3*s+a4*s-a1*D1-a2*D1-a5*D1+a1*D2+a2*D2+a5*D2
#!    -a3*D3-a4*D3-a5*D3+a3*D4+a4*D4+a5*D4+(D)*D1+(-D)*D2+(D)*D3+(-D)*D4 ]|
ViewList( DecomposeInMonomials( ibp3 ) );
#! [ [ |[ -a1-a2-a5+(D) ]|, |[ D1 ]| ],
#!   [ |[ a1+a2+a5+(-D) ]|, |[ D2 ]| ],
#!   [ |[ -a3-a4-a5+(D) ]|, |[ D3 ]| ],
#!   [ |[ a3+a4+a5+(-D) ]|, |[ D4 ]| ],
#!   [ |[ -a1*s+a2*s-a3*s+a4*s ]|, |[ 1 ]| ] ]
IBPRelation( S[3], LD, [ 1, 1, 1, 1, 1 ] );
#! |[ (D-3)*D1+(-D+3)*D2+(D-3)*D3+(-D+3)*D4 ]|
ibp4 := ShiftOperator( S[4], LD );
#! |[ a5*D1*D4-a5*D3*D4-a4*s*D5-a4*D3*D5-2*a3*D4*D5-a4*D4*D5-a5*D4*D5
#!    -D1*D4+D3*D4+s*D5+D3*D5+(D+2)*D4*D5 ]|
ViewList( DecomposeInMonomials( ibp4 ) );
#! [ [ |[ a5-1 ]|, |[ D1*D4 ]| ],
#!   [ |[ -a5+1 ]|, |[ D3*D4 ]| ],
#!   [ |[ -a4+1 ]|, |[ D3*D5 ]| ],
#!   [ |[ -2*a3-a4-a5+(D+2) ]|, |[ D4*D5 ]| ],
#!   [ |[ -a4*s+s ]|, |[ D5 ]| ] ]
IBPRelation( S[4], LD, [ 1, 1, 1, 1, 1 ] );
#! |[ (D-2)*D4*D5 ]|
#! @EndExample
