#! @Chunk Kite

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1..2", "p", 4 );
#! <A loop diagram with loop momenta [ l1, l2 ] & external momenta [ p ]>
s := p^2;;
SetAbbreviation( s, "s" );
SetRelationsOfExternalMomenta( LD, [ ] );
SetIndependentLorentzInvariants( LD, [ l1^2, l2^2, l1*l2, l1*p, l2*p, s ] );
SetPropagators( LD, -[ l1^2, (l1+p)^2, l2^2, (l2-p)^2, (l1+l2)^2 ] );
SetNumerators( LD, -[ ] );
SetExtraLorentzInvariants( LD, [ s ] );
e12 := PairOfMatricesOfLoopDiagramInIndependentLorentzInvariants( LD );
#! [ <An unevaluated non-zero 5 x 6 matrix over an external ring>,
#!   <An unevaluated non-zero 5 x 5 matrix over an external ring> ]
Display( e12[1] );
#! -2*l1l1,-2*l1l2,-2*l1p,      0,      0,      0,
#! _[2,1], _[2,2], -2*l1p-2*s,  0,      0,      0,
#! 0,      0,      0,           -2*l1l2,-2*l2l2,-2*l2p,
#! 0,      0,      0,           _[4,4], _[4,5], -2*l2p+2*s,
#! _[5,1], _[5,2], -2*l1p-2*l2p,_[5,4], _[5,5], -2*l1p-2*l2p
Display( e12[2] );
#! -l1l1,0,            0,    0,            0,
#! 0,    -l1l1-2*l1p-s,0,    0,            0,
#! 0,    0,            -l2l2,0,            0,
#! 0,    0,            0,    -l2l2+2*l2p-s,0,
#! 0,    0,            0,    0,            -l1l1-l2l2-2*l1l2
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <A 5 x 6 matrix over an external ring>,
#!   <A 5 x 5 matrix over an external ring> ]
Display( E12[1] );
#! 2*D1,    -D1-D3+D5,  s-D1+D2,     0,          0,        0,
#! s+D1+D2, -s-D1-D4+D5,-s-D1+D2,    0,          0,        0,
#! 0,       0,          0,           -D1-D3+D5,  2*D3,     -s+D3-D4,
#! 0,       0,          0,           -s-D2-D3+D5,s+D3+D4,  s+D3-D4,
#! D1-D3+D5,-D1+D3+D5,  -D1+D2+D3-D4,D1-D3+D5,   -D1+D3+D5,-D1+D2+D3-D4
Display( E12[2] );
#! D1,0, 0, 0, 0,
#! 0, D2,0, 0, 0,
#! 0, 0, D3,0, 0,
#! 0, 0, 0, D4,0,
#! 0, 0, 0, 0, D5
S := SyzygiesOfColumns( E12 );
#! <A non-zero 6 x 16 matrix over an external ring>
Display( EntriesOfHomalgMatrix( CertainColumns( S, [ 1 ] ) ) );
#! [ D1-D2, 0, D1, 0, D3-D4, -D3 ]
Display( EntriesOfHomalgMatrix( CertainColumns( S, [ 2 ] ) ) );
#! [ 2*D2, 0, 0, s+D3+D4, s+D1+D2+2*D4-2*D5, -D1+D3+D5 ]
Display( EntriesOfHomalgMatrix( CertainColumns( S, [ 3 ] ) ) );
#! [ s+2*D2+D3+D4-2*D5, s+D1+D2, -D1+D3-D5, 0, 2*D4, 0 ]
Display( EntriesOfHomalgMatrix( CertainColumns( S, [ 4 ] ) ) );
#! [ 0, 0, 0, 0, D4*D5, 0 ]
Sibp1 := IBPRelation( CertainColumns( S, [ 1 ] ), LD );
#! |[ -s*a1+s*a2-s*a3+s*a4+d*D1-a1*D1-a2*D1-a5*D1-d*D2+a1*D2+a2*D2+a5*D2\
#!    +d*D3-a3*D3-a4*D3-a5*D3-d*D4+a3*D4+a4*D4+a5*D4 ]|
ViewList( DecomposeInMonomials( Sibp1 ) );
#! [ [ |[ d-a1-a2-a5 ]|, |[ D1 ]| ],
#!   [ |[ -d+a1+a2+a5 ]|, |[ D2 ]| ],
#!   [ |[ d-a3-a4-a5 ]|, |[ D3 ]| ],
#!   [ |[ -d+a3+a4+a5 ]|, |[ D4 ]| ],
#!   [ |[ -s*a1+s*a2-s*a3+s*a4 ]|, |[ 1 ]| ] ]
sibp1 := IBPRelation( CertainColumns( S, [ 1 ] ), LD, [ 1, 1, 1, 1, 1 ] );
#! |[ d*D1-d*D2+d*D3-d*D4-3*D1+3*D2-3*D3+3*D4 ]|
ViewList( DecomposeInMonomials( sibp1 ) );
#! [ [ |[ d-3 ]|, |[ D1 ]| ],
#!   [ |[ -d+3 ]|, |[ D2 ]| ],
#!   [ |[ d-3 ]|, |[ D3 ]| ],
#!   [ |[ -d+3 ]|, |[ D4 ]| ] ]
Sibp2 := IBPRelation( CertainColumns( S, [ 2 ] ), LD );
#! |[ d*s-2*s*a2-2*s*a4-2*s*a5+d*D1-2*a2*D1-2*a4*D1-2*a5*D1\
#!    +3*d*D2-4*a1*D2-2*a2*D2-2*a3*D2-4*a5*D2+2*d*D4-2*a3*D4-2*a4*D4-2*a5*D4\
#!    -2*d*D5+2*a3*D5+2*a4*D5+2*a5*D5+2*s+2*D1+2*D2 ]|
ViewList( DecomposeInMonomials( Sibp2 ) );
#! [ [ |[ d-2*a2-2*a4-2*a5+2 ]|, |[ D1 ]| ],
#!   [ |[ 3*d-4*a1-2*a2-2*a3-4*a5+2 ]|, |[ D2 ]| ],
#!   [ |[ 2*d-2*a3-2*a4-2*a5 ]|, |[ D4 ]| ],
#!   [ |[ -2*d+2*a3+2*a4+2*a5 ]|, |[ D5 ]| ],
#!   [ |[ d*s-2*s*a2-2*s*a4-2*s*a5+2*s ]|, |[ 1 ]| ] ]
sibp2 := IBPRelation( CertainColumns( S, [ 2 ] ), LD, [ 1, 1, 1, 1, 1 ] );
#! |[ d*s+d*D1+3*d*D2+2*d*D4-2*d*D5-4*s-4*D1-10*D2-6*D4+6*D5 ]|
ViewList( DecomposeInMonomials( sibp2 ) );
#! [ [ |[ d-4 ]|, |[ D1 ]| ],
#!   [ |[ 3*d-10 ]|, |[ D2 ]| ],
#!   [ |[ 2*d-6 ]|, |[ D4 ]| ],
#!   [ |[ -2*d+6 ]|, |[ D5 ]| ],
#!   [ |[ d*s-4*s ]|, |[ 1 ]| ] ]
Sibp3 := IBPRelation( CertainColumns( S, [ 3 ] ), LD );
#! |[ d*s-2*s*a2-2*s*a4-2*s*a5+2*d*D2-2*a1*D2-2*a2*D2-2*a5*D2\
#!    +d*D3-2*a2*D3-2*a4*D3-2*a5*D3+3*d*D4-2*a1*D4-4*a3*D4-2*a4*D4-4*a5*D4\
#!    -2*d*D5+2*a1*D5+2*a2*D5+2*a5*D5+2*s+2*D3+2*D4 ]|
ViewList( DecomposeInMonomials( Sibp3 ) );
#! [ [ |[ 2*d-2*a1-2*a2-2*a5 ]|, |[ D2 ]| ],
#!   [ |[ d-2*a2-2*a4-2*a5+2 ]|, |[ D3 ]| ],
#!   [ |[ 3*d-2*a1-4*a3-2*a4-4*a5+2 ]|, |[ D4 ]| ],
#!   [ |[ -2*d+2*a1+2*a2+2*a5 ]|, |[ D5 ]| ],
#!   [ |[ d*s-2*s*a2-2*s*a4-2*s*a5+2*s ]|, |[ 1 ]| ] ]
sibp3 := IBPRelation( CertainColumns( S, [ 3 ] ), LD, [ 1, 1, 1, 1, 1 ] );
#! |[ d*s+2*d*D2+d*D3+3*d*D4-2*d*D5-4*s-6*D2-4*D3-10*D4+6*D5 ]|
ViewList( DecomposeInMonomials( sibp3 ) );
#! [ [ |[ 2*d-6 ]|, |[ D2 ]| ],
#!   [ |[ d-4 ]|, |[ D3 ]| ],
#!   [ |[ 3*d-10 ]|, |[ D4 ]| ],
#!   [ |[ -2*d+6 ]|, |[ D5 ]| ],
#!   [ |[ d*s-4*s ]|, |[ 1 ]| ] ]
Sibp4 := IBPRelation( CertainColumns( S, [ 4 ] ), LD );
#! |[ a5*D1*D4-a5*D3*D4-s*a4*D5-a4*D3*D5+d*D4*D5-2*a3*D4*D5-a4*D4*D5-a5*D4*D5\
#!    -D1*D4+D3*D4+s*D5+D3*D5+2*D4*D5 ]|
ViewList( DecomposeInMonomials( Sibp4 ) );
#! [ [ |[ a5-1 ]|, |[ D1*D4 ]| ],
#!   [ |[ -a5+1 ]|, |[ D3*D4 ]| ],
#!   [ |[ -a4+1 ]|, |[ D3*D5 ]| ],
#!   [ |[ d-2*a3-a4-a5+2 ]|, |[ D4*D5 ]| ],
#!   [ |[ -s*a4+s ]|, |[ D5 ]| ] ]
sibp4 := IBPRelation( CertainColumns( S, [ 4 ] ), LD, [ 1, 1, 1, 1, 1 ] );
#! |[ d*D4*D5-2*D4*D5 ]|
ViewList( DecomposeInMonomials( sibp4 ) );
#! [ [ |[ d-2 ]|, |[ D4*D5 ]| ] ]
#! @EndExample
