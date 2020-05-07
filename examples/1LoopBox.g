#! @Chunk 1LoopBox

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1", "k1..2,k4", 2 );
#! <A loop diagram with loop momenta [ l1 ] & external momenta [ k1, k2, k4 ]>
s12 := 2*k1*k2;;
SetAbbreviation( s12, "s12" );
s14 := 2*k1*k4;;
SetAbbreviation( s14, "s14" );
rel1 := List( ExternalMomenta( LD ), k -> k^2 );;
rel2 := [ (k1+k2+k4)^2 ];;
SetRelationsOfMomenta( LD, Concatenation( rel1, rel2 ) );
SetIndependentLorentzInvariants( LD,
        [ l1^2, l1*k1, l1*k2, l1*k4, s12, s14 ] );
SetPropagators( LD, -[ l1^2, (l1-k1)^2, (l1-k1-k2)^2, (l1+k4)^2 ] );
SetNumerators( LD, -[ ] );
SetExtraLorentzInvariants( LD, [ s12, s14 ] );
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <A 4 x 4 matrix over an external ring>,
#!   <A 4 x 4 matrix over an external ring> ]
Display( E12[1] );
#! 2*D1,      D1+D2,    s12+D1+D3, D1+D4,
#! D1-D2,     D1-D2,    s12+D1-D2, -s14+D1-D2,
#! -s12+D2-D3,D2-D3,    D2-D3,     s14+D2-D3,
#! -D1+D4,    s14-D1+D4,-s12-D1+D4,-D1+D4
S := SyzygiesOfRows( E12 );
#! <A non-zero 12 x 4 matrix over an external ring>
Sred := ReducedBasisOfRowModule( S );
#! <A non-zero 6 x 4 matrix over an external ring>
Display( Sred );
#! D2-D4,                D4,                        0,          D2,
#! D1-D3,                -D1,                       -D1,        0,
#! s12*D4+2*D3*D4-2*D4^2,-s12*D4+D2*D4-D3*D4+2*D4^2,D1*D4+D2*D4,2*D2*D4,
#! s14*D3-2*D3^2+2*D3*D4,-D1*D3-D3*D4,              -2*D1*D3,   -D1*D3-D2*D3,
#! -D3*D4^2,             D3*D4^2,                   0,          D2*D3*D4,
#! D3^2*D4,              0,                         D1*D3*D4,   0
Display( EntriesOfHomalgMatrixAsListList( Sred{[ 1 .. 3 ]} ) );
#! [ [ D2-D4, D4, 0, D2 ],
#!   [ D1-D3, -D1, -D1, 0 ],
#!   [ s12*D4+2*D3*D4-2*D4^2, -s12*D4+D2*D4-D3*D4+2*D4^2, D1*D4+D2*D4, 2*D2*D4 ] ]
Sibp1 := ShiftOperator( S[1], LD );
#! |[ -s14*a2+s14*a4-a1*D2-a2*D2-a3*D2-a4*D2
#!    +a1*D4+a2*D4+a3*D4+a4*D4+(D)*D2+(-D)*D4 ]|
ViewList( DecomposeInMonomials( Sibp1 ) );
#! [ [ |[ -a1-a2-a3-a4+(D) ]|, |[ D2 ]| ],
#!   [ |[ a1+a2+a3+a4+(-D) ]|, |[ D4 ]| ],
#!   [ |[ -s14*a2+s14*a4 ]|, |[ 1 ]| ] ]
sibp1 := IBPRelation( S[1], LD, [ 1, 1, 1, 1 ] );
#! |[ (D-4)*D2+(-D+4)*D4 ]|
Sibp2 := ShiftOperator( S[2], LD );
#! |[ -s12*a1+s12*a3-a1*D1-a2*D1-a3*D1-a4*D1
#!    +a1*D3+a2*D3+a3*D3+a4*D3+(D)*D1+(-D)*D3 ]|
ViewList( DecomposeInMonomials( Sibp2 ) );
#! [ [ |[ -a1-a2-a3-a4+(D) ]|, |[ D1 ]| ],
#!   [ |[ a1+a2+a3+a4+(-D) ]|, |[ D3 ]| ],
#!   [ |[ -s12*a1+s12*a3 ]|, |[ 1 ]| ] ]
sibp2 := IBPRelation( S[2], LD, [ 1, 1, 1, 1 ] );
#! |[ (D-4)*D1+(-D+4)*D3 ]|
Sibp3 := ShiftOperator( S[3], LD );;
ViewList( DecomposeInMonomials( Sibp3 ) );
#! [ [ |[ -2*a1-2*a2-2*a3-2*a4+(2*D+2) ]|, |[ D3*D4 ]| ],
#!   [ |[ 2*a1+2*a2+2*a3+2*a4+(-2*D-2) ]|, |[ D4^2 ]| ],
#!   [ |[ -s14*a4+s14 ]|, |[ D1 ]| ],
#!   [ |[ -s12*a4+s12 ]|, |[ D2 ]| ],
#!   [ |[ -s14*a4+s14 ]|, |[ D3 ]| ],
#!   [ |[ -2*s12*a2-2*s14*a2-2*s12*a3-s12*a4+2*s14*a4+(D+1)*s12-2*s14 ]|,
#!   |[ D4 ]| ],
#!   [ |[ -s12*s14*a4+s12*s14 ]|, |[ 1 ]| ] ]
sibp3 := IBPRelation( S[3], LD, [ 1, 1, 1, 1 ] );
#! |[ (D-4)*s12*D4-2*s14*D4+(2*D-6)*D3*D4+(-2*D+6)*D4^2 ]|
#! @EndExample
