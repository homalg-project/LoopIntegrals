#! @Chunk 2LoopBox

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1..2", "k1..2,k4" );
#! <A loop diagram with loop momenta [ l1, l2 ] &
#!  external momenta [ k1, k2, k4 ]>
s12 := 2*k1*k2;;
SetAbbreviation( s12, "s12" );
s14 := 2*k1*k4;;
SetAbbreviation( s14, "s14" );
rel1 := List( ExternalMomenta( LD ), k -> k^2 );;
rel2 := [ (k1+k2+k4)^2 ];;
SetRelationsOfExternalMomenta( LD, Concatenation( rel1, rel2 ) );
SetIndependentLorentzInvariants( LD,
        [ l1^2, l1*l2, l2^2, l1*k1, l1*k2, l1*k4,
          l2*k1, l2*k2, l2*k4, s12, s14 ] );
SetPropagators( LD,
        -[ l1^2, (l1-k1)^2, (l1-k1-k2)^2,
           l2^2, (l2-k4)^2, (l2+k1+k2)^2,
           (l1+l2)^2 ] );
SetNumerators( LD, -[ l1*k4, l2*k1 ] );
SetExtraLorentzInvariants( LD, [ s12, s14 ] );
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <A 7 x 10 matrix over an external ring>,
#!   <A 7 x 7 matrix over an external ring> ]
S := SyzygiesOfColumns( E12 );
#! <A non-zero 10 x 95 matrix over an external ring>
Display( EntriesOfHomalgMatrix( CertainColumns( S, [ 1 ] ) ) );
#! [ -D2+D3+2*N8, 0, -2*N8, D1, -D2, 0, -2*D4+D5+D6-2*N9, -D4+D5, -D4, D4+2*N9 ]
Display( EntriesOfHomalgMatrix( CertainColumns( S, [ 2 ] ) ) );
#! [ D1-D2+2*N8, 0, -D1-2*N8, 0, -D2, 0, -D4+D5-2*N9, D5, 0, D4+2*N9 ]
Display( EntriesOfHomalgMatrix( CertainColumns( S, [ 3 ] ) ) );
#! [ 0, 0, 0, 0, 0, 0, D5*D7-D6*D7-2*D7*N9, D4*D7+D5*D7, D4*D7, D4*D7+2*D7*N9 ]
Sibp1 := IBPRelation( CertainColumns( S, [ 1 ] ), LD );;
ViewList( DecomposeInMonomials( Sibp1 ) );
#! [ [ |[ -1/2*s12*a8-1/2*s14*a8 ]|, |[ D1*N8_ ]| ],
#!   [ |[ -1/2*s12*a9+1/2*s14*a9 ]|, |[ D4*N9_ ]| ],
#!   [ |[ -d+a1+a2+a3+a7+a8 ]|, |[ D2 ]| ],
#!   [ |[ d-a1-a2-a3-a7-a8 ]|, |[ D3 ]| ],
#!   [ |[ -2*d+2*a4+2*a5+2*a6+2*a7+2*a9 ]|, |[ D4 ]| ],
#!   [ |[ d-a4-a5-a6-a7-a9 ]|, |[ D5 ]| ],
#!   [ |[ d-a4-a5-a6-a7-a9 ]|, |[ D6 ]| ],
#!   [ |[ 2*d-2*a1-2*a2-2*a3-2*a7-2*a8 ]|, |[ N8 ]| ],
#!   [ |[ -2*d+2*a4+2*a5+2*a6+2*a7+2*a9 ]|, |[ N9 ]| ],
#!   [ |[ s12*a1+s14*a2-s12*a3+s12*a4-s14*a5-s12*a6-s14*a8+s14*a9 ]|, |[ 1 ]| ] ]
sibp1 := IBPRelation( CertainColumns( S, [ 1 ] ), LD, [ 1, 1, 1, 1, 1, 1, 1, 1, 1 ] );
#! |[ -1/2*s12*D1*N8_-1/2*s14*D1*N8_-1/2*s12*D4*N9_+1/2*s14*D4*N9_-d*D2+d*D3\
#!    -2*d*D4+d*D5+d*D6+2*d*N8-2*d*N9+5*D2-5*D3+10*D4-5*D5-5*D6-10*N8+10*N9 ]|
ViewList( DecomposeInMonomials( sibp1 ) );
#! [ [ |[ -1/2*s12-1/2*s14 ]|, |[ D1*N8_ ]| ],
#!   [ |[ -1/2*s12+1/2*s14 ]|, |[ D4*N9_ ]| ],
#!   [ |[ -d+5 ]|, |[ D2 ]| ],
#!   [ |[ d-5 ]|, |[ D3 ]| ],
#!   [ |[ -2*d+10 ]|, |[ D4 ]| ],
#!   [ |[ d-5 ]|, |[ D5 ]| ],
#!   [ |[ d-5 ]|, |[ D6 ]| ],
#!   [ |[ 2*d-10 ]|, |[ N8 ]| ],
#!   [ |[ -2*d+10 ]|, |[ N9 ]| ] ]
Sibp2 := IBPRelation( CertainColumns( S, [ 2 ] ), LD );;
ViewList( DecomposeInMonomials( Sibp2 ) );
#! [ [ |[ -1/2*s14*a8 ]|, |[ D1*N8_ ]| ],
#!   [ |[ 1/2*s14*a9 ]|, |[ D4*N9_ ]| ],
#!   [ |[ d-a1-a2-a3-a7-a8 ]|, |[ D1 ]| ],
#!   [ |[ -d+a1+a2+a3+a7+a8 ]|, |[ D2 ]| ],
#!   [ |[ -d+a4+a5+a6+a7+a9 ]|, |[ D4 ]| ],
#!   [ |[ d-a4-a5-a6-a7-a9 ]|, |[ D5 ]| ],
#!   [ |[ 2*d-2*a1-2*a2-2*a3-2*a7-2*a8 ]|, |[ N8 ]| ],
#!   [ |[ -2*d+2*a4+2*a5+2*a6+2*a7+2*a9 ]|, |[ N9 ]| ],
#!   [ |[ s14*a2-s14*a5-s14*a8+s14*a9 ]|, |[ 1 ]| ] ]
sibp2 := IBPRelation( CertainColumns( S, [ 2 ] ), LD, [ 1, 1, 1, 1, 1, 1, 1, 1, 1 ] );
#! |[ -1/2*s14*D1*N8_+1/2*s14*D4*N9_+d*D1-d*D2-d*D4+d*D5\
#!    +2*d*N8-2*d*N9-5*D1+5*D2+5*D4-5*D5-10*N8+10*N9 ]|
ViewList( DecomposeInMonomials( sibp2 ) );
#! [ [ |[ -1/2*s14 ]|, |[ D1*N8_ ]| ],
#!   [ |[ 1/2*s14 ]|, |[ D4*N9_ ]| ],
#!   [ |[ d-5 ]|, |[ D1 ]| ],
#!   [ |[ -d+5 ]|, |[ D2 ]| ],
#!   [ |[ -d+5 ]|, |[ D4 ]| ],
#!   [ |[ d-5 ]|, |[ D5 ]| ],
#!   [ |[ 2*d-10 ]|, |[ N8 ]| ],
#!   [ |[ -2*d+10 ]|, |[ N9 ]| ] ]
Sibp3 := IBPRelation( CertainColumns( S, [ 3 ] ), LD );;
ViewList( DecomposeInMonomials( Sibp3 ) );
#! [ [ |[ 1/2*s12*a9+1/2*s14*a9 ]|, |[ D4*D7*N9_ ]| ],
#!   [ |[ -a7+1 ]|, |[ D1*D4 ]| ],
#!   [ |[ a7-1 ]|, |[ D3*D4 ]| ],
#!   [ |[ a7-1 ]|, |[ D2*D5 ]| ],
#!   [ |[ -a7+1 ]|, |[ D1*D6 ]| ],
#!   [ |[ d-a4-a5-a6-a7-a9+1 ]|, |[ D5*D7 ]| ],
#!   [ |[ -d+a4+a5+a6+a7+a9-1 ]|, |[ D6*D7 ]| ],
#!   [ |[ -2*a7+2 ]|, |[ D4*N8 ]| ],
#!   [ |[ -2*a7+2 ]|, |[ D1*N9 ]| ],
#!   [ |[ -2*d+2*a4+2*a5+2*a6+2*a7+2*a9-2 ]|, |[ D7*N9 ]| ],
#!   [ |[ -4*a7+4 ]|, |[ N8*N9 ]| ],
#!   [ |[ -s12*a4-s14*a5+s12*a6+s14*a9 ]|, |[ D7 ]| ] ]
sibp3 := IBPRelation( CertainColumns( S, [ 3 ] ), LD, [ 1, 1, 1, 1, 1, 1, 1, 1, 1 ] );
#! |[ 1/2*s12*D4*D7*N9_+1/2*s14*D4*D7*N9_+d*D5*D7-d*D6*D7\
#!    -2*d*D7*N9-4*D5*D7+4*D6*D7+8*D7*N9 ]|
ViewList( DecomposeInMonomials( sibp3 ) );
#! [ [ |[ 1/2*s12+1/2*s14 ]|, |[ D4*D7*N9_ ]| ],
#!   [ |[ d-4 ]|, |[ D5*D7 ]| ],
#!   [ |[ -d+4 ]|, |[ D6*D7 ]| ],
#!   [ |[ -2*d+8 ]|, |[ D7*N9 ]| ] ]
gen := GeneratorsOfScalelessSectors( LD );
#! <A 1 x 10 matrix over an external ring>
EntriesOfHomalgMatrix( gen );
#! [ D6*D7*N8*N9, D4*D7*N8*N9, D3*D7*N8*N9, D1*D7*N8*N9, D4*D5*D6*N8*N9,
#!   D3*D5*D6*N8*N9, D2*D3*D6*N8*N9, D1*D4*D5*N8*N9, D1*D2*D4*N8*N9,
#!   D1*D2*D3*N8*N9 ]
#! @EndExample
