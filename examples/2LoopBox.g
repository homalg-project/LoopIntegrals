#! @Chunk 2LoopBox

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1..2", "k1..2,k4", 1 );
#! <A loop diagram with loop momenta [ l1, l2 ] &
#!  external momenta [ k1, k2, k4 ]>
s12 := 2*k1*k2;;
SetAbbreviation( s12, "s12" );
s14 := 2*k1*k4;;
SetAbbreviation( s14, "s14" );
rel1 := List( ExternalMomenta( LD ), k -> k^2 );;
rel2 := [ (k1+k2+k4)^2 ];;
SetRelationsOfMomenta( LD, Concatenation( rel1, rel2 ) );
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
#! [ <A 10 x 7 matrix over an external ring>,
#!   <A 7 x 7 matrix over an external ring> ]
S := SyzygiesOfRows( E12 );
#! <A non-zero 85 x 10 matrix over an external ring>
Display( EntriesOfHomalgMatrix( S[1] ) );
#! [ -D2+D3+2*N8, 0, -2*N8, D1, -D2, 0, -2*D4+D5+D6-2*N9, -D4+D5, -D4, D4+2*N9 ]
Display( EntriesOfHomalgMatrix( S[2] ) );
#! [ D1-D2+2*N8, 0, -D1-2*N8, 0, -D2, 0, -D4+D5-2*N9, D5, 0, D4+2*N9 ]
Display( EntriesOfHomalgMatrix( S[3] ) );
#! [ 0, 0, 0, 0, 0, 0, D5*D7-D6*D7+2*D7*N9-D7*s12, 0, D5*D7, D6*D7-2*D7*N9+D7*s12 ]
#! @EndExample
