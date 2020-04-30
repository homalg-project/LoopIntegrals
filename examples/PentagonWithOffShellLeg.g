#! @Chunk PentagonWithOffShellLeg

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "k1", "p1, p2, p3, p4", 1 );
#! <A loop diagram with loop momenta [ k1 ] &
#!  external momenta [ p1, p2, p3, p4 ]>
pp5:= 2*p1*p2 + 2*p1*p3 + 2*p2*p3 + 2*p1*p4 + 2*p2*p4 + 2*p3*p4;;
SetAbbreviation( pp5, "pp5" );
s12:= 2*p1*p2;;
SetAbbreviation( s12, "s12" );
s23:= 2*p2*p3;;
SetAbbreviation( s23, "s23" );
s34:= 2*p3*p4;;
SetAbbreviation( s34, "s34" );
s45:= 2*p1*p2 + 2*p1*p3 + 2*p2*p3;;
SetAbbreviation( s45, "s45" );
s51:= 2*p2*p3 + 2*p2*p4 + 2*p3*p4;;
SetAbbreviation( s51, "s51" );
rel := [ p1^2, p2^2, p3^2, p4^2 ];;
SetRelationsOfMomenta( LD, rel );
SetIndependentLorentzInvariants( LD,
        [ k1^2, k1*p1, k1*p2, k1*p3, k1*p4, pp5, s12, s23, s34, s45, s51 ] );
SetPropagators( LD,
        -[ k1^2, (k1 + p1)^2, (k1 + p1 + p2)^2, (k1 + p1 + p2 + p3)^2,
           (k1 + p1 + p2 + p3 + p4)^2 ] );
SetNumerators( LD, -[  ] );
SetExtraLorentzInvariants( LD, [ pp5, s12, s23, s34, s45, s51 ] );
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <An unevaluated non-zero 5 x 5 matrix over an external ring>,
#!   <An unevaluated non-zero 5 x 5 matrix over an external ring> ]
Display( E12[1] );
#! 2*D1,          D1+D2,         D1+D3+s12, D1+D4+s45,     D1+D5+pp5,
#! -D1+D2,        -D1+D2,        -D1+D2-s12,-D1+D2+s23-s45,-D1+D2-pp5+s51,
#! -D2+D3+s12,    -D2+D3,        -D2+D3,    -D2+D3-s23,    -D2+D3+s34-s51,
#! -D3+D4-s12+s45,-D3+D4+s23,    -D3+D4,    -D3+D4,        -D3+D4-s34,
#! -D4+D5+pp5-s45,-D4+D5-s23+s51,-D4+D5+s34,-D4+D5,        -D4+D5
S := SyzygiesOfRows( E12 );
#! <A non-zero 56 x 5 matrix over an external ring>
#! @EndExample
