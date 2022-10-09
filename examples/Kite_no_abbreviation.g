#! @Chunk Kite_no_abbreviation

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1..2", "p", 4 : abbreviation := false );
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
#! -2*x1,     -2*x3,     -2*x4,     0,         0,         0,
#! -2*x1-2*x4,-2*x3-2*x5,-2*x4-2*x6,0,         0,         0,
#! 0,         0,         0,         -2*x3,     -2*x2,     -2*x5,
#! 0,         0,         0,         -2*x3+2*x4,-2*x2+2*x5,-2*x5+2*x6,
#! -2*x1-2*x3,-2*x2-2*x3,-2*x4-2*x5,-2*x1-2*x3,-2*x2-2*x3,-2*x4-2*x5
Display( e12[2] );
#! -x1,0,          0,  0,          0,
#! 0,  -x1-2*x4-x6,0,  0,          0,
#! 0,  0,          -x2,0,          0,
#! 0,  0,          0,  -x2+2*x5-x6,0,
#! 0,  0,          0,  0,          -x1-x2-2*x3
#! @EndExample
