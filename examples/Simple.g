#! @Chunk Simple

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "k", "": parameters := "m" );
#! <A loop diagram with loop momenta [ k ] & external momenta [ ]>
SetRelationsOfMomenta( LD, [ ] );
SetIndependentLorentzInvariants( LD, [ k^2 ] );
SetPropagators( LD, [ k^2 - m^2 ] );
SetNumerators( LD, [ ] );
SetExtraLorentzInvariants( LD, [ ] );
R := RingOfLoopDiagram( LD );
#! Q[m,D][D1]
ibps := MatrixOfIBPRelations( LD );
#! <A 1 x 1 matrix over a residue class ring>
ibp := ibps[1,1];
#! |[ -2*m^2*a1*D1_+D-2*a1 ]|
ViewList( DecomposeInMonomials( ibp ) );
#! [ [ |[ -2*m^2*a1 ]|, |[ D1_ ]| ],
#!   [ |[ D-2*a1 ]|, |[ 1 ]| ] ]
Y := HomalgRing( ibp );
#! Q[m,D][a1]<D1,D1_>/( D1*D1_-1 )
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <A 1 x 1 matrix over an external ring>,
#!   <A 1 x 1 matrix over an external ring> ]
Display( E12[1] );
#! 2*m^2+2*D1
Display( E12[2] );
#! D1
S := SyzygiesOfRows( E12 );
#! <A non-zero 1 x 1 matrix over an external ring>
Display( S );
#! D1
Sibp := IBPRelation( S[1], LD );
#! |[ -2*m^2*a1+2*m^2+D*D1-2*a1*D1+2*D1 ]|
ViewList( DecomposeInMonomials( Sibp ) );
#! [ [ |[ D-2*a1+2 ]|, |[ D1 ]| ],
#!   [ |[ -2*m^2*a1+2*m^2 ]|, |[ 1 ]| ] ]
Sibps := MatrixOfSpecialIBPRelations( LD );
#! <A 1 x 1 matrix over a residue class ring>
x := RightDivide( Sibps, ibps );
#! <An unevaluated 1 x 1 matrix over a residue class ring>
x * ibps = Sibps;
#! true
Display( x );
#! D1
#! 
#! modulo [ D1*D1_-1 ]
y := RightDivide( ibps, Sibps );
#! <An unevaluated 1 x 1 matrix over a residue class ring>
y * Sibps = ibps;
#! true
Display( y );
#! D1_
#! 
#! modulo [ D1*D1_-1 ]
bas := BasisOfIBPRelations( LD );
#! <A non-zero 1 x 1 matrix over a residue class ring>
Sbas := BasisOfSpecialIBPRelations( LD );
#! <A non-zero 1 x 1 matrix over a residue class ring>
Sbas = bas;
#! true
r := "2 * m^2 * ( a1 + 4 - 2 ) * D1_^(4-1)" / Y;
#! |[ 2*m^2*a1*D1_^3+4*m^2*D1_^3 ]|
d := DecideZero( r, bas );
#! |[ D*D1_^2-2*a1*D1_^2-4*D1_^2 ]|
DecideZero( r - d, bas );
#! |[ 0 ]|
s := "D1_^(4 - 2) * ( D - 2 * a1 )" / Y;
#! |[ D*D1_^2-2*a1*D1_^2-4*D1_^2 ]|
DecideZero( r - s, bas );
#! |[ 0 ]|
#! @EndExample
