#! @Chunk 1LoopTadpole

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "k", "" : masses := "m" );
#! <A loop diagram with loop momenta [ k ] & external momenta [ ] &
#!  masses [ m ]>
SetRelationsOfExternalMomenta( LD, [ ] );
SetIndependentLorentzInvariants( LD, [ k^2 ] );
SetPropagators( LD, [ k^2 - m^2 ] );
SetNumerators( LD, [ ] );
SetExtraLorentzInvariants( LD, [ ] );
R := RingOfLoopDiagram( LD );
#! Q[m,d][D1]
ibps := MatrixOfIBPRelations( LD );
#! <A 1 x 1 matrix over a residue class ring>
ibp := ibps[1,1];
#! |[ -2*m^2*a1*D1_+d-2*a1 ]|
ViewList( DecomposeInMonomials( ibp ) );
#! [ [ |[ -2*m^2*a1 ]|, |[ D1_ ]| ],
#!   [ |[ d-2*a1 ]|, |[ 1 ]| ] ]
Ypol := HomalgRing( ibp );
#! Q[m,d][a1]<D1,D1_>/( D1*D1_-1 )
ibpws := MatrixOfIBPRelationsInWeylAlgebra( LD );
#! <A 1 x 1 matrix over an external ring>
ibpw := MatElm( ibpws, 1, 1 );
#! -2*m^2*A1+d-2*D1*A1-2
W := HomalgRing( ibpw );
#! Q[m,d][D1]<A1>
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <A 1 x 1 matrix over an external ring>,
#!   <A 1 x 1 matrix over an external ring> ]
Display( E12[1] );
#! 2*m^2+2*D1
Display( E12[2] );
#! D1
S := SyzygiesOfColumns( E12 );
#! <A non-zero 1 x 1 matrix over an external ring>
Display( S );
#! D1
Sibp := IBPRelation( CertainColumns( S, [ 1 ] ), LD );
#! |[ -2*m^2*a1+2*m^2+d*D1-2*a1*D1+2*D1 ]|
ViewList( DecomposeInMonomials( Sibp ) );
#! [ [ |[ d-2*a1+2 ]|, |[ D1 ]| ],
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
ExportVariables( Ypol );
#! [ |[ m ]|, |[ d ]|, |[ a1 ]|, |[ D1 ]|, |[ D1_ ]| ]
r := 2 * m^2 * ( a1 + 4 - 2 ) * D1_^(4-1);
#! |[ 2*m^2*a1*D1_^3+4*m^2*D1_^3 ]|
n := DecideZero( r, bas );
#! |[ d*D1_^2-2*a1*D1_^2-4*D1_^2 ]|
DecideZero( n - r, bas );
#! |[ 0 ]|
s := D1_^(4 - 2) * ( d - 2 * a1 );
#! |[ d*D1_^2-2*a1*D1_^2-4*D1_^2 ]|
s = n;
#! true
t := ( d - 2 * a1 - 2 * ( 4 - 2 ) ) * D1_^(4-2);
#! |[ d*D1_^2-2*a1*D1_^2-4*D1_^2 ]|
t = n;
#! true
w := ( d - 2 * 1 - 2 * ( 4 - 2 ) ) * D1_^(4-2);
#! |[ d*D1_^2-6*D1_^2 ]|
DecideZero( w, bas );
#! |[ d*D1_^2-6*D1_^2 ]|
#! @EndExample
