
using HomalgProject

LoadPackage( "LoopIntegrals" )

LD = LoopDiagram( "k", "", masses = "m" )

SetRelationsOfExternalMomenta( LD, [ ] )

SetIndependentLorentzInvariants( LD, [ k^2 ] )

SetPropagators( LD, [ k^2 - m^2 ] )

SetNumerators( LD, [ ] )

SetExtraLorentzInvariants( LD, [ ] )

R = RingOfLoopDiagram( LD )

ibps = MatrixOfIBPRelations( LD )

ibp = MatElm( ibps, 1, 1 )

ViewList( DecomposeInMonomials( ibp ) )

Y = HomalgRing( ibp )

ibpws = MatrixOfIBPRelationsInWeylAlgebra( LD )

ibpw = MatElm( ibpws, 1, 1 )

W = HomalgRing( ibpw )

E12 = PairOfMatricesOfLoopDiagramInPropagators( LD )

Display( E12[1] )

Display( E12[2] )

S = SyzygiesOfColumns( E12 )

Display( S )

Sibp = IBPRelation( CertainColumns( S, [ 1 ] ), LD )

ViewList( DecomposeInMonomials( Sibp ) )

Sibps = MatrixOfSpecialIBPRelations( LD )

x = RightDivide( Sibps, ibps )

x * ibps == Sibps

Display( x )

y = RightDivide( ibps, Sibps )

y * Sibps == ibps

Display( y )

bas = BasisOfIBPRelations( LD )

Sbas = BasisOfSpecialIBPRelations( LD )

Sbas == bas

ExportVariablesToJulia( Y )

D1_

r = 2 * m^2 * ( a1 + 4 - 2 ) * D1_^(4-1)

d = DecideZero( r, bas )

DecideZero( d - r, bas )

s = D1_^(4 - 2) * ( d - 2 * a1 )

s == d

t = ( d - 2 * a1 - 2 * ( 4 - 2 ) ) * D1_^(4-2)

t == d

w = ( d - 2 * 1 - 2 * ( 4 - 2 ) ) * D1_^(4-2)

DecideZero( w - r, bas )


