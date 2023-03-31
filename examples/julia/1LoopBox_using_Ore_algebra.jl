using CapAndHomalg

LoadPackage( "LoopIntegrals" )

LD = LoopDiagram( "l1", "k1..2,k4" )

rel1 = List( ExternalMomenta( LD ), k -> k^2 )

rel2 = GapObj([ (k1+k2+k4)^2 ]);

SetRelationsOfExternalMomenta( LD, Concatenation( rel1, rel2 ) )

s12 = 2*k1*k2

SetAbbreviation( s12, "s12" )

s14 = 2*k1*k4

SetAbbreviation( s14, "s14" )

SetPropagators( LD, -[ l1^2, (l1-k1)^2, (l1-k1-k2)^2, (l1+k4)^2 ] )

SetNumerators( LD, -[ ] )

SetIndependentLorentzInvariants( LD, [ l1^2, l1*k1, l1*k2, l1*k4, s12, s14 ] )

SetExtraLorentzInvariants( LD, [ s12, s14 ] )

R = RingOfLoopDiagram( LD )

Ypol = DoubleShiftAlgebra( R )

R = RingOfLoopDiagram( LD )

ibps = MatrixOfIBPRelations( LD )

r1 = ibps[1,1]

ViewList( DecomposeInMonomials( r1 ) )

Gpol = BasisOfRows( ibps )

NormalForm( "a1*D1_" /  Ypol, Gpol )

gen = GeneratorsOfScalelessSectors( LD )

Display( gen )

E12 = PairOfMatricesOfLoopDiagramInPropagators( LD )

Display( E12[1] )

S = SyzygiesOfColumns( E12 )

Sred = ReducedBasisOfColumnModule( S )

Display( Sred )

Display( EntriesOfHomalgMatrixAsListList( CertainColumns( Sred, Array( 1:3 ) ) ) )

Sibp1 = IBPRelation( CertainColumns( Sred, [ 1 ] ), LD )

ViewList( DecomposeInMonomials( Sibp1 ) )

Sibp2 = IBPRelation( CertainColumns( Sred, [ 2 ] ), LD )

ViewList( DecomposeInMonomials( Sibp2 ) )

Sibp3 = IBPRelation( CertainColumns( Sred, [ 3 ] ), LD )

ViewList( DecomposeInMonomials( Sibp3 ) )

SGpol = BasisOfSpecialIBPRelations( LD )

Gpol == SGpol

Y = RationalDoubleShiftAlgebra( R )

ribps = Y * ibps

G = BasisOfRows( ribps )

NormalForm( "a1*D1_" /  Y, G )

NormalForm( "a2*D2_" /  Y, G )

NormalForm( "a3*D3_" /  Y, G )

NormalForm( "a4*D4_" /  Y, G )

NormalForm( "1" /  Y, G )

NormalForm( "D1" /  Y, G )

NormalForm( "D2" /  Y, G )

NormalForm( "D3" /  Y, G )

NormalForm( "D4" /  Y, G )

NormalFormWrtInitialIntegral( "D1_" / Y, G )

NormalFormWrtInitialIntegral( "D3_" / Y, G )

NormalFormWrtInitialIntegral( "D1*D2" / Y, G )
