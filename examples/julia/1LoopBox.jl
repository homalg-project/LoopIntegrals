
using HomalgProject

LoadPackage( "LoopIntegrals" )

LD = LoopDiagram( "l1", "k1..2,k4" )

s12 = 2*k1*k2

SetAbbreviation( s12, "s12" )

s14 = 2*k1*k4

SetAbbreviation( s14, "s14" )

rel1 = List( ExternalMomenta( LD ), k -> k^2 )

rel2 = @gap([ (k1+k2+k4)^2 ]);

SetRelationsOfMomenta( LD, Concatenation( rel1, rel2 ) )

SetIndependentLorentzInvariants( LD, [ l1^2, l1*k1, l1*k2, l1*k4, s12, s14 ] )

SetPropagators( LD, -[ l1^2, (l1-k1)^2, (l1-k1-k2)^2, (l1+k4)^2 ] )

SetNumerators( LD, -[ ] )

SetExtraLorentzInvariants( LD, [ s12, s14 ] )

R = RingOfLoopDiagram( LD )

id = HomalgIdentityMatrix( DimensionOfCoefficientsVector( LD ), R )

ibp1 = IBPRelation( id[1], LD )

ViewList( DecomposeInMonomials( ibp1 ) )

HomalgRing( ibp1 )

E12 = PairOfMatricesOfLoopDiagramInPropagators( LD )

Display( E12[1] )

S = SyzygiesOfRows( E12 )

Sred = ReducedBasisOfRowModule( S )

Display( Sred )

Display( EntriesOfHomalgMatrixAsListList( CertainRows( Sred, Array( 1:3 ) ) ) )

Sibp1 = IBPRelation( Sred[1], LD )

ViewList( DecomposeInMonomials( Sibp1 ) )

sibp1 = IBPRelation( Sred[1], LD, [ 1, 1, 1, 1 ] )

ViewList( DecomposeInMonomials( sibp1 ) )

Sibp2 = IBPRelation( Sred[2], LD )

ViewList( DecomposeInMonomials( Sibp2 ) )

sibp2 = IBPRelation( Sred[2], LD, [ 1, 1, 1, 1 ] )

ViewList( DecomposeInMonomials( sibp2 ) )

Sibp3 = IBPRelation( Sred[3], LD )

ViewList( DecomposeInMonomials( Sibp3 ) )

sibp3 = IBPRelation( Sred[3], LD, [ 1, 1, 1, 1 ] )

ViewList( DecomposeInMonomials( sibp2 ) )
