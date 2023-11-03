#! @Chunk 1LoopBubble

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1", "k1" );
#! <A loop diagram with loop momenta [ l1 ] & external momenta [ k1 ]>
s := k1*k1;
#! k1^2
SetAbbreviation( s, "s" );
rel := [ ];;
SetRelationsOfExternalMomenta( LD, rel );
SetIndependentLorentzInvariants( LD, [ ] );
SetPropagators( LD, -[ l1^2, ( l1 + k1 )^2 ] );
SetNumerators( LD, -[ ] );
SetExtraLorentzInvariants( LD, [ s ] );
R := RingOfLoopDiagram( LD );
#! Q[d,s][D1,D2]
ibps := MatrixOfIBPRelations( LD );
#! <A 2 x 1 matrix over a residue class ring>
ibp1 := ibps[1,1];
#! |[ -s*a2*D2_-a2*D1*D2_+d-2*a1-a2 ]|
ViewList( DecomposeInMonomials( ibp1 ) );
#! [ [ |[ -a2 ]|, |[ D1*D2_ ]| ],
#!   [ |[ -s*a2 ]|, |[ D2_ ]| ],
#!   [ |[ d-2*a1-a2 ]|, |[ 1 ]| ] ]
Ypol := HomalgRing( ibp1 );
#! Q[d,s][a1,a2]<D1,D1_,D2,D2_>/( D2*D2_-1, D1*D1_-1 )
bas := BasisOfRows( ibps );
#! <A non-zero 4 x 1 matrix over a residue class ring>
gen := GeneratorsOfScalelessSectors( LD );
#! <A 1 x 2 matrix over an external ring>
Display( gen );
#! D2,D1
gen2 := GeneratorsOfScalelessSectors( LD, [ 2, 2 ] );
#! <An unevaluated 1 x 2 matrix over an external ring>
Display( gen2 );
#! D1*D2^2,D1^2*D2
prel2 := ColumnReversedMatrixOfCoefficientsOfParametricIBPs( LD, 2 );
#! [ <A non-zero 5 x 6 matrix over an external ring>,
#!   [ D1_^2, D1_*D2_, D2_^2, D1_, D2_, 1 ] ]
EntriesOfHomalgMatrix( prel2[1][4] );
#! [ 0, 0, 0,
#!   (d*s*a1-2*s*a1^2-2*s*a1), 0,
#!   (-d^2+3*d*a1+3*d*a2+d-2*a1^2-4*a1*a2-2*a1-2*a2^2-2*a2) ]
EntriesOfHomalgMatrix( prel2[1][5] );
#! [ 0, 0, 0,
#!   0, (d*s*a2-2*s*a2^2-2*s*a2),
#!   (-d^2+3*d*a1+3*d*a2+d-2*a1^2-4*a1*a2-2*a1-2*a2^2-2*a2) ]
#! @EndExample
