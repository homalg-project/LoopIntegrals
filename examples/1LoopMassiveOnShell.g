#! @Chunk 1LoopMassiveOnShell

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "k", "q" : masses := "m" );
#! <A loop diagram with loop momenta [ k ] & external momenta [ q ] &
#!  masses [ m ]>
rel := [ q^2 - m^2 ];
#! [ q^2+(-m^2) ]
SetRelationsOfExternalMomenta( LD, rel );
SetIndependentLorentzInvariants( LD, [  ] );
SetPropagators( LD, -[ k^2, k^2 + 2*q*k ] );
SetNumerators( LD, -[ ] );
SetExtraLorentzInvariants( LD, [ ] );
R := RingOfLoopDiagram( LD );
#! Q[m,d][D1,D2]
ibps := MatrixOfIBPRelations( LD );
#! <A 2 x 1 matrix over a residue class ring>
Display( ibps );
#! -a2*D1*D2_+d-2*a1-a2,
#! 2*m^2*a2*D2_-a1*D1_*D2+a2*D1*D2_+a1-a2
#! 
#! modulo [ D2*D2_-1, D1*D1_-1 ]
ibp1 := ibps[1,1];
#! |[ -a2*D1*D2_+d-2*a1-a2 ]|
ViewList( DecomposeInMonomials( ibp1 ) );
#! [ [ |[ -a2 ]|, |[ D1*D2_ ]| ],
#!   [ |[ d-2*a1-a2 ]|, |[ 1 ]| ] ]
ibp2 := ibps[2,1];
#! |[ 2*m^2*a2*D2_-a1*D1_*D2+a2*D1*D2_+a1-a2 ]|
ViewList( DecomposeInMonomials( ibp2 ) );
#! [ [ |[ -a1 ]|, |[ D1_*D2 ]| ],
#!   [ |[ a2 ]|, |[ D1*D2_ ]| ],
#!   [ |[ 2*m^2*a2 ]|, |[ D2_ ]| ],
#!   [ |[ a1-a2 ]|, |[ 1 ]| ] ]
ibp2_red := DecideZero( ibp2, ibps{[ 1 ]} );
#! |[ 2*m^2*a2*D2_-a1*D1_*D2+d-a1-2*a2 ]|
ViewList( DecomposeInMonomials( ibp2_red ) );
#! [ [ |[ -a1 ]|, |[ D1_*D2 ]| ],
#!   [ |[ 2*m^2*a2 ]|, |[ D2_ ]| ],
#!   [ |[ d-a1-2*a2 ]|, |[ 1 ]| ] ]
Ypol := HomalgRing( ibp1 );
#! Q[m,d][a1,a2]<D1,D1_,D2,D2_>/( D2*D2_-1, D1*D1_-1 )
gen := GeneratorsOfScalelessSectors( LD );
#! <A 1 x 1 matrix over an external ring>
Display( gen );
#! D2
#! @EndExample
