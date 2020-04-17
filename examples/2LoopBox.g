#! @Chunk 2LoopBox

LoadPackage( "LoopIntegrals" );

#! @Example
LOOP_INTEGRALS.Dimension := 2;
#! 2
LD := LoopDiagram( "l1..2", "k1..2,k4" );
#! <A loop diagram with loop momenta [ l1, l2 ] &
#!  external momenta [ k1, k2, k4 ]>
SetPropagators( LD,
        [ l1^2, (l1-k1)^2, (l1-k1-k2)^2,
          l2^2, (l2-k4)^2, (l2+k1+k2)^2,
          (l1+l2)^2 ] );
SetNumerators( LD, [ l1*k4, l2*k1 ] );
rel1 := List( ExternalMomenta( LD ), k -> k^2 );;
rel2 := [ (k1+k2+k4)^2 ];;
SetRelationsOfMomenta( LD, Concatenation( rel1, rel2 ) );
SetIndependetLorentzInvariants( LD,
        [ l1^2, l1*l2, l2^2, l1*k1, l1*k2, l1*k4,
          l2*k1, l2*k2, l2*k4, 2*k1*k2, 2*k1*k4 ] );
SetExtraLorentzInvariants( LD, [ 2*k1*k2, 2*k1*k4 ] );
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <A non-zero 10 x 7 matrix over an external ring>,
#!   <A non-zero 7 x 7 matrix over an external ring> ]
S := SyzygiesOfRows( E12[1], E12[2] );
#! <A non-zero 85 x 10 matrix over an external ring>
#! @EndExample
