#! @Chunk LorentzVector

LoadPackage( "LoopIntegrals" );

#! @Example
LOOP_INTEGRALS.Dimension := 4;
#! 4
LD := LoopDiagram( "l1..2", "p1..2" );
#! <A loop diagram with loop momenta [ l1, l2 ] & external momenta [ p1, p2 ]>
l1;
#! l1
Display( l1 );
#! l1_0,
#! l1_1,
#! l1_2,
#! l1_3
-l1;
#! <A 4-vector>
Display( -l1 );
#! -l1_0,
#! -l1_1,
#! -l1_2,
#! -l1_3
l1^0;
#! 1
l1^2;
#! l1_0^2-l1_1^2-l1_2^2-l1_3^2
l1*l1*p1;
#! <A 4-vector>
Display( l1*l1*p1 );
#! l1_0^2*p1_0-l1_1^2*p1_0-l1_2^2*p1_0-l1_3^2*p1_0,
#! l1_0^2*p1_1-l1_1^2*p1_1-l1_2^2*p1_1-l1_3^2*p1_1,
#! l1_0^2*p1_2-l1_1^2*p1_2-l1_2^2*p1_2-l1_3^2*p1_2,
#! l1_0^2*p1_3-l1_1^2*p1_3-l1_2^2*p1_3-l1_3^2*p1_3
l1*l1*p1*p1 = (l1^2)*(p1^2);
#! true
l1*(l1*p1)*p1 = (l1*(l1*p1))*p1;
#! true
#! @EndExample
