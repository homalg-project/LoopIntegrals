#
# LoopIntegrals: Loop integrals
#
# Implementations
#

InstallValue( LOOP_INTEGRALS,
        rec(
            PropagatorSymbol := "D",
            NumeratorSymbol := "N",
            LorentzSymbol := "x"
           ) );

##
InstallMethod( LorentzVector,
        [ IsHomalgMatrix ],
        
  function( M )
    local vector;
    
    vector := rec( );
    
    ObjectifyWithAttributes(
            vector, TheTypeLorentzVector,
            Dimension, NrRows( M ),
            UnderlyingMatrix, M );
    
    return vector;
    
end );
    
##
InstallMethod( LorentzVector,
        [ IsString and IsStringRep, IsInt ],
        
  function( str, dim )
    local vector;
    
    vector := rec( symbols := List( [ 0 .. dim - 1 ], i -> Concatenation( str, "_", String( i ) ) ) );
    
    ObjectifyWithAttributes(
            vector, TheTypeLorentzVector,
            Dimension, dim,
            Name, ShallowCopy( str ) );
    
    if ValueOption( "IsLoopMomentum" ) = true then
        SetIsLoopMomentum( vector, true );
    fi;
    
    if ValueOption( "IsExternalMomentum" ) = true then
        SetIsExternalMomentum( vector, true );
    fi;
    
    return vector;
    
end );
    
##
InstallMethod( LoopDiagram,
        [ IsList, IsList, IsInt ],
        
  function( L, K, dim )
    local LD, LorentzVectors, C, momenta, R;
    
    LD := rec( );
    
    if IsStringRep( L ) then
        L := [ L ];
    fi;
    
    L := Concatenation( List( List( L, l -> SplitString( l, "," ) ), ParseListOfIndeterminates ) );
    
    L := List( L, l -> LorentzVector( l, dim : IsLoopMomentum := true ) );
    
    Perform( L,
      function( l )
        local name;
        name := Name( l );
        if IsBoundGlobal( name ) then
            MakeReadWriteGlobal( name );
            UnbindGlobal( name );
        fi;
        BindGlobal( name, l );
    end );
    
    if IsStringRep( K ) then
        K := [ K ];
    fi;
    
    K := Concatenation( List( List( K, k -> SplitString( k, "," ) ), ParseListOfIndeterminates ) );
    
    K := List( K, k -> LorentzVector( k, dim : IsExternalMomentum := true ) );
    
    Perform( K,
      function( k )
        local name;
        name := Name( k );
        if IsBoundGlobal( name ) then
            MakeReadWriteGlobal( name );
            UnbindGlobal( name );
        fi;
        BindGlobal( name, k );
    end );
    
    if ValueOption( "C" ) = fail then
        C := HomalgFieldOfRationalsInSingular( );
    fi;
    
    momenta := Concatenation( List( [ L, K ], M -> Concatenation( List( M, m -> m!.symbols ) ) ) );
    
    R := C * momenta;
    
    ObjectifyWithAttributes(
            LD, TheTypeLoopDiagram,
            LoopMomenta, L,
            ExternalMomenta, K,
            UnderlyingRing, R );
    
    Perform( L, function( l ) SetLoopDiagram( l, LD ); end );
    Perform( K, function( k ) SetLoopDiagram( k, LD ); end );
    
    return LD;
    
end );

##
InstallMethod( Components,
        [ IsLorentzVector and HasLoopDiagram ],
        
  function( x )
    local R;
    
    R := UnderlyingRing( LoopDiagram( x ) );
    
    return List( x!.symbols, c -> c / R );
    
end );
    
##
InstallMethod( UnderlyingMatrix,
        [ IsLorentzVector ],
        
  function( x )
    
    return HomalgMatrix( Components( x ), Dimension( x ), 1, HomalgRing( x[1] ) );
    
end );
    
##
InstallMethod( Components,
        [ IsLorentzVector and HasUnderlyingMatrix ],
        
  function( x )
    
    return EntriesOfHomalgMatrix( UnderlyingMatrix( x ) );
    
end );
    
##
InstallMethod( \[\],
        [ IsLorentzVector, IsInt ],
        
  function( x, i )
    
    return Components( x )[i + 1];
    
end );
    
##
InstallMethod( \+,
        [ IsLorentzVector and HasLoopDiagram, IsLorentzVector ],
        
  function( x, y )
    local xy;
    
    xy := LorentzVector( UnderlyingMatrix( x ) + UnderlyingMatrix( y ) );
    
    SetLoopDiagram( xy, LoopDiagram( x ) );
    
    return xy;
    
end );

##
InstallMethod( \-,
        [ IsLorentzVector and HasLoopDiagram, IsLorentzVector ],
        
  function( x, y )
    local xy;
    
    xy := LorentzVector( UnderlyingMatrix( x ) - UnderlyingMatrix( y ) );
    
    SetLoopDiagram( xy, LoopDiagram( x ) );
    
    return xy;
    
end );

##
InstallMethod( AdditiveInverse,
        [ IsLorentzVector and HasLoopDiagram ],
        
  function( x )
    local mx;
    
    mx := LorentzVector( -UnderlyingMatrix( x ) );
    
    SetLoopDiagram( mx, LoopDiagram( x ) );
    
    return mx;
    
end );

##
InstallMethod( \*,
        [ IsLorentzVector, IsLorentzVector ],
        
  function( x, y )
    local cx, cy, xy;
    
    cx := Components( x );
    cy := Components( y );
    
    xy := cx[1] *cy[1] - Sum( [ 2 .. Dimension( x ) ], i -> cx[i] * cy[i] );
    
    return xy;
    
end );

##
InstallMethod( \*,
        [ IsRingElement, IsLorentzVector and HasLoopDiagram ],
        
  function( r, x )
    local rx;
    
    rx :=  LorentzVector( r * UnderlyingMatrix( x ) );

    SetLoopDiagram( rx, LoopDiagram( x ) );

    return rx;
    
end );

##
InstallMethod( \*,
        [ IsLorentzVector and HasLoopDiagram, IsRingElement ],
        
  function( x, r )
    local rx;
    
    rx :=  LorentzVector( r * UnderlyingMatrix( x ) );
    
    SetLoopDiagram( rx, LoopDiagram( x ) );
    
    return rx;
    
end );

##
InstallMethod( POW,
        [ IsLorentzVector and HasLoopDiagram, IsInt ],
        
  function( x, power )
    
    if power = 0 then
        return One( UnderlyingRing( LoopDiagram( x ) ) );
    fi;
    
    return Product( ListWithIdenticalEntries( power, x ) );
    
end );

##
InstallMethod( RelationsMatrixOfMomenta,
        [ IsLoopDiagram and HasRelationsOfMomenta ],
        
  function( LD )
    local rel;
    
    rel := RelationsOfMomenta( LD );
    
    return HomalgMatrix( rel, Length( rel ), 1, UnderlyingRing( LD ) );
    
end );

##
InstallMethod( ReductionMatrixOfIndependetLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators ],
        
  function( LD )
    local symbol, I, M, invariants, R, red;
    
    symbol := ValueOption( "symbol" );
    
    if symbol = fail then
        symbol := LOOP_INTEGRALS.LorentzSymbol;
    fi;
    
    I := IndependetLorentzInvariants( LD );
    
    M := Length( I );
    
    invariants := List( [ 1 .. M ], i -> Concatenation( symbol, String( i ) ) );
    
    R := UnderlyingRing( LD );
    
    R := CoefficientsRing( R ) * invariants * List( Indeterminates( R ), String );
    
    R := PolynomialRingWithProductOrdering( R );
    
    invariants := List( invariants, i -> i / R );
    
    invariants := ListN( I, invariants, {a,b} -> a / R - b );
    
    invariants := HomalgMatrix( invariants, M, 1, R );
    
    red := UnionOfRows( invariants, R * RelationsMatrixOfMomenta( LD ) );
    
    return red;
    
end );

##
InstallMethod( ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    local symbolD, D, N, propagators, symbolN, Z, A, numerators,
          symbolK, K, M, invariants, indets, R, red;
    
    symbolD := ValueOption( "symbolD" );
    
    if symbolD = fail then
        symbolD := LOOP_INTEGRALS.PropagatorSymbol;
    fi;
    
    D := Propagators( LD );
    
    N := Length( D );
    
    propagators := List( [ 1 .. N ], i -> Concatenation( symbolD, String( i ) ) );
    
    symbolN := ValueOption( "symbolN" );
    
    if symbolN = fail then
        symbolN := LOOP_INTEGRALS.NumeratorSymbol;
    fi;
    
    Z := Numerators( LD );
    
    A := Length( Z );
    
    numerators := List( [ N + 1 .. N + A ], i -> Concatenation( symbolN, String( i ) ) );
    
    symbolK := ValueOption( "symbolK" );
    
    if symbolK = fail then
        symbolK := LOOP_INTEGRALS.LorentzSymbol;
    fi;
    
    K := ExtraLorentzInvariants( LD );
    
    M := Length( K );
    
    invariants := List( [ N + A + 1 .. N + A + M ], i -> Concatenation( symbolK, String( i ) ) );
    
    indets := Concatenation( propagators, numerators, invariants );
    
    R := UnderlyingRing( LD );
    
    R := CoefficientsRing( R ) * indets * List( Indeterminates( R ), String );
    
    R := PolynomialRingWithProductOrdering( R );
    
    indets := List( indets, p -> p / R );
    
    indets := ListN( Concatenation( D, Z, K ), indets, {a,b} -> a / R - b );
    
    indets := HomalgMatrix( indets, N + A + M, 1, R );
    
    red := UnionOfRows( indets, R * RelationsMatrixOfMomenta( LD ) );
    
    return red;
    
end );

##
InstallMethod( PairOfOriginalMatricesOfLoopDiagram,
        [ IsLoopDiagram and HasRelationsOfMomenta ],
        
  function( LD )
    local l, R, propagators, E1, E2;
    
    l := UnionOfRows( List( LoopMomenta( LD ), UnderlyingMatrix ) );
    
    R := UnderlyingRing( LD );
    
    propagators := Propagators( LD );
    
    E1 := UnionOfColumns( List( propagators, p -> Diff( l, HomalgMatrix( p, 1, 1, R ) ) ) );
    
    E2 := HomalgDiagonalMatrix( propagators );
    
    return [ E1, E2 ];
    
end );

##
InstallMethod( PairOfMatricesOfLoopDiagram,
        [ IsLoopDiagram and HasRelationsOfMomenta ],
        
  function( LD )
    local pair, l, k, vectors, rel;
    
    pair := PairOfOriginalMatricesOfLoopDiagram( LD );
    
    l := List( LoopMomenta( LD ), UnderlyingMatrix );
    k := List( ExternalMomenta( LD ), UnderlyingMatrix );
    
    vectors := Involution( UnionOfColumns( Concatenation( l, k ) ) );
    
    vectors := DiagMat( ListWithIdenticalEntries( Length( l ), vectors ) );
    
    pair := [ vectors * pair[1], pair[2] ];
    
    rel := RelationsMatrixOfMomenta( LD );
    
    rel := DiagMat( ListWithIdenticalEntries( NrColumns( pair[1] ), rel ) );
    
    pair[1] := DecideZeroRows( pair[1], BasisOfRows( rel ) );
    
    return pair;
    
end );

##
InstallMethod( PairOfMatricesOfLoopDiagramInLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfMomenta ],
        
  function( LD )
    local pair, col, red, R;
    
    pair := PairOfMatricesOfLoopDiagram( LD );
    
    col := NrColumns( pair[1] );
    
    red := ReductionMatrixOfIndependetLorentzInvariants( LD );
    
    red := BasisOfRows( red );
    
    red := DiagMat( ListWithIdenticalEntries( col, red ) );
    
    R := HomalgRing( red );
    
    return List( pair, mat -> DecideZeroRows( R * mat, red ) );
    
end );

##
InstallMethod( PairOfMatricesOfLoopDiagramInPropagators,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasExtraLorentzInvariants ],
        
  function( LD )
    local pair, col, red, R;
    
    pair := PairOfMatricesOfLoopDiagram( LD );
    
    col := NrColumns( pair[1] );
    
    red := ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants( LD );
    
    red := BasisOfRows( red );
    
    red := DiagMat( ListWithIdenticalEntries( col, red ) );
    
    R := HomalgRing( red );
    
    return List( pair, mat -> DecideZeroRows( R * mat, red ) );
    
end );

##
InstallMethod( ViewObj,
        [ IsLorentzVector ],
        
  function( x )
    
    Print( "<A ", Dimension( x ), "-vector>" );
    
end );

##
InstallMethod( Display,
        [ IsLorentzVector ],
        
  function( x )
    
    Display( UnderlyingMatrix( x ) );
    
end );

##
InstallMethod( ViewObj,
        [ IsLoopDiagram ],
        
  function( LD )
    
    Print( "<A loop diagram with loop momenta ", LoopMomenta( LD ),
           " & external momenta ", ExternalMomenta( LD ), ">" );
    
end );
