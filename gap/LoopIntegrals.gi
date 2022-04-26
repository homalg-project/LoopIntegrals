# SPDX-License-Identifier: GPL-2.0-or-later
# LoopIntegrals: Compute master integrals using commutative and noncommutative methods from computational algebraic geometry
#
# Implementations
#

InstallValue( LOOP_INTEGRALS,
        rec(
            ConstructorOfDefaultField := HomalgFieldOfRationalsInSingular,
            LorentzSymbol := "x",
            DualSymbol := "z",
            Dimension := 1,
            DimensionSymbol := "D",
            PropagatorSymbol := "D",
            NumeratorSymbol := "N",
            ExponentSymbol := "a",
            WeylSymbol := "A",
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

    if dim = 1 then
        vector := rec( symbols := [ str ] );
    else
        vector := rec( symbols := List( [ 0 .. dim - 1 ], i -> Concatenation( str, "_", String( i ) ) ) );
    fi;
    
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
    local LD, LorentzVectors, masses, C, momenta, R, symbol, symbolD, symbolN;
    
    LD := rec( );
    
    if IsStringRep( L ) then
        L := ShallowCopy( L );
        RemoveCharacters( L, " \n\t\r" );
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
        K := ShallowCopy( K );
        RemoveCharacters( K, " \n\t\r" );
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
    
    masses := ValueOption( "masses" );
    
    if IsStringRep( masses ) then
        masses := SplitString( masses, "," );
    elif not IsList( masses ) then
        masses := [ ];
    fi;
    
    if IsIdenticalObj( ValueOption( "C" ), fail ) then
        if masses = "" then
            C := LOOP_INTEGRALS.ConstructorOfDefaultField( );
        else
            C := LOOP_INTEGRALS.ConstructorOfDefaultField( JoinStringsWithSeparator( masses ) );
        fi;
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
    
    LD!.DimensionSymbol := LOOP_INTEGRALS.DimensionSymbol;
    
    Perform( masses,
      function( m )
        if IsBoundGlobal( m ) then
            MakeReadWriteGlobal( m );
            UnbindGlobal( m );
        fi;
        BindGlobal( m, m / R );
        SetString( ValueGlobal( m ), m );
    end );
    
    LD!.masses := List( masses, ValueGlobal );
    
    symbol := ValueOption( "symbol" );
    
    if IsStringRep( symbol ) then
        LD!.LorentzSymbol := symbol;
    else
        LD!.LorentzSymbol := LOOP_INTEGRALS.LorentzSymbol;
    fi;
    
    symbolD := ValueOption( "symbolD" );
    
    if IsStringRep( symbolD ) then
        LD!.PropagatorSymbol := symbolD;
    else
        LD!.PropagatorSymbol := LOOP_INTEGRALS.PropagatorSymbol;
    fi;
    
    symbolN := ValueOption( "symbolN" );
    
    if IsStringRep( symbolN ) then
        LD!.NumeratorSymbol := symbolN;
    else
        LD!.NumeratorSymbol := LOOP_INTEGRALS.NumeratorSymbol;
    fi;
    
    if IsIdenticalObj( ValueOption( "abbreviation" ), false ) then
        LD!.abbreviation := false;
    else
        LD!.abbreviation := true;
    fi;
    
    ## do not treat the extra Lorentz invariants as rational parameters
    ## as this slows down the syzygies computations in Singular significantly
    if IsIdenticalObj( ValueOption( "rational" ), true ) then
        LD!.rational := true;
    else
        LD!.rational := false;
    fi;
    
    return LD;
    
end );

##
InstallMethod( LoopDiagram,
        [ IsList, IsList ],
        
  function( L, K )
    
    return LoopDiagram( L, K, LOOP_INTEGRALS.Dimension );
    
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
    
    return HomalgMatrix( Components( x ), Dimension( x ), 1, HomalgRing( x[0] ) );
    
end );
    
##
InstallMethod( Components,
        [ IsLorentzVector and HasUnderlyingMatrix ],
        
  function( x )
    
    return EntriesOfHomalgMatrix( UnderlyingMatrix( x ) );
    
end );
    
## L * ( L + E )
InstallMethod( DimensionOfCoefficientsVector,
        [ IsLoopDiagram and HasLoopMomenta and HasExternalMomenta ],
        
  function( LD )
    local L;
    
    L := Length( LoopMomenta( LD ) );
    
    return L * ( L + Length( ExternalMomenta( LD ) ) );
    
end );

##
InstallMethod( RelationsMatrixOfExternalMomenta,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta ],
        
  function( LD )
    local rel;
    
    rel := RelationsOfExternalMomenta( LD );
    
    return HomalgMatrix( rel, Length( rel ), 1, UnderlyingRing( LD ) );
    
end );

##
InstallMethod( OriginalJacobianOfPropagators,
        [ IsLoopDiagram and HasPropagators ],
        
  function( LD )
    local l, R;
    
    l := UnionOfRows( List( LoopMomenta( LD ), UnderlyingMatrix ) );
    
    R := UnderlyingRing( LD );
    
    return UnionOfColumns( List( Propagators( LD ), p -> Diff( l, HomalgMatrix( p, 1, 1, R ) ) ) );
    
end );

##
InstallMethod( OriginalJacobianOfNumerators,
        [ IsLoopDiagram and HasNumerators ],
        
  function( LD )
    local l, R, numerators;
    
    l := UnionOfRows( List( LoopMomenta( LD ), UnderlyingMatrix ) );
    
    R := UnderlyingRing( LD );
    
    numerators := Numerators( LD );
    
    if IsEmpty( numerators ) then
        return HomalgZeroMatrix( NrRows( l ), 0, R );
    fi;
    
    return UnionOfColumns( List( numerators, p -> Diff( l, HomalgMatrix( p, 1, 1, R ) ) ) );
    
end );

##
InstallMethod( OriginalTaylorOfPropagators,
        [ IsLoopDiagram and HasPropagators ],
        
  function( LD )
    local l, lt, R, propagators, jac, symbols, sum, hess, taylor;
    
    l := UnionOfRows( List( LoopMomenta( LD ), UnderlyingMatrix ) );
    lt := TransposedMatrix( l );
    
    R := UnderlyingRing( LD );
    
    propagators := Propagators( LD );
    
    jac := List( propagators, p -> Diff( l, HomalgMatrix( p, 1, 1, R ) ) );
    
    symbols := List( [ 1 .. Length( propagators ) ], i -> Concatenation( LOOP_INTEGRALS.DualSymbol, String( i ) ) );
    
    R := R * symbols;
    
    symbols := List( symbols, a -> a / R );
    
    sum := Sum( ListN( symbols, propagators, {x,p} -> x * ( p / R ) ) );
    
    hess := -Sum( ListN( symbols, jac, {x,p} -> 1/2 * x * ( R * Diff( lt, p ) ) ) );
    
    taylor := [ hess ];
    
    hess := -( ( R * lt ) * hess * ( R * l ) )[1,1];
    
    sum := sum - hess;
    
    jac := 1/2 * Diff( ( R * l ), HomalgMatrix( [ sum ], 1, 1, R ) );
    
    Add( taylor, jac );
    
    jac := ( ( R * lt ) * jac )[1,1];
    
    sum := sum - 2 * jac;
    
    Add( taylor, sum );
    
    R := HomalgRing( ReductionMatrixOfExtraLorentzInvariants( LD ) );
    
    return [ R * taylor[1], R * taylor[2], taylor[3] / R ];
    
end );

##
InstallMethod( MatrixOfMomenta,
        [ IsLoopDiagram ],
        
  function( LD )
    local l, k, vectors;
    
    l := List( LoopMomenta( LD ), UnderlyingMatrix );
    k := List( ExternalMomenta( LD ), UnderlyingMatrix );
    
    vectors := Involution( UnionOfColumns( Concatenation( l, k ) ) );
    
    return DiagMat( ListWithIdenticalEntries( Length( l ), vectors ) );
    
end );

##
InstallMethod( JacobianOfPropagators,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators ],
        
  function( LD )
    local jac, rel;
    
    jac := OriginalJacobianOfPropagators( LD );
    
    jac := MatrixOfMomenta( LD ) * jac;
    
    rel := RelationsMatrixOfExternalMomenta( LD );
    
    jac := List( [ 1 .. NrColumns( jac ) ],
                   j -> DecideZeroRows( CertainColumns( jac, [ j ] ), rel ) );
    
    return UnionOfColumns( jac );
    
end );

##
InstallMethod( JacobianOfNumerators,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasNumerators ],
        
  function( LD )
    local jac, rel;
    
    jac := OriginalJacobianOfNumerators( LD );
    
    jac := MatrixOfMomenta( LD ) * jac;
    
    if NrColumns( jac ) = 0 then
        return jac;
    fi;
    
    rel := RelationsMatrixOfExternalMomenta( LD );
    
    jac := List( [ 1 .. NrColumns( jac ) ],
                   j -> DecideZeroRows( CertainColumns( jac, [ j ] ), rel ) );
    
    return UnionOfColumns( jac );
    
end );

##
InstallMethod( PairOfOriginalMatricesOfLoopDiagram,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators ],
        
  function( LD )
    
    return [ OriginalJacobianOfPropagators( LD ),
             HomalgDiagonalMatrix( Propagators( LD ) ) ];
    
end );

##
InstallMethod( PairOfMatricesOfLoopDiagram,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators ],
        
  function( LD )
    
    return [ JacobianOfPropagators( LD ),
             HomalgDiagonalMatrix( Propagators( LD ) ) ];
    
end );

#########################
#
# Operations
#
#########################

##
InstallMethod( SetAbbreviation,
        [ IsHomalgRingElement, IsString ],
        
  function( xy, str )

    xy!.Abbreviation := str;
    
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
    
    if HasName( x ) and HasName( y ) then
        xy!.Abbreviation := Concatenation( Name( x ), Name( y ) );
    fi;
    
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
InstallMethod( RingOfExtraLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasExtraLorentzInvariants and HasPropagators ],
        
  function( LD )
    local I, M, abbreviation, symbol, invariants, R, symbols;
    
    I := ExtraLorentzInvariants( LD );
    
    M := Length( I );
    
    abbreviation := LD!.abbreviation;
    
    symbol := LD!.LorentzSymbol;
    
    invariants := List( [ 1 .. M ],
                        function( i )
                          if IsBound( I[i]!.Abbreviation ) and abbreviation then
                              return I[i]!.Abbreviation;
                          fi;
                          return Concatenation( symbol, String( i ) );
                      end );
    
    R := CoefficientsRing( UnderlyingRing( LD ) ) * invariants;
    
    symbols := List( [ 1 .. Length( Propagators( LD ) ) ], i -> Concatenation( LOOP_INTEGRALS.DualSymbol, String( i ) ) );
    
    return R * symbols;
    
end );

##
InstallMethod( ReductionMatrixOfExtraLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasExtraLorentzInvariants ],
        
  function( LD )
    local R, invariants, I, M, red;
    
    R := RingOfExtraLorentzInvariants( LD );
    
    invariants := Indeterminates( BaseRing( R ) );
    
    R := R * Concatenation( List( ExternalMomenta( LD ), p -> p!.symbols ) );
    
    R := PolynomialRingWithProductOrdering( R );
    
    invariants := List( invariants, i -> i / R );
    
    I := ExtraLorentzInvariants( LD );
    
    M := Length( I );
    
    invariants := ListN( I, invariants, {a,b} -> a / R - b );
    
    invariants := HomalgMatrix( invariants, M, 1, R );
    
    red := UnionOfRows( invariants, R * RelationsMatrixOfExternalMomenta( LD ) );
    
    return BasisOfRows( red );
    
end );

##
InstallMethod( ExpressInExtraLorentzInvariants,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfExternalMomenta and HasExtraLorentzInvariants and HasPropagators ],
        
  function( mat, LD )
    local red, R, col;
    
    red := ReductionMatrixOfExtraLorentzInvariants( LD );
    
    R := HomalgRing( red );
    
    col := NrColumns( mat );
    
    mat := R * mat;
    
    if NrColumns( mat ) = 0 then
        return mat;
    fi;
    
    mat := List( [ 1 .. col ],
                 j -> DecideZeroRows( CertainColumns( mat, [ j ] ), red ) );
    
    return UnionOfColumns( mat );
    
end );

##
InstallMethod( SymanzikPolynomials,
        [ IsLoopDiagram and HasPropagators ],
        
  function( LD )
    local taylor, U, R, F;
    
    taylor := OriginalTaylorOfPropagators( LD );
    
    U := Determinant( taylor[1] );
    
    R := HomalgRing( U );
    
    F := ( ( TransposedMatrix( taylor[2] ) * AdjunctMatrix( taylor[1] ) * taylor[2] )[1,1] + U * taylor[3] );
    
    F := HomalgMatrix( [ F ], 1, 1, R );
    
    F := ExpressInExtraLorentzInvariants( F, LD );
    
    if not IsEmpty( ExternalMomenta( LD ) ) then
        R := BaseRing( R );
    fi;
    
    return [ U / R, F[1,1] / R ];
    
end );

##
InstallMethod( RingOfIndependentLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasIndependentLorentzInvariants ],
        
  function( LD )
    local I, M, abbreviation, symbol, invariants;
    
    I := IndependentLorentzInvariants( LD );
    
    M := Length( I );
    
    abbreviation := LD!.abbreviation;
    
    symbol := LD!.LorentzSymbol;
    
    invariants := List( [ 1 .. M ],
                        function( i )
                          if IsBound( I[i]!.Abbreviation ) and abbreviation then
                              return I[i]!.Abbreviation;
                          fi;
                          return Concatenation( symbol, String( i ) );
                      end );
    
    return CoefficientsRing( UnderlyingRing( LD ) ) * invariants;
    
end );

##
InstallMethod( ReductionMatrixOfIndependentLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasIndependentLorentzInvariants ],
        
  function( LD )
    local R, invariants, I, M, red;
    
    R := RingOfIndependentLorentzInvariants( LD );
    
    invariants := Indeterminates( R );
    
    R := R * List( Indeterminates( UnderlyingRing( LD ) ), String );
    
    R := PolynomialRingWithProductOrdering( R );
    
    invariants := List( invariants, i -> i / R );
    
    I := IndependentLorentzInvariants( LD );
    
    M := Length( I );
    
    invariants := ListN( I, invariants, {a,b} -> a / R - b );
    
    invariants := HomalgMatrix( invariants, M, 1, R );
    
    red := UnionOfRows( invariants, R * RelationsMatrixOfExternalMomenta( LD ) );
    
    return BasisOfRows( red );
    
end );

##
InstallMethod( ExpressInIndependentLorentzInvariants,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfExternalMomenta and HasIndependentLorentzInvariants ],
        
  function( mat, LD )
    local red, R, col;
    
    red := ReductionMatrixOfIndependentLorentzInvariants( LD );
    
    R := HomalgRing( red );
    
    col := NrColumns( mat );
    
    mat := R * mat;
    
    if NrColumns( mat ) = 0 then
        return mat;
    fi;
    
    mat := List( [ 1 .. col ],
                 j -> DecideZeroRows( CertainColumns( mat, [ j ] ), red ) );
    
    return UnionOfColumns( mat );
    
end );

##
InstallMethod( RingOfPropagatorsAndNumeratorsAndExtraLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    local symbolD, D, N, propagators, symbolN, Z, A, numerators,
          symbolK, abbreviation, K, M, invariants, indets;
    
    symbolD := LD!.PropagatorSymbol;
    
    D := Propagators( LD );
    
    N := Length( D );
    
    propagators := List( [ 1 .. N ], i -> Concatenation( symbolD, String( i ) ) );
    
    symbolN := LD!.NumeratorSymbol;
    
    Z := Numerators( LD );
    
    A := Length( Z );
    
    numerators := List( [ N + 1 .. N + A ], i -> Concatenation( symbolN, String( i ) ) );
    
    K := ExtraLorentzInvariants( LD );
    
    M := Length( K );
    
    symbolK := LD!.LorentzSymbol;
    
    abbreviation := LD!.abbreviation;
    
    invariants := List( [ 1 .. M ],
                        function( i )
                          if IsBound( K[i]!.Abbreviation ) and abbreviation then
                              return K[i]!.Abbreviation;
                          fi;
                          return Concatenation( symbolK, String( N + A + i ) );
                      end );
    
    indets := Concatenation( propagators, numerators );
    
    return CoefficientsRing( UnderlyingRing( LD ) ) * invariants * indets;
    
end );

##
InstallMethod( RingOfLoopDiagram,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    local K, D, Z, R, indets, invariants, propagators, numerators, masses, S;
    
    K := ExtraLorentzInvariants( LD );
    D := Propagators( LD );
    Z := Numerators( LD );
    
    R := RingOfPropagatorsAndNumeratorsAndExtraLorentzInvariants( LD );
    
    indets := Indeterminates( R );
    indets := List( indets, String );
    
    invariants := indets{[ 1 .. Length( K ) ]};
    propagators := indets{[ Length( K ) + 1 .. Length( K ) + Length( D ) ]};
    numerators := indets{[ Length( K ) + Length( D ) + 1 .. Length( indets ) ]};
    
    indets := Concatenation( propagators, numerators );
    
    masses := List( LD!.masses, String );
    
    invariants := Concatenation( masses, [ LD!.DimensionSymbol ], invariants );
    
    if LD!.rational then
        
        invariants := JoinStringsWithSeparator( invariants );
        
        S := LOOP_INTEGRALS.ConstructorOfDefaultField( invariants, R );
        
        S := S * indets;
        
    else
        
        S := LOOP_INTEGRALS.ConstructorOfDefaultField( R ) * invariants;
        
        S := S * indets;
        
    fi;
    
    return S;
    
end );

##
InstallMethod( ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    local K, D, Z, R, indets, red, S;
    
    K := ExtraLorentzInvariants( LD );
    D := Propagators( LD );
    Z := Numerators( LD );
    
    R := RingOfPropagatorsAndNumeratorsAndExtraLorentzInvariants( LD );
    
    indets := Indeterminates( R );
    indets := List( indets, String );
    
    R := R * List( Indeterminates( UnderlyingRing( LD ) ), String );
    
    R := PolynomialRingWithProductOrdering( R );
    
    indets := List( indets, p -> p / R );
    
    R!.MatrixOfPropagatorsAndNumerators := HomalgMatrix( indets{[ Length( K ) + 1 .. Length( indets ) ]}, 1, Length( D ) + Length( Z ), R );
    
    indets := ListN( Concatenation( K, D, Z ), indets, {a,b} -> a / R - b );
    
    indets := HomalgMatrix( indets, Length( indets ), 1, R );
    
    red := UnionOfRows( indets, R * RelationsMatrixOfExternalMomenta( LD ) );
    
    S := RingOfLoopDiagram( LD );
    
    return red;
    
end );

##
InstallMethod( ExpressInPropagatorsAndNumeratorsAndExtraLorentzInvariants,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( mat, LD )
    local red, R, col, S, Smat;
    
    red := ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants( LD );
    
    R := HomalgRing( red );
    
    col := NrColumns( mat );
    
    S := RingOfLoopDiagram( LD );
    
    if NrColumns( mat ) = 0 then
        return S * mat;
    fi;
    
    mat := R * mat;
    
    mat := List( [ 1 .. col ],
                 j -> DecideZeroRows( CertainColumns( mat, [ j ] ), red ) );
    
    mat := UnionOfColumns( mat );
    
    Smat := S * mat;
    
    if not R * Smat = mat then
        Error( "expressing the matrix mat in terms of the propagators, numerators, and extra Lorentz invariants failed\n" );
    fi;
    
    return Smat;
    
end );

##
InstallMethod( JacobianOfPropagatorsInIndependentLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasIndependentLorentzInvariants ],
        
  function( LD )
    
    return ExpressInIndependentLorentzInvariants(
                   JacobianOfPropagators( LD ), LD );
    
end );

##
InstallMethod( JacobianOfPropagatorsInPropagators,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    
    return ExpressInPropagatorsAndNumeratorsAndExtraLorentzInvariants(
                   JacobianOfPropagators( LD ), LD );
    
end );

##
InstallMethod( JacobianOfNumeratorsInPropagators,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    
    return ExpressInPropagatorsAndNumeratorsAndExtraLorentzInvariants(
                   JacobianOfNumerators( LD ), LD );
    
end );

##
InstallMethod( PairOfMatricesOfLoopDiagramInIndependentLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasIndependentLorentzInvariants ],
        
  function( LD )

    return [ JacobianOfPropagatorsInIndependentLorentzInvariants( LD ),
             ExpressInIndependentLorentzInvariants(
                     HomalgDiagonalMatrix( Propagators( LD ) ), LD ) ];
end );

##
InstallMethod( PairOfMatricesOfLoopDiagramInPropagators,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )

    return [ JacobianOfPropagatorsInPropagators( LD ),
             ExpressInPropagatorsAndNumeratorsAndExtraLorentzInvariants(
                     HomalgDiagonalMatrix( Propagators( LD ) ), LD ) ];
    
end );

##
InstallMethod( JacobianOfLoopDiagramInPropagators,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    
    return UnionOfColumns(
                   JacobianOfPropagatorsInPropagators( LD ),
                   JacobianOfNumeratorsInPropagators( LD ) );
    
end );

##
InstallMethod( JacobianOfCoefficientsVectorInPropagators,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( vec, LD )
    local red, R, S, indets;
    
    red := ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants( LD );
    
    R := HomalgRing( red );
    
    indets := R!.MatrixOfPropagatorsAndNumerators;
    
    S := RingOfLoopDiagram( LD );
    
    indets := S * indets;
    
    return Diff( indets, Involution( vec ) );
    
end );

##
InstallMethod( DivergenceOfCoefficientsVectorOfLoopDiagram,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( vec, LD )
    local l, k, jacLD, jacLV, trace, sum, D;
    
    l := Length( LoopMomenta( LD ) );
    
    k := Length( ExternalMomenta( LD ) );
    
    Assert( 0, NrRows( vec ) = 1 and NrColumns( vec ) = l * ( l + k ) );
    
    jacLD := JacobianOfLoopDiagramInPropagators( LD );
    jacLV := JacobianOfCoefficientsVectorInPropagators( vec, LD );
    
    trace := Sum( [ 1 .. NrRows( jacLD ) ], i -> ( jacLD[i] * Involution( jacLV[i] ) )[1,1] );
    
    sum := Sum( [ 1 .. l ], p -> vec[1, (p - 1) * (l + k) + p] );
    
    D := ( LD!.DimensionSymbol / HomalgRing( sum ) );
    
    return trace + D * sum;
    
end );

##
InstallMethod( DoubleShiftAlgebra,
        [ IsHomalgRing ],
        
  function( R )
    local Ds, D_s, c, exponents, B, A, shifts, pairs, Y;
    
    if IsBound( R!.DoubleShiftAlgebra ) then
        return R!.DoubleShiftAlgebra;
    fi;
    
    Ds := RelativeIndeterminatesOfPolynomialRing( R );
    
    Ds := List( Ds, String );
    
    D_s := List( Ds, D -> Concatenation( D, "_" ) );
    
    c := Length( Ds );
    
    exponents := List( [ 1 .. c ], i -> Concatenation( LOOP_INTEGRALS.ExponentSymbol, String( i ) ) );
    
    B := BaseRing( R );
    
    A := B * JoinStringsWithSeparator( exponents );
    
    if IsIdenticalObj( ValueOption( "pairs" ), false ) then
        shifts := Concatenation( Ds, D_s );
        pairs := false;
    else
        shifts := Concatenation( ListN( Ds, D_s, {d, d_} -> [ d, d_ ] ) );
        pairs := true;
    fi;
    
    Y := DoubleShiftAlgebra( A, shifts : steps := -1, pairs := pairs );
    
    Y!.Ds := Ds;
    Y!.D_s := D_s;
    
    AmbientRing( Y )!.Ds := Ds;
    AmbientRing( Y )!.D_s := D_s;
    
    R!.DoubleShiftAlgebra := Y;
    
    return Y;
    
end );

##
InstallMethod( RationalDoubleShiftAlgebra,
        [ IsHomalgRing ],
        
  function( R )
    local Q, r, Ds, D_s, c, exponents, B, A, shifts, pairs, Y;
    
    if IsBound( R!.RationalDoubleShiftAlgebra ) then
        return R!.RationalDoubleShiftAlgebra;
    fi;
    
    Q := HomalgFieldOfRationalsInMaple();
    
    B := Q * List( Indeterminates( BaseRing( R ) ), String );
    
    Ds := RelativeIndeterminatesOfPolynomialRing( R );
    
    Ds := List( Ds, String );
    
    D_s := List( Ds, D -> Concatenation( D, "_" ) );
    
    c := Length( Ds );
    
    exponents := List( [ 1 .. c ], i -> Concatenation( LOOP_INTEGRALS.ExponentSymbol, String( i ) ) );
    
    A := B * JoinStringsWithSeparator( exponents );
    
    if IsIdenticalObj( ValueOption( "pairs" ), false ) then
        shifts := Concatenation( Ds, D_s );
        pairs := false;
    else
        shifts := Concatenation( ListN( Ds, D_s, {d, d_} -> [ d, d_ ] ) );
        pairs := true;
    fi;
    
    Y := RationalDoubleShiftAlgebra( A, shifts : steps := -1, pairs := pairs );
    
    Y!.Ds := Ds;
    Y!.D_s := D_s;
    
    AmbientRing( Y )!.Ds := Ds;
    AmbientRing( Y )!.D_s := D_s;
    
    R!.RationalDoubleShiftAlgebra := Y;
    
    return Y;
    
end );

##
InstallMethod( AssociatedWeylAlgebra,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    local R, c, oper;
    
    ## Q[m,s,D][a1,...,a_s]<D1,D1_,...,Ds,Ds_>
    ## IBPRelation is a construction which is linear only over Q[m,s,D]
    R := RingOfLoopDiagram( LD );
    
    c := Length( RelativeIndeterminatesOfPolynomialRing( R ) );
    
    oper := List( [ 1 .. c ], i -> Concatenation( LOOP_INTEGRALS.WeylSymbol, String( i ) ) );
    
    return RingOfDerivations( R, oper );
    
end );

##
InstallMethod( IBPRelation,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( vec, LD )
    local R, Y, exponents, c, D_s, oper, div, jacLD;
    
    ## Q[D,s12,s14][D1,D2,D3,D4,D5,D6,D7,N8,N9]
    R := HomalgRing( vec );
    
    ## Q[D,s12,s14][a1,a2,a3,a4,a5,a6,a7,a8,a9]<D1,D1_,D2,D2_,D3,D3_,D4,D4_,D5,D5_,D6,D6_,D7,D7_,N8,N8_,N9,N9_>/( D1*D1_-1, D2*D2_-1, D3*D3_-1, D4*D4_-1, D5*D5_-1, D6*D6_-1, D7*D7_-1, N8*N8_-1, N9*N9_-1 )
    Y := DoubleShiftAlgebra( R );
    
    ## [ "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9" ]
    exponents := RelativeIndeterminatesOfPolynomialRing( BaseRing( Y ) );
    exponents := List( exponents, String );
    
    ## n + z
    c := Length( exponents );
    
    ## [ "D1_", "D2_", "D3_", "D4_", "D5_", "D6_", "D7_", "N8_", "N9_" ]
    D_s := IndeterminateShiftsOfDoubleShiftAlgebra( Y ){List( [ 1 .. c ], i -> 2 * i )};
    D_s := List( D_s, String );
    
    oper := List( [ 1 .. c ], i -> Concatenation( exponents[i], "*", D_s[i] ) );
    
    oper := Concatenation( "[", JoinStringsWithSeparator( oper ), "]" );
    
    ## -[ ... -a_i D_i ... -a_{j+n} N_{j+n} ... | i = 1 .. L, j = 1 .. E ]
    oper := -HomalgMatrix( oper, 1, c, Y );
    
    div := DivergenceOfCoefficientsVectorOfLoopDiagram( vec, LD );
    
    jacLD := JacobianOfLoopDiagramInPropagators( LD );
    
    ## * Usually we would act on the space of integrals as usual from the left.
    ## * This action is only linear w.r.t. the constants (= BaseRing( R ) ).
    ## * The IBP relations would then generate a right ideal over Y which is much larger
    ##   than the set of IBP relations parametrized by vec over R.
    ## * However, the GB engines only support GB of left ideals in the noncommutative setup.
    ## * This forces us to act from the right in order for the above mentioned ideal
    ##   to be a left ideal.
    return div / Y + ( oper * TransposedMatrix( Y * ( vec * jacLD ) ) )[1,1];
    
end );

##
InstallMethod( IBPRelation,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants, IsList ],
        
  function( vec, LD, exponents )
    local ibp, Y, P, S, T, phi, c;
    
    ibp := IBPRelation( vec, LD );
    
    Y := HomalgRing( ibp );

    P := AmbientRing( Y );
    
    S := P!.CommutativeDoubleShiftAlgebra;
    T := P!.LaurentAlgebra;
    
    ibp := ibp / S;
    
    phi := Concatenation(
                   [ Indeterminates( BaseRing( T ) ),
                     exponents,
                     RelativeIndeterminatesOfPolynomialRing( T ) ] );
    
    phi := List( phi, String );
    
    c := Length( phi );
    
    phi := Concatenation( "[", JoinStringsWithSeparator( phi ), "]" );
    
    phi := HomalgMatrix( phi, 1, c, T );
    
    phi := RingMap( phi, AmbientRing( S ), T );
    
    ## specify the exponents
    return Pullback( phi, ibp );
    
end );

##
InstallMethod( MatrixOfIBPRelations,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( mat, LD )
    local ibps;
    
    ibps := List( [ 1 .. NrRows( mat ) ], i -> IBPRelation( mat[i], LD ) );
    
    return HomalgMatrix( ibps, Length( ibps ), 1, HomalgRing( ibps[1] ) );
    
end );

##
InstallMethod( MatrixOfIBPRelations,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants, IsList ],
        
  function( mat, LD, exponents )
    local ibps;
    
    ibps := List( [ 1 .. NrRows( mat ) ], i -> IBPRelation( mat[i], LD, exponents ) );
    
    return HomalgMatrix( ibps, Length( ibps ), 1, HomalgRing( ibps[1] ) );
    
end );

##
InstallMethod( MatrixOfIBPRelations,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    local id;
    
    id := HomalgIdentityMatrix( DimensionOfCoefficientsVector( LD ), RingOfLoopDiagram( LD ) );
    
    return MatrixOfIBPRelations( id, LD );
    
end );

##
InstallMethod( MatrixOfIBPRelations,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants, IsList ],
        
  function( LD, exponents )
    local id;
    
    id := HomalgIdentityMatrix( DimensionOfCoefficientsVector( LD ), RingOfLoopDiagram( LD ) );
    
    return MatrixOfIBPRelations( id, LD, exponents );
    
end );

##
InstallMethod( BasisOfIBPRelations,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )

    return BasisOfRows( MatrixOfIBPRelations( LD ) );
    
end );

##
InstallMethod( MatrixOfSpecialIBPRelations,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    local syz;
    
    ## Q[m,s,D][D1,...,Ds]
    syz := SyzygiesOfRows( PairOfMatricesOfLoopDiagramInPropagators( LD ) );
    
    #syz := ReducedBasisOfRowModule( syz );
    
    return MatrixOfIBPRelations( syz, LD );
    
end );

##
InstallMethod( MatrixOfSpecialIBPRelations,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants, IsList ],
        
  function( LD, exponents )
    local syz;
    
    syz := SyzygiesOfRows( PairOfMatricesOfLoopDiagramInPropagators( LD ) );
    
    return MatrixOfIBPRelations( syz, LD, exponents );
    
end );

##
InstallMethod( BasisOfSpecialIBPRelations,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    
    return BasisOfRows( MatrixOfSpecialIBPRelations( LD ) );
    
end );

##
InstallMethod( IBPRelationInWeylAlgebra,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( vec, LD )
    local W, oper, div, jacLD;
    
    W := AssociatedWeylAlgebra( LD );
    
    oper := IndeterminateDerivationsOfRingOfDerivations( W );
    
    oper := -HomalgMatrix( oper, 1, Length( oper ), W );
    
    div := DivergenceOfCoefficientsVectorOfLoopDiagram( vec, LD );
    
    jacLD := JacobianOfLoopDiagramInPropagators( LD );
    
    return div / W + ( oper * TransposedMatrix( W * ( vec * jacLD ) ) )[1,1];
    
end );

##
InstallMethod( MatrixOfIBPRelationsInWeylAlgebra,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( mat, LD )
    local ibps;
    
    ibps := List( [ 1 .. NrRows( mat ) ], i -> IBPRelationInWeylAlgebra( mat[i], LD ) );
    
    return HomalgMatrix( ibps, Length( ibps ), 1, HomalgRing( ibps[1] ) );
    
end );

##
InstallMethod( MatrixOfIBPRelationsInWeylAlgebra,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    local id;
    
    id := HomalgIdentityMatrix( DimensionOfCoefficientsVector( LD ), RingOfLoopDiagram( LD ) );
    
    return MatrixOfIBPRelationsInWeylAlgebra( id, LD );
    
end );

##
InstallMethod( BasisOfIBPRelationsInWeylAlgebra,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )

    return BasisOfRows( MatrixOfIBPRelationsInWeylAlgebra( LD ) );
    
end );

##
InstallMethod( MatrixOfSpecialIBPRelationsInWeylAlgebra,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    local syz, id, ibps;
    
    syz := SyzygiesOfRows( PairOfMatricesOfLoopDiagramInPropagators( LD ) );
    
    ibps := List( [ 1 .. NrRows( syz ) ], i -> IBPRelationInWeylAlgebra( syz[i], LD ) );
    
    return HomalgMatrix( ibps, Length( ibps ), 1, HomalgRing( ibps[1] ) );
    
end );

##
InstallMethod( BasisOfSpecialIBPRelationsInWeylAlgebra,
        [ IsLoopDiagram and HasRelationsOfExternalMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    
    return BasisOfRows( MatrixOfSpecialIBPRelationsInWeylAlgebra( LD ) );
    
end );

##
InstallMethod( SymanzikPolynomials,
        [ IsLoopDiagram and HasPropagators, IsList ],
        
  function( LD, list_of_ones )
    local UF, R, var, Tvar, B, T, map, Bvar;
    
    UF := SymanzikPolynomials( LD );
    
    R := HomalgRing( UF[1] );
    
    if HasRelativeIndeterminatesOfPolynomialRing( R ) then
        var := RelativeIndeterminatesOfPolynomialRing( R );
    else
        var := Indeterminates( R );
    fi;
    
    Tvar := var{list_of_ones};
    
    B := BaseRing( R );
    
    T := B * List( Tvar, String );
    
    map := ListWithIdenticalEntries( Length( var ), 0 );
    
    map{list_of_ones} := Tvar;
    
    Bvar := Indeterminates( B );
    
    map := Concatenation( Bvar, map );
    
    map := List( map, a -> a / T );
    
    map := RingMap( map, R, T );
    
    return List( UF, s -> Pullback( map, s ) );
    
end );

##
InstallMethod( DegreesOfMonomialsOfProductOfSymanzikPolynomials,
        [ IsLoopDiagram and HasPropagators, IsList ],
        
  function( LD, list_of_ones )
    local s, degree, R, l, b, monomials;
    
    s := Product( SymanzikPolynomials( LD, list_of_ones ) );
    
    R := HomalgRing( s );
    
    l := Length( list_of_ones );
    b := Length( Indeterminates( BaseRing( R ) ) );
    
    degree := ListWithIdenticalEntries( b + l, ListWithIdenticalEntries( l, 0 ) );
    
    degree{b + [ 1 .. l ]} := IdentityMat( l );
    
    degree := DegreeOfRingElementFunction( HomalgRing( s ), degree );
    
    monomials := Coefficients( s )!.monomials;
    
    return List( monomials, degree );
    
end );

##
InstallMethod( IsScalelessLoopIntegral,
        [ IsLoopDiagram and HasPropagators, IsList ],
        
  function( LD, list_of_ones )
    local degrees, v1, l, rank;
    
    if IsEmpty( list_of_ones ) then
        return true;
    fi;
    
    degrees := DegreesOfMonomialsOfProductOfSymanzikPolynomials( LD, list_of_ones );
    
    if IsEmpty( degrees ) then
        return true;
    fi;
    
    v1 := degrees[1];
    
    l := Length( degrees );
    
    if l = 1 then
        degrees := [ v1 ];
    else
        degrees := List( degrees{[ 2 .. Length( degrees ) ]}, v -> v1 - v );
    fi;
    
    rank := Rank( degrees );
    
    return rank < Length( v1 ) - 1 and rank > 0;
    
end );

##
InstallMethod( GeneratorsOfScalelessSectors,
        [ IsLoopDiagram and HasPropagators ],
        
  function( LD )
    local m, n, generators, k, iter, comb, Y, shifts, mults;
    
    m := Length( Propagators( LD ) );
    
    n := [ 1 .. m ];
    
    generators := [ ];
    
    for k in Reversed( [ 0 .. Length( n ) - 1 ] ) do
        
        iter := IteratorOfCombinations( n, k );
        
        for comb in iter do
            
            if ForAny( generators, gen -> IsSubset( gen, comb ) ) then
                continue;
            fi;
            
            if IsScalelessLoopIntegral( LD, comb ) then
                Add( generators, comb );
            fi;
            
        od;
        
    od;
    
    shifts :=
      Concatenation(
              List( n, i -> Concatenation( LOOP_INTEGRALS.PropagatorSymbol, String( i ) ) ),
              List( Length( n ) + [ 1 .. Length( Numerators( LD ) ) ], i -> Concatenation( LOOP_INTEGRALS.NumeratorSymbol, String( i ) ) ) );
    
    n := [ 1 .. Length( shifts ) ];
    
    Y := DoubleShiftAlgebra( RingOfLoopDiagram( LD ) );
    
    Y := AmbientRing( Y );
    
    shifts := List( shifts, D -> D / Y );
    
    generators := List( generators, gen -> Product( shifts{Difference( n , gen )} ) );
    
    return HomalgMatrix( generators, 1, Length( generators ), Y );
    
end );

##
InstallMethod( GeneratorsOfScalelessSectors,
        [ IsLoopDiagram and HasPropagators, IsList ],
        
  function( LD, exponents )
    local n, shifts, Y, a, gen;
    
    if not Length( exponents ) = Length( Propagators( LD ) ) + Length( Numerators( LD ) ) then
        Error( "the length of the list exponents must be equal to the sum of the number of propagators and numerators\n" );
    fi;
    
    exponents := List( exponents, a -> a - 1 );
    
    n := [ 1 .. Length( Propagators( LD ) ) ];
    
    shifts :=
      Concatenation(
              List( n, i -> Concatenation( LOOP_INTEGRALS.PropagatorSymbol, String( i ) ) ),
              List( Length( n ) + [ 1 .. Length( Numerators( LD ) ) ], i -> Concatenation( LOOP_INTEGRALS.NumeratorSymbol, String( i ) ) ) );
    
    gen := GeneratorsOfScalelessSectors( LD );
    
    Y := HomalgRing( gen );
    
    shifts := List( shifts, D -> D / Y );
    
    a := Product( ListN( shifts, exponents, {D,a} -> D^a ) );
    
    return a * gen;
    
end );

##
InstallMethod( GeneratorsOfScalelessSectorsInWeylAlgebra,
        [ IsLoopDiagram and HasPropagators ],
        
  function( LD )
    
    return AssociatedWeylAlgebra( LD ) * GeneratorsOfScalelessSectors( LD );
    
end );

##
InstallMethod( NormalForm,
        [ IsHomalgRingElement, IsHomalgMatrix ],
        
  function( operator, G )
    
    return DecideZero( operator, G );
    
end );

##
InstallMethod( NormalForm,
        [ IsHomalgRingElement, IsHomalgMatrix ],
        
  function( operator, G )
    local Y, P, NF;
    
    Y := HomalgRing( operator );
    
    P := AmbientRing( Y );
    
    if not IsHomalgExternalRingInMapleRep( P ) then
        TryNextMethod( );
    fi;
    
    NF := DecideZero( operator, G );
    
    NF := homalgSendBlocking(
                  [ "map(factor, collect( ", EvalRingElement( NF ), ", ",
                    String( IndeterminateShiftsOfRationalPseudoDoubleShiftAlgebra( P ) ), " ) )" ],
                  P, "define" );
    
    NF := HomalgExternalRingElement( NF, P );
    
    return NF / Y;
    
end );

##
InstallMethod( NormalForm,
        [ IsHomalgRingElement, IsHomalgMatrix, IsList ],
        
  function( operator, G, initial_integral )
    local Y, P, NF;
    
    Y := HomalgRing( operator );
    
    P := AmbientRing( Y );
    
    if not IsHomalgExternalRingInMapleRep( P ) then
        TryNextMethod( );
    fi;
    
    NF := NormalForm( operator, G );
    
    NF := homalgSendBlocking(
                  [ "map(factor, collect( simplify( subs( [ ",
                    JoinStringsWithSeparator(
                            ListN( RelativeParametersOfRationalPseudoDoubleShiftAlgebra( P ), initial_integral,
                                   { a, v } -> Concatenation( String( a ), " = ", String( v ), " " ) ) ),
                    "], ", EvalRingElement( NF ), " ) ), ",
                    String( IndeterminateShiftsOfRationalPseudoDoubleShiftAlgebra( P ) ), " ) )" ],
                  P, "define" );
    
    NF := HomalgExternalRingElement( NF, P );
    
    return NF / Y;
    
end );

##
InstallMethod( NormalFormWrtInitialIntegral,
        [ IsHomalgRingElement, IsHomalgMatrix ],
        
  function( operator, G )
    local P, initial_integral;
    
    P := AmbientRing( HomalgRing( operator ) );
    
    initial_integral :=
      Concatenation(
              ListWithIdenticalEntries( Length( Filtered( P!.Ds, Di -> Di{[1]} = LOOP_INTEGRALS.PropagatorSymbol ) ), 1 ),
              ListWithIdenticalEntries( Length( Filtered( P!.Ds, Di -> Di{[1]} = LOOP_INTEGRALS.NumeratorSymbol ) ), 0 ) );
    
    return NormalForm( operator, G, initial_integral );
    
end );

#########################
#
# View & Display methods
#
#########################

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
           " & external momenta ", ExternalMomenta( LD ) );
    
    if IsBound( LD!.masses ) and not IsEmpty( LD!.masses ) then
        Print( " & masses ", LD!.masses );
    fi;
    
    Print( ">" );
    
end );
