#
# LoopIntegrals: Loop integrals
#
# Implementations
#

InstallValue( LOOP_INTEGRALS,
        rec(
            ConstructorOfDefaultField := HomalgFieldOfRationalsInSingular,
            LorentzSymbol := "x",
            Dimension := 1,
            DimensionSymbol := "D",
            PropagatorSymbol := "D",
            NumeratorSymbol := "N",
            ExponentSymbol := "a",
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
    local LD, LorentzVectors, parameters, C, momenta, R;
    
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
    
    parameters := ValueOption( "parameters" );
    
    if not IsStringRep( parameters ) then
        parameters := "";
    fi;
    
    LD!.parameters := parameters;
    
    if IsIdenticalObj( ValueOption( "C" ), fail ) then
        if parameters = "" then
            C := LOOP_INTEGRALS.ConstructorOfDefaultField( );
        else
            C := LOOP_INTEGRALS.ConstructorOfDefaultField( parameters );
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
    
    Perform( SplitString( parameters, "," ),
      function( m )
        local name;
        if IsBoundGlobal( m ) then
            MakeReadWriteGlobal( m );
            UnbindGlobal( m );
        fi;
        BindGlobal( m, m / R );
    end );
    
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
    
##
InstallMethod( DimensionOfCoefficientsVector,
        [ IsLoopDiagram and HasLoopMomenta and HasExternalMomenta ],
        
  function( LD )
    local L;
    
    L := Length( LoopMomenta( LD ) );
    
    return L * ( L + Length( ExternalMomenta( LD ) ) );
    
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
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators ],
        
  function( LD )
    local jac, rel;
    
    jac := OriginalJacobianOfPropagators( LD );

    jac := MatrixOfMomenta( LD ) * jac;
    
    rel := RelationsMatrixOfMomenta( LD );
    
    jac := List( [ 1 .. NrColumns( jac ) ],
                   j -> DecideZeroRows( CertainColumns( jac, [ j ] ), rel ) );
    
    return UnionOfColumns( jac );
    
end );

##
InstallMethod( JacobianOfNumerators,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasNumerators ],
        
  function( LD )
    local jac, rel;
    
    jac := OriginalJacobianOfNumerators( LD );
    
    jac := MatrixOfMomenta( LD ) * jac;
    
    if NrColumns( jac ) = 0 then
        return jac;
    fi;
    
    rel := RelationsMatrixOfMomenta( LD );
    
    jac := List( [ 1 .. NrColumns( jac ) ],
                   j -> DecideZeroRows( CertainColumns( jac, [ j ] ), rel ) );
    
    return UnionOfColumns( jac );
    
end );

##
InstallMethod( PairOfOriginalMatricesOfLoopDiagram,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators ],
        
  function( LD )
    
    return [ OriginalJacobianOfPropagators( LD ),
             HomalgDiagonalMatrix( Propagators( LD ) ) ];
    
end );

##
InstallMethod( PairOfMatricesOfLoopDiagram,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators ],
        
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
InstallMethod( RingOfIndependentLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasIndependentLorentzInvariants ],
        
  function( LD )
    local abbreviation, I, M, symbol, invariants, R;
    
    abbreviation := ValueOption( "abbreviation" );
    
    if not IsIdenticalObj( abbreviation, false ) then
        abbreviation := true;
    fi;
    
    if abbreviation then
        if IsBound( LD!.RingOfIndependentLorentzInvariants ) then
            return LD!.RingOfIndependentLorentzInvariants;
        fi;
    else
        if IsBound( LD!.RingOfIndependentLorentzInvariants_noabbreviation ) then
            return LD!.RingOfIndependentLorentzInvariants_noabbreviation;
        fi;
    fi;
    
    I := IndependentLorentzInvariants( LD );
    
    M := Length( I );
    
    symbol := ValueOption( "symbol" );
    
    if not IsStringRep( symbol ) then
        symbol := LOOP_INTEGRALS.LorentzSymbol;
    fi;
    
    invariants := List( [ 1 .. M ],
                        function( i )
                          if IsBound( I[i]!.Abbreviation ) and abbreviation then
                              return I[i]!.Abbreviation;
                          fi;
                          return Concatenation( symbol, String( i ) );
                      end );
    
    R := CoefficientsRing( UnderlyingRing( LD ) ) * invariants;
    
    if abbreviation then
        LD!.RingOfIndependentLorentzInvariants := R;
    else
        LD!.RingOfIndependentLorentzInvariants_noabbreviation := R;
    fi;
    
    return R;
    
end );

##
InstallMethod( ReductionMatrixOfIndependentLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasIndependentLorentzInvariants ],
        
  function( LD )
    local abbreviation, R, invariants, I, M, red;
    
    abbreviation := ValueOption( "abbreviation" );
    
    if not IsIdenticalObj( abbreviation, false ) then
        abbreviation := true;
    fi;
    
    if abbreviation then
        if IsBound( LD!.ReductionMatrixOfIndependentLorentzInvariants ) then
            return LD!.ReductionMatrixOfIndependentLorentzInvariants;
        fi;
    else
        if IsBound( LD!.ReductionMatrixOfIndependentLorentzInvariants_noabbreviation ) then
            return LD!.ReductionMatrixOfIndependentLorentzInvariants_noabbreviation;
        fi;
    fi;
    
    R := RingOfIndependentLorentzInvariants( LD );
    
    invariants := Indeterminates( R );
    
    R := R * List( Indeterminates( UnderlyingRing( LD ) ), String );
    
    R := PolynomialRingWithProductOrdering( R );
    
    invariants := List( invariants, i -> i / R );
    
    I := IndependentLorentzInvariants( LD );
    
    M := Length( I );
    
    invariants := ListN( I, invariants, {a,b} -> a / R - b );
    
    invariants := HomalgMatrix( invariants, M, 1, R );
    
    red := UnionOfRows( invariants, R * RelationsMatrixOfMomenta( LD ) );
    
    red := BasisOfRows( red );
    
    if abbreviation then
        LD!.ReductionMatrixOfIndependentLorentzInvariants := red;
    else
        LD!.ReductionMatrixOfIndependentLorentzInvariants_noabbreviation := red;
    fi;
    
    return red;
    
end );

##
InstallMethod( ExpressInIndependentLorentzInvariants,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfMomenta and HasIndependentLorentzInvariants ],
        
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
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    local symbolD, D, N, propagators, symbolN, Z, A, numerators,
          symbolK, K, M, invariants, indets, R;
    
    if IsBound( LD!.RingOfPropagatorsAndNumeratorsAndExtraLorentzInvariants ) then
        return LD!.RingOfPropagatorsAndNumeratorsAndExtraLorentzInvariants;
    fi;
    
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
    
    K := ExtraLorentzInvariants( LD );
    
    M := Length( K );
    
    symbolK := ValueOption( "symbolK" );
    
    if not IsStringRep( symbolK ) then
        symbolK := LOOP_INTEGRALS.LorentzSymbol;
    fi;
    
    invariants := List( [ 1 .. M ],
                        function( i )
                          if IsBound( K[i]!.Abbreviation ) then
                              return K[i]!.Abbreviation;
                          fi;
                          return Concatenation( symbolK, String( N + A + i ) );
                      end );
    
    indets := Concatenation( propagators, numerators );
    
    R := CoefficientsRing( UnderlyingRing( LD ) ) * invariants * indets;
    
    LD!.RingOfPropagatorsAndNumeratorsAndExtraLorentzInvariants := R;
    
    return R;
    
end );

##
InstallMethod( RingOfLoopDiagram,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    local rational, K, D, Z, R, indets, parameters, invariants, propagators, numerators, S;
    
    ## do not treat the extra Lorentz invariants as rational parameters
    ## as this slows down the syzygies computations in Singular significantly
    rational := IsIdenticalObj( ValueOption( "rational" ), true );
    
    if rational then
        if IsBound( LD!.RingOfLoopDiagram_rational ) then
            return LD!.RingOfLoopDiagram_rational;
        fi;
    else
        if IsBound( LD!.RingOfLoopDiagram ) then
            return LD!.RingOfLoopDiagram;
        fi;
    fi;
    
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
    
    parameters := SplitString( LD!.parameters, "," );
    
    invariants := Concatenation( parameters, [ LD!.DimensionSymbol ], invariants );
    
    if rational then
        
        invariants := JoinStringsWithSeparator( invariants );
        
        S := LOOP_INTEGRALS.ConstructorOfDefaultField( invariants, R );
        
        S := S * indets;
        
        LD!.RingOfLoopDiagram_rational := S;
        
    else
        
        S := LOOP_INTEGRALS.ConstructorOfDefaultField( R ) * invariants;
        
        S := S * indets;
        
        LD!.RingOfLoopDiagram := S;
        
    fi;
    
    return S;
    
end );

##
InstallMethod( ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    local rational, K, D, Z, R, indets, red, S;

    ## do not treat the extra Lorentz invariants as rational parameters
    ## as this slows down the syzygies computations in Singular significantly
    rational := IsIdenticalObj( ValueOption( "rational" ), true );
    
    if rational then
        if IsBound( LD!.ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants_rational ) then
            return LD!.ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants_rational;
        fi;
    else
        if IsBound( LD!.ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants ) then
            return LD!.ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants;
        fi;
    fi;
    
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
    
    red := UnionOfRows( indets, R * RelationsMatrixOfMomenta( LD ) );
    
    red := BasisOfRows( red );
    
    S := RingOfLoopDiagram( LD );
    
    if rational then
        LD!.ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants_rational := red;
        R!.RingAfterSuccessfulReduction_rational := S;
    else
        LD!.ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants := red;
        R!.RingAfterSuccessfulReduction := S;
    fi;
    
    return red;
    
end );

##
InstallMethod( ExpressInPropagatorsAndNumeratorsAndExtraLorentzInvariants,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( mat, LD )
    local rational, red, R, col, S, Smat;
    
    ## do not treat the extra Lorentz invariants as rational parameters
    ## as this slows down the syzygies computations in Singular significantly
    rational := IsIdenticalObj( ValueOption( "rational" ), true );
    
    red := ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants( LD );
    
    R := HomalgRing( red );
    
    col := NrColumns( mat );
    
    if rational then
        S := R!.RingAfterSuccessfulReduction_rational;
    else
        S := R!.RingAfterSuccessfulReduction;
    fi;
    
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
        [ IsLoopDiagram and HasRelationsOfMomenta and HasIndependentLorentzInvariants ],
        
  function( LD )
    
    return ExpressInIndependentLorentzInvariants(
                   JacobianOfPropagators( LD ), LD );
    
end );

##
InstallMethod( JacobianOfPropagatorsInPropagators,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    
    return ExpressInPropagatorsAndNumeratorsAndExtraLorentzInvariants(
                   JacobianOfPropagators( LD ), LD );
    
end );

##
InstallMethod( JacobianOfNumeratorsInPropagators,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    
    return ExpressInPropagatorsAndNumeratorsAndExtraLorentzInvariants(
                   JacobianOfNumerators( LD ), LD );
    
end );

##
InstallMethod( PairOfMatricesOfLoopDiagramInIndependentLorentzInvariants,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasIndependentLorentzInvariants ],
        
  function( LD )

    return [ JacobianOfPropagatorsInIndependentLorentzInvariants( LD ),
             ExpressInIndependentLorentzInvariants(
                     HomalgDiagonalMatrix( Propagators( LD ) ), LD ) ];
end );

##
InstallMethod( PairOfMatricesOfLoopDiagramInPropagators,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )

    return [ JacobianOfPropagatorsInPropagators( LD ),
             ExpressInPropagatorsAndNumeratorsAndExtraLorentzInvariants(
                     HomalgDiagonalMatrix( Propagators( LD ) ), LD ) ];
    
end );

##
InstallMethod( JacobianOfLoopDiagramInPropagators,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    
    return UnionOfColumns(
                   JacobianOfPropagatorsInPropagators( LD ),
                   JacobianOfNumeratorsInPropagators( LD ) );
    
end );

##
InstallMethod( JacobianOfCoefficientsVectorInPropagators,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( vec, LD )
    local red, R, S, indets;
    
    red := ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants( LD );
    
    R := HomalgRing( red );
    
    indets := R!.MatrixOfPropagatorsAndNumerators;
    
    S := R!.RingAfterSuccessfulReduction;
    
    indets := S * indets;
    
    return Diff( indets, Involution( vec ) );
    
end );

##
InstallMethod( DivergenceOfCoefficientsVectorOfLoopDiagram,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( vec, LD )
    local l, k, jacLD, jacLV, trace, sum, D;
    
    l := Length( LoopMomenta( LD ) );
    
    k := Length( ExternalMomenta( LD ) );
    
    Assert( 0, NrRows( vec ) = 1 and NrColumns( vec ) = l * ( l + k ) );
    
    jacLD := JacobianOfLoopDiagramInPropagators( LD );
    jacLV := JacobianOfCoefficientsVectorInPropagators( vec, LD );
    
    trace := Sum( [ 1 .. NrRows( jacLD ) ],  i -> ( jacLD[i] * Involution( jacLV[i] ) )[1,1] );
    
    sum := Sum( [ 1 .. l ], p -> vec[1, (p - 1) * (l + k) + p] );
    
    D := ( LD!.DimensionSymbol / HomalgRing( sum ) );
    
    return trace + D * sum;
    
end );

##
InstallMethod( DoubleShiftAlgebra,
        [ IsHomalgRing ],
        
  function( R )
    local Ds, D_s, c, exponents, B, A, shifts, Y;
    
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
    
    if IsIdenticalObj( ValueOption( "pairs" ), true ) then
        shifts := Concatenation( ListN( Ds, D_s, {d, d_} -> [ d, d_ ] ) );
    else
        shifts := Concatenation( Ds, D_s );
    fi;
    
    Y := DoubleShiftAlgebra( A, shifts : steps := -1 );
    
    R!.DoubleShiftAlgebra := Y;
    
    return Y;
    
end );

##
InstallMethod( IBPRelation,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( vec, LD )
    local R, Y, exponents, c, D_s, oper, S, div, jacLD;
    
    R := HomalgRing( vec );
    
    Y := DoubleShiftAlgebra( R : pairs := true );
    
    exponents := RelativeIndeterminatesOfPolynomialRing( BaseRing( Y ) );
    exponents := List( exponents, String );
    
    c := Length( exponents );
    
    D_s := IndeterminateShiftsOfDoubleShiftAlgebra( Y ){List( [ 1 .. c ], i -> 2 * i )};
    D_s := List( D_s, String );
    
    oper := List( [ 1 .. c ], i -> Concatenation( "-", exponents[i], "*", D_s[i]  ) );
    
    oper := Concatenation( "[", JoinStringsWithSeparator( oper ), "]" );
    
    S := Y!.CommutativeDoubleShiftAlgebra;
    
    oper := HomalgMatrix( oper, c, 1, S );
    
    div := DivergenceOfCoefficientsVectorOfLoopDiagram( vec, LD );
    
    jacLD := JacobianOfLoopDiagramInPropagators( LD );
    
    oper := DecideZero( ( div / S ) + ( ( S * vec ) * ( S * jacLD ) * oper )[1,1] );
    
    return oper / Y;
    
end );

##
InstallMethod( IBPRelation,
        [ IsHomalgMatrix, IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants, IsList ],
        
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
    
    return Pullback( phi, ibp );
    
end );

##
InstallMethod( MatrixOfIBPRelations,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    local id, ibps;
    
    id := HomalgIdentityMatrix( DimensionOfCoefficientsVector( LD ), RingOfLoopDiagram( LD ) );
    
    ibps := List( [ 1 .. NrRows( id ) ], i -> IBPRelation( id[i], LD ) );
    
    return HomalgMatrix( ibps, Length( ibps ), 1, HomalgRing( ibps[1] ) );

end );

##
InstallMethod( BasisOfIBPRelations,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )

    return BasisOfRows(  MatrixOfIBPRelations( LD ) );
    
end );

##
InstallMethod( MatrixOfSpecialIBPRelations,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    local syz, id, ibps;
    
    syz := SyzygiesOfRows( PairOfMatricesOfLoopDiagramInPropagators( LD ) );
    
    ibps := List( [ 1 .. NrRows( syz ) ], i -> IBPRelation( syz[i], LD ) );
    
    return HomalgMatrix( ibps, Length( ibps ), 1, HomalgRing( ibps[1] ) );
    
end );

##
InstallMethod( BasisOfSpecialIBPRelations,
        [ IsLoopDiagram and HasRelationsOfMomenta and HasPropagators and HasNumerators and HasExtraLorentzInvariants ],
        
  function( LD )
    
    return BasisOfRows( MatrixOfSpecialIBPRelations( LD ) );
    
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
           " & external momenta ", ExternalMomenta( LD ), ">" );
    
end );
