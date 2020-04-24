#
# LoopIntegrals: Loop integrals
#
# Declarations
#

#! @Chapter Loop integrals

####################################
#
#! @Section Global variables
#
####################################

#!
DeclareGlobalVariable( "LOOP_INTEGRALS" );

####################################
#
#! @Section GAP categories
#
####################################

DeclareCategory( "IsLorentzVector",
        IsAttributeStoringRep );

BindGlobal( "TheFamilyOfLorentzVectors",
        NewFamily( "TheFamilyOfLorentzVector" ) );

BindGlobal( "TheTypeLorentzVector",
        NewType( TheFamilyOfLorentzVectors,
                IsLorentzVector ) );

DeclareCategory( "IsLoopDiagram",
        IsAttributeStoringRep );

BindGlobal( "TheFamilyOfLoopDiagrams",
        NewFamily( "TheFamilyOfLoopDiagrams" ) );

BindGlobal( "TheTypeLoopDiagram",
        NewType( TheFamilyOfLoopDiagrams,
                IsLoopDiagram ) );

####################################
#
#! @Section Properties
#
####################################

#!
DeclareProperty( "IsLoopMomentum",
        IsLoopDiagram );

#!
DeclareProperty( "IsExternalMomentum",
        IsLoopDiagram );

####################################
#
#! @Section Attributes
#
####################################

#!
DeclareAttribute( "Dimension",
        IsLorentzVector );

#!
DeclareAttribute( "LoopDiagram",
        IsLorentzVector );

#!
DeclareAttribute( "LoopMomenta",
        IsLoopDiagram );

#!
DeclareAttribute( "ExternalMomenta",
        IsLoopDiagram );

#!
DeclareAttribute( "RelationsOfMomenta",
        IsLoopDiagram );

#!
DeclareAttribute( "RelationsMatrixOfMomenta",
        IsLoopDiagram );

#!
DeclareAttribute( "Propagators",
        IsLoopDiagram );

#!
DeclareAttribute( "Numerators",
        IsLoopDiagram );

#!
DeclareAttribute( "IndependetLorentzInvariants",
        IsLoopDiagram );

#!
DeclareAttribute( "ReductionMatrixOfIndependetLorentzInvariants",
        IsLoopDiagram );

#!
DeclareAttribute( "ExtraLorentzInvariants",
        IsLoopDiagram );

#!
DeclareAttribute( "ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants",
        IsLoopDiagram );

#!
DeclareAttribute( "PairOfOriginalMatricesOfLoopDiagram",
        IsLoopDiagram );

#!
DeclareAttribute( "PairOfMatricesOfLoopDiagram",
        IsLoopDiagram );

#!
DeclareAttribute( "PairOfMatricesOfLoopDiagramInLorentzInvariants",
        IsLoopDiagram );

#!
DeclareAttribute( "PairOfMatricesOfLoopDiagramInPropagators",
        IsLoopDiagram );

#!
DeclareAttribute( "UnderlyingRing",
        IsLoopDiagram );

#!
DeclareAttribute( "UnderlyingMatrix",
        IsLorentzVector );

#!
DeclareAttribute( "Components",
        IsLorentzVector );

####################################
#
#! @Section Constructors
#
####################################

#!
DeclareOperation( "LorentzVector",
        [ IsString, IsInt ] );

#!
DeclareOperation( "LorentzVector",
        [ IsHomalgMatrix ] );

#!
DeclareOperation( "LoopDiagram",
        [ IsList, IsList, IsInt ] );

####################################
#
#! @Section Operations
#
####################################

#!
DeclareOperation( "+",
        [ IsLorentzVector, IsLorentzVector ] );
    
#!
DeclareOperation( "-",
        [ IsLorentzVector, IsLorentzVector ] );
    
#!
DeclareAttribute( "AdditiveInverse",
        IsLorentzVector );
    
#!
DeclareOperation( "[]",
        [ IsLorentzVector, IsInt ] );
    
#!
DeclareOperation( "*",
        [ IsLorentzVector, IsLorentzVector ] );
    
#!
DeclareOperation( "*",
        [ IsRingElement, IsLorentzVector ] );
    
#!
DeclareOperation( "*",
        [ IsLorentzVector, IsRingElement ] );
    
DeclareOperation( "POW",
        [ IsLorentzVector, IsInt ] );
    
####################################
#
#! @Section Tools
#
####################################

