# SPDX-License-Identifier: GPL-2.0-or-later
# LoopIntegrals: Compute master integrals using commutative and noncommutative methods from computational algebraic geometry
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
DeclareAttribute( "UnderlyingMatrix",
        IsLorentzVector );

#!
DeclareAttribute( "LoopDiagram",
        IsLorentzVector );

#!
DeclareAttribute( "Components",
        IsLorentzVector );

#!
DeclareAttribute( "LoopMomenta",
        IsLoopDiagram );

#!
DeclareAttribute( "ExternalMomenta",
        IsLoopDiagram );

#!
DeclareAttribute( "DimensionOfCoefficientsVector",
        IsLoopDiagram );

#!
DeclareAttribute( "UnderlyingRing",
        IsLoopDiagram );

#!
DeclareAttribute( "RelationsOfExternalMomenta",
        IsLoopDiagram );

#!
DeclareAttribute( "RelationsMatrixOfExternalMomenta",
        IsLoopDiagram );

#!
DeclareAttribute( "IndependentLorentzInvariants",
        IsLoopDiagram );

#!
DeclareAttribute( "Propagators",
        IsLoopDiagram );

#!
DeclareAttribute( "Numerators",
        IsLoopDiagram );

#!
DeclareAttribute( "ExtraLorentzInvariants",
        IsLoopDiagram );

#!
DeclareAttribute( "OriginalJacobianOfPropagators",
        IsLoopDiagram );

#!
DeclareAttribute( "OriginalJacobianOfNumerators",
        IsLoopDiagram );

#!
DeclareAttribute( "OriginalTaylorOfPropagators",
        IsLoopDiagram );

#!
DeclareAttribute( "MatrixOfMomenta",
        IsLoopDiagram );

#!
DeclareAttribute( "JacobianOfPropagators",
        IsLoopDiagram );

#!
DeclareAttribute( "JacobianOfNumerators",
        IsLoopDiagram );

#!
DeclareAttribute( "PairOfOriginalMatricesOfLoopDiagram",
        IsLoopDiagram );

#!
DeclareAttribute( "PairOfMatricesOfLoopDiagram",
        IsLoopDiagram );

#!
DeclareAttribute( "SymanzikPolynomials",
        IsLoopDiagram );

#!
DeclareAttribute( "GeneratorsOfScalelessSectors",
        IsLoopDiagram );

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

#! @Group LoopDiagram
DeclareOperation( "LoopDiagram",
        [ IsList, IsList, IsInt ] );

#! @Group LoopDiagram
DeclareOperation( "LoopDiagram",
        [ IsList, IsList ] );

####################################
#
#! @Section Operations
#
####################################

#!
DeclareOperation( "SetAbbreviation",
        [ IsHomalgRingElement, IsString ] );
    
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
    
#!
DeclareOperation( "RingOfExtraLorentzInvariants",
        [ IsLoopDiagram ] );

#!
DeclareOperation( "ReductionMatrixOfExtraLorentzInvariants",
        [ IsLoopDiagram ] );

#!
DeclareOperation( "ExpressInExtraLorentzInvariants",
        [ IsHomalgMatrix, IsLoopDiagram ] );

#!
DeclareOperation( "RingOfIndependentLorentzInvariants",
        [ IsLoopDiagram ] );

#!
DeclareOperation( "ReductionMatrixOfIndependentLorentzInvariants",
        [ IsLoopDiagram ] );

#!
DeclareOperation( "ExpressInIndependentLorentzInvariants",
        [ IsHomalgMatrix, IsLoopDiagram ] );

#!
DeclareOperation( "RingOfPropagatorsAndNumeratorsAndExtraLorentzInvariants",
        [ IsLoopDiagram ] );

#!
DeclareOperation( "RingOfLoopDiagram",
        [ IsLoopDiagram ] );

#!
DeclareOperation( "ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants",
        [ IsLoopDiagram ] );

#!
DeclareOperation( "ExpressInPropagatorsAndNumeratorsAndExtraLorentzInvariants",
        [ IsHomalgMatrix, IsLoopDiagram ] );

#!
DeclareOperation( "JacobianOfPropagatorsInIndependentLorentzInvariants",
        [ IsLoopDiagram ] );

#!
DeclareOperation( "JacobianOfPropagatorsInPropagators",
        [ IsLoopDiagram ] );

#!
DeclareOperation( "JacobianOfNumeratorsInPropagators",
        [ IsLoopDiagram ] );

#!
DeclareOperation( "PairOfMatricesOfLoopDiagramInIndependentLorentzInvariants",
        [ IsLoopDiagram ] );

#!
DeclareOperation( "PairOfMatricesOfLoopDiagramInPropagators",
        [ IsLoopDiagram ] );

#!
DeclareOperation( "JacobianOfLoopDiagramInPropagators",
        [ IsLoopDiagram ] );

#! @Arguments vec, LD
DeclareOperation( "JacobianOfCoefficientsVectorInPropagators",
        [ IsHomalgMatrix, IsLoopDiagram ] );

#! @Arguments vec, LD
DeclareOperation( "DivergenceOfCoefficientsVectorOfLoopDiagram",
        [ IsHomalgMatrix, IsLoopDiagram ] );

#!
DeclareOperation( "DoubleShiftAlgebra",
        [ IsHomalgRing ] );

#! @Arguments vec, LD
DeclareOperation( "IBPRelation",
        [ IsHomalgMatrix, IsLoopDiagram ] );

#! @Arguments vec, LD, exponents
DeclareOperation( "IBPRelation",
        [ IsHomalgMatrix, IsLoopDiagram, IsList ] );

#! @Arguments mat, LD
DeclareOperation( "MatrixOfIBPRelations",
        [ IsHomalgMatrix, IsLoopDiagram ] );

#! @Arguments mat, LD, exponents
DeclareOperation( "MatrixOfIBPRelations",
        [ IsHomalgMatrix, IsLoopDiagram, IsList ] );

#!
DeclareOperation( "MatrixOfIBPRelations",
        [ IsLoopDiagram ] );

#! @Arguments LD, exponents
DeclareOperation( "MatrixOfIBPRelations",
        [ IsLoopDiagram, IsList ] );

#!
DeclareOperation( "BasisOfIBPRelations",
        [ IsLoopDiagram ] );

#!
DeclareOperation( "MatrixOfSpecialIBPRelations",
        [ IsLoopDiagram ] );

#! @Arguments LD, exponents
DeclareOperation( "MatrixOfSpecialIBPRelations",
        [ IsLoopDiagram, IsList ] );

#!
DeclareOperation( "BasisOfSpecialIBPRelations",
        [ IsLoopDiagram ] );

#! @Arguments LD, list_of_ones
DeclareOperation( "SymanzikPolynomials",
        [ IsLoopDiagram, IsList ] );

#! @Arguments LD, list_of_ones
DeclareOperation( "DegreesOfMonomialsOfProductOfSymanzikPolynomials",
        [ IsLoopDiagram, IsList ] );

#! @Arguments LD, list_of_ones
DeclareOperation( "IsScalelessLoopIntegral",
        [ IsLoopDiagram, IsList ] );

#! @Arguments LD, exponents
DeclareOperation( "GeneratorsOfScalelessSectors",
        [ IsLoopDiagram, IsList ] );

####################################
#
#! @Section Tools
#
####################################

