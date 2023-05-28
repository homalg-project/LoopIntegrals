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

#! @Arguments LD
DeclareAttribute( "LoopMomenta",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "ExternalMomenta",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "DimensionOfCoefficientsVector",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "UnderlyingRing",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "RelationsOfExternalMomenta",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "RelationsMatrixOfExternalMomenta",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "IndependentLorentzInvariants",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "Propagators",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "Numerators",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "ExtraLorentzInvariants",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "OriginalEulerMatrixOfPropagators",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "OriginalEulerMatrixOfNumerators",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "OriginalTaylorOfPropagators",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "MatrixOfMomenta",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "EulerMatrixOfPropagators",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "EulerMatrixOfNumerators",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "PairOfMatricesOfLoopDiagram",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "RingOfExtraLorentzInvariants",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "ReductionMatrixOfExtraLorentzInvariants",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "RingOfIndependentLorentzInvariants",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "ReductionMatrixOfIndependentLorentzInvariants",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "RingOfPropagatorsAndNumeratorsAndExtraLorentzInvariants",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "RingOfLoopDiagram",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "ReductionMatrixOfPropagatorsAndNumeratorsAndExtraLorentzInvariants",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "EulerMatrixOfPropagatorsInIndependentLorentzInvariants",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "EulerMatrixOfPropagatorsInPropagators",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "EulerMatrixOfNumeratorsInPropagators",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "PairOfMatricesOfLoopDiagramInIndependentLorentzInvariants",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "PairOfMatricesOfLoopDiagramInPropagators",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "EulerMatrixOfLoopDiagramInPropagators",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "MatrixOfIBPRelations",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "BasisOfIBPRelations",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "MatrixOfSpecialIBPRelations",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "BasisOfSpecialIBPRelations",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "MatrixOfIBPRelationsInWeylAlgebra",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "BasisOfIBPRelationsInWeylAlgebra",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "MatrixOfSpecialIBPRelationsInWeylAlgebra",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "BasisOfSpecialIBPRelationsInWeylAlgebra",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "SymanzikPolynomials",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "AssociatedWeylAlgebra",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "GeneratorsOfScalelessSectors",
        IsLoopDiagram );

#! @Arguments LD
DeclareAttribute( "GeneratorsOfScalelessSectorsInWeylAlgebra",
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
DeclareOperation( "ExpressInExtraLorentzInvariants",
        [ IsHomalgMatrix, IsLoopDiagram ] );

#!
DeclareOperation( "ExpressInIndependentLorentzInvariants",
        [ IsHomalgMatrix, IsLoopDiagram ] );

#!
DeclareOperation( "ExpressInPropagatorsAndNumeratorsAndExtraLorentzInvariants",
        [ IsHomalgMatrix, IsLoopDiagram ] );

#! @Arguments vec, LD
DeclareOperation( "JacobianOfCoefficientsVectorInPropagators",
        [ IsHomalgMatrix, IsLoopDiagram ] );

#! @Arguments vec, LD
DeclareOperation( "DivergenceOfCoefficientsVectorOfLoopDiagram",
        [ IsHomalgMatrix, IsLoopDiagram ] );

#! @Arguments R
DeclareOperation( "DoubleShiftAlgebra",
        [ IsHomalgRing ] );

#! @Arguments R
DeclareOperation( "RationalDoubleShiftAlgebra",
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

#! @Arguments LD, exponents
DeclareOperation( "MatrixOfIBPRelations",
        [ IsLoopDiagram, IsList ] );

#! @Arguments LD, exponents
DeclareOperation( "MatrixOfSpecialIBPRelations",
        [ IsLoopDiagram, IsList ] );

#! @Arguments vec, LD
DeclareOperation( "IBPRelationInWeylAlgebra",
        [ IsHomalgMatrix, IsLoopDiagram ] );

#! @Arguments mat, LD
DeclareOperation( "MatrixOfIBPRelationsInWeylAlgebra",
        [ IsHomalgMatrix, IsLoopDiagram ] );

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

#! @Arguments operator, G
DeclareOperation( "NormalForm",
        [ IsHomalgRingElement, IsHomalgMatrix ] );

#! @Arguments operator, G, initial_integral
DeclareOperation( "NormalForm",
        [ IsHomalgRingElement, IsHomalgMatrix, IsList ] );

#! @Arguments operator, G
DeclareOperation( "NormalFormWrtInitialIntegral",
        [ IsHomalgRingElement, IsHomalgMatrix ] );

####################################
#
#! @Section Tools
#
####################################

