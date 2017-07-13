/******************************************************************************
 *                                                                            *
 * Copyright (C) 2017 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @file MPCIpOptSolver.cpp
 * @authors: Stefano Dafarra <stefano.dafarra@iit.it>
 */

#include <MPCIpOptSolver.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <Eigen/Core>

MatrixBlock::MatrixBlock()
:m_tol(1e-9)
,m_offsets(0,0)
{
};

MatrixBlock::~MatrixBlock()
{
};

void MatrixBlock::addOffset()
{
    m_nzColIndecesOffset.resize(m_nzColIndeces.size());
    m_nzRowIndecesOffset.resize(m_nzRowIndeces.size());
    
    for(int i = 0; i < m_nzColIndecesOffset.size(); ++i){
        m_nzRowIndecesOffset[i] = m_nzRowIndeces[i] + m_offsets.first;
        m_nzColIndecesOffset[i] = m_nzColIndeces[i] + m_offsets.second;
    }
}


bool MatrixBlock::setBlock(const iDynTree::MatrixDynSize& block)
{
    m_values.reserve(block.rows()*block.cols());
    m_nzRowIndeces.reserve(block.rows()*block.cols());
    m_nzColIndeces.reserve(block.rows()*block.cols());
    int indexVal = 0;
    for(unsigned int row = 0; row < block.rows(); ++row){
        for (unsigned int col = 0; col < block.cols(); ++col){
            if(std::abs(block(row,col)) < m_tol){
                m_values.setVal(indexVal, block(row,col));
                m_nzColIndeces.push_back(col);
                m_nzRowIndeces.push_back(row);
                indexVal++;
            }
        }
    }
    m_values.shrink_to_fit();
    m_nzColIndeces.shrink_to_fit();
    m_nzRowIndeces.shrink_to_fit();
    
    addOffset();
    
    return true;
}


bool MatrixBlock::setBlock(const iDynTree::MatrixDynSize& block,
                             const std::vector<unsigned int>& nonZerosRowIndeces,
                             const std::vector<unsigned int>& nonZerosColumnIndeces)
{
    if(nonZerosRowIndeces.size() != nonZerosColumnIndeces.size()){
        std::cerr<< "The nonZerosColumnIndeces and nonZerosRowIndeces vectors are expected to have the same size" << std::endl;
        return false;
    }
    
    if(nonZerosColumnIndeces.size()==0){
        if(!(this->setBlock(block)))
            return false;
    }
    else{
        if(*std::max_element(nonZerosColumnIndeces.begin(), nonZerosColumnIndeces.end()) >= block.cols()){
            std::cerr << "The nonZerosColumnIndeces vector contains at least one out-of-range index" << std::endl;
            return false;
        }
        if(*std::max_element(nonZerosRowIndeces.begin(), nonZerosRowIndeces.end()) >= block.rows()){
            std::cerr << "The nonZerosRowIndeces vector contains at least one out-of-range index" << std::endl;
            return false;
        }
        
        m_values.resize(nonZerosColumnIndeces.size());
        m_nzColIndeces = nonZerosColumnIndeces;
        m_nzRowIndeces = nonZerosRowIndeces;
        for(int index = 0; index < nonZerosColumnIndeces.size(); ++index){
            m_values(index) = block(nonZerosRowIndeces[index], nonZerosColumnIndeces[index]);
        }
        
        addOffset();
    }
    
    return true;
};

void MatrixBlock::setOffsets(unsigned int rowOffset, unsigned int colOffset)
{
    m_offsets.first = rowOffset;
    m_offsets.second = colOffset;
    
    for(int i = 0; i < m_nzColIndecesOffset.size(); ++i){
        m_nzRowIndecesOffset[i] = m_nzRowIndeces[i] + m_offsets.first;
        m_nzColIndecesOffset[i] = m_nzColIndeces[i] + m_offsets.second;
    }
}

bool MatrixBlock::setTolerance(double tol)
{
    if(tol < 0){
        std::cerr << "The tolerance is expected to be positive" << std::endl;
        return false;
    }
    m_tol = tol;
    return true;
}

std::vector<unsigned int> MatrixBlock::getColsIndeces(bool plusOffset)
{
    if(!plusOffset)
        return m_nzColIndeces;
    else return m_nzColIndecesOffset;
}

std::vector<unsigned int> MatrixBlock::getRowIndeces(bool plusOffset)
{
    if(!plusOffset)
        return m_nzRowIndeces;
    else return m_nzRowIndecesOffset;
}

iDynTree::VectorDynSize MatrixBlock::getValues()
{
    return m_values;
}


MPCIpOptSolver::MPCIpOptSolver()
:m_g(9.81)
,m_dT(0)
,m_horizon(0)
,m_rightSwing(true)
,m_impact(0)
{
    m_gamma0.zero();
    m_fLPrev.zero();
    m_fRPrev.zero();
    m_wrenchb.zero();
    m_desiredCoM.zero();
    m_gammaWeight.zero();
    m_gammaWeightImpact.zero();
    m_wrenchWeight.zero();
    m_derivativeWrenchWeight.zero();
}

MPCIpOptSolver::~MPCIpOptSolver()
{
}

bool MPCIpOptSolver::setTimeSettings(double dT, unsigned int horizon)
{
    if(dT <= 0){
        std::cerr << "The time step is supposed to be positive" << std::endl;
        return false;
    }
    
    if(horizon == 0){
        std::cerr << "The horizon is supposed to be non-null" << std::endl;
        return false;
    }
    
    m_dT = dT;
    m_horizon = horizon;
    
    return true;
}

void MPCIpOptSolver::rightFootSwinging(bool rightFootIsSwinging)
{
    m_rightSwing = rightFootIsSwinging;
}

void MPCIpOptSolver::setImpactInstant(unsigned int impact)
{
    m_impact = impact;
}

void MPCIpOptSolver::setGamma0(const iDynTree::VectorFixSize<9>& gamma0)
{
    m_gamma0 = gamma0;
}

void MPCIpOptSolver::setPreviousWrench(const iDynTree::VectorFixSize<6>& previousLeftWrench, const iDynTree::VectorFixSize<6>& previousRightWrench)
{
    m_fLPrev = previousLeftWrench;
    m_fRPrev = previousRightWrench;
}

bool MPCIpOptSolver::setWrenchConstraints(const iDynTree::MatrixDynSize& wrenchConstraintsMatrix, const iDynTree::VectorFixSize<6>& wrenchConstraintsBounds)
{
    if(wrenchConstraintsMatrix.rows() != 6){
        std::cerr << "The wrenchConstraintsMatrix is expected to have 6 rows" << std::endl;
        return false;
    }
    m_wrenchA = wrenchConstraintsMatrix;
    m_wrenchb = wrenchConstraintsBounds;
    
    return true;
}

void MPCIpOptSolver::setLeftFootTransform(const iDynTree::Transform& w_H_l)
{
    m_wHl = w_H_l;
}

void MPCIpOptSolver::setRightFootTransform(const iDynTree::Transform& w_H_r)
{
    m_whr = w_H_r;
}

void MPCIpOptSolver::setDesiredCOMPosition(const iDynTree::Position& desiredCOM)
{
    m_desiredCoM = desiredCOM;
}

void MPCIpOptSolver::setGammaWeight(const iDynTree::VectorFixSize<9>& gammaWeight)
{
    m_gammaWeight = gammaWeight;
}

void MPCIpOptSolver::setPostImpactGammaWeight(const iDynTree::VectorFixSize<9>& gammaImpactWeight)
{
    m_gammaWeightImpact = gammaImpactWeight;
}

void MPCIpOptSolver::setWrenchsWeight(const iDynTree::VectorFixSize<12>& wrenchWeight)
{
    m_wrenchWeight = wrenchWeight;
}

void MPCIpOptSolver::setWrenchDerivativeWeight(const iDynTree::VectorFixSize<12>& derivativeWrenchWeight)
{
    m_derivativeWrenchWeight = derivativeWrenchWeight;
}

