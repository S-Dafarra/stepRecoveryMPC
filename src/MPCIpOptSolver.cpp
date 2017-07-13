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
