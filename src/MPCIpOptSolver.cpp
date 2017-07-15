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

MPCIpOptSolver::MPCIpOptSolver()
:m_g(9.81)
,m_mass(0.0)
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
    m_EvGamma.resize(9,9);
    m_EvGamma.zero();
    m_EvGammaSparse.resize(9,9);
    m_FGamma.resize(9,12);
    m_FGamma.zero();
    m_FGammaSparse.resize(9,12);
    m_skewBuffer.resize(3,3);
    m_bias.zero();
    m_minusIdentity.resize(9,9);
    
    iDynTree::Triplets values;
    values.addDiagonalMatrix(0,0,-1,9);
    m_minusIdentity.setFromTriplets(values);
    
    m_wrenchTransform.resize(6,6);
    m_wrenchTransform.zero();
    
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
    m_modelConstraintsJacobian.resize(3*m_horizon);
    
    return true;
}

bool MPCIpOptSolver::setRobotMass(const double mass)
{
    if(mass <= 0){
        std::cerr << "The mass is expected to be a positive number." << std::endl;
        return false;
    }
    m_mass = mass;
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
    
    if(wrenchConstraintsMatrix.rows() != wrenchConstraintsBounds.size()){
        std::cerr << "Unbalanced dimensions between the matrix and the vector of constraints." << std::endl;
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
    m_wHr = w_H_r;
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

bool MPCIpOptSolver::computeModelMatrices()
{
    if(m_mass==0){
        std::cerr << "First you have to define the mass of the robot." << std::endl;
        return false;
    }
    
    if((m_dT == 0)||(m_horizon == 0)){
        std::cerr<< "First you have to specify the time settings." << std::endl;
        return false;
    }
    
    //EV_gamma
    
    iDynTree::iDynTreeEigenMatrixMap map_EvGamma = iDynTree::toEigen(m_EvGamma);
    Eigen::MatrixXd identity9;
    identity9.resize(9,9);
    identity9.setIdentity();
    iDynTree::Triplets valuesEV;
    valuesEV.reserve(21);
    
    map_EvGamma = identity9;
    valuesEV.addDiagonalMatrix(0, 0, 1, 9);
    
    map_EvGamma.block<3,3>(0,3) = m_dT*identity9.block<3,3>(0,0);
    valuesEV.addDiagonalMatrix(0, 3, m_dT, 3);

    Eigen::Map <Eigen::VectorXd> fl_map(m_fLPrev.data(), 6);
    Eigen::Map <Eigen::VectorXd> fr_map (m_fRPrev.data(), 6);
    Eigen::Vector3d temp = fl_map.head(3) + fr_map.head(3);
    iDynTree::toEigen(m_skewBuffer) = m_dT * iDynTree::skew(temp);
    map_EvGamma.block<3,3>(6,0) = iDynTree::toEigen(m_skewBuffer);
    valuesEV.addSubMatrix(6,0,m_skewBuffer);
    
    m_EvGammaSparse.setFromTriplets(valuesEV);
    
    //F_gamma
    iDynTree::Triplets valuesF;
    double tempRatio = m_dT/m_mass;
    valuesF.reserve(30);
    iDynTree::iDynTreeEigenMatrixMap map_FGamma = iDynTree::toEigen(m_FGamma);
    
    map_FGamma.block<3,3>(3,0) = tempRatio*identity9.block<3,3>(0,0);
    valuesF.addDiagonalMatrix(3,0,tempRatio,3);
    
    map_FGamma.block<3,3>(3,6) = tempRatio*identity9.block<3,3>(0,0);
    valuesF.addDiagonalMatrix(3,6,tempRatio,3);
    
    map_FGamma.block<3,3>(6,3) = m_dT*identity9.block<3,3>(0,0);
    valuesF.addDiagonalMatrix(6, 3, m_dT, 3);
    
    map_FGamma.block<3,3>(6,9) = m_dT*identity9.block<3,3>(0,0);
    valuesF.addDiagonalMatrix(6, 9, m_dT, 3);
    
    Eigen::Map <Eigen::VectorXd> gamma0_map(m_gamma0.data(), 9);
    Eigen::Map <const Eigen::VectorXd> xl_map(m_wHl.getPosition().data(), 3);
    Eigen::Map <const Eigen::VectorXd> xr_map(m_wHr.getPosition().data(), 3);
    
    temp = xl_map - gamma0_map.head<3>();
    iDynTree::toEigen(m_skewBuffer) = m_dT * iDynTree::skew(temp);
    map_FGamma.block<3,3>(6,0) = iDynTree::toEigen(m_skewBuffer);
    valuesF.addSubMatrix(6,0,m_skewBuffer);
    
    temp = xr_map - gamma0_map.head<3>();
    iDynTree::toEigen(m_skewBuffer) = m_dT * iDynTree::skew(temp);
    map_FGamma.block<3,3>(6,6) = iDynTree::toEigen(m_skewBuffer);
    valuesF.addSubMatrix(6,6,m_skewBuffer);
    
    m_FGammaSparse.setFromTriplets(valuesF);
    
    return true;
}

bool MPCIpOptSolver::computeModelBias()
{
    Eigen::Map <Eigen::VectorXd> bias_map (m_bias.data(), 9);
    
    bias_map[5] = -m_g*m_dT;
    Eigen::Map <Eigen::VectorXd> fl_map(m_fLPrev.data(), 6);
    Eigen::Map <Eigen::VectorXd> fr_map (m_fRPrev.data(), 6);
    Eigen::Map <Eigen::VectorXd> gamma0_map(m_gamma0.data(), 9);
    
    Eigen::Vector3d temp = fl_map.head<3>() + fr_map.head<3>();
    iDynTree::toEigen(m_skewBuffer) = m_dT * iDynTree::skew(temp);
    bias_map.tail<3>() = -m_dT*iDynTree::toEigen(m_skewBuffer)*gamma0_map.head(3);
    
    return true;
}

bool MPCIpOptSolver::computeModelConstraintsJacobian()
{
    MatrixBlock templateEv, templateF, templateId;
    
    templateEv.blockPtr.reset(&m_EvGammaSparse);
    templateF.blockPtr.reset(&m_FGammaSparse);
    templateId.blockPtr.reset(&m_minusIdentity);
    
    m_modelConstraintsJacobian.resize(3*m_horizon);
    
    unsigned int row = 0;
    unsigned int col = 0;
    
    templateId.rowOffset = row;
    templateId.colOffset = col;
    m_modelConstraintsJacobian[0] = templateId;
    col += 9;
    
    templateF.rowOffset = row;
    templateF.colOffset = col;
    m_modelConstraintsJacobian[1] = templateF;
    row += 9;
    col = 0;
    
    int index = 2;
    for(int t = 1; t < m_horizon; ++t){
        
        templateEv.rowOffset = row;
        templateEv.colOffset = col;
        m_modelConstraintsJacobian[index] = templateEv;
        ++index;
        col += 21;
        
        templateId.rowOffset = row;
        templateId.colOffset = col;
        m_modelConstraintsJacobian[index] = templateId;
        ++index;
        col += 9;
        
        templateF.rowOffset = row;
        templateF.colOffset = col;
        m_modelConstraintsJacobian[index] = templateF;
        ++index;
        
        row += 9;
        col = 21*(t-1);
    }

    return true;
}

bool MPCIpOptSolver::computeWrenchConstraints()
{
    iDynTree::iDynTreeEigenMatrixMap Al_map = iDynTree::toEigen(m_wrenchAl);
    iDynTree::iDynTreeEigenMatrixMap Ar_map = iDynTree::toEigen(m_wrenchAr);
    iDynTree::iDynTreeEigenMatrixMap A_map = iDynTree::toEigen(m_wrenchA);
    iDynTree::iDynTreeEigenMatrixMap wrenchTransform_map = iDynTree::toEigen(m_wrenchTransform);
    iDynTree::iDynTreeEigenConstMatrixMap rotationL_map (m_wHl.getRotation().data(), 3, 3);
    iDynTree::iDynTreeEigenConstMatrixMap rotationR_map (m_wHr.getRotation().data(), 3, 3);
    
    wrenchTransform_map.block<3,3>(0,0) = rotationL_map.transpose();
    wrenchTransform_map.block<3,3>(3,3) = rotationL_map.transpose();
    
    Al_map = A_map*wrenchTransform_map;
    
    wrenchTransform_map.block<3,3>(0,0) = rotationR_map.transpose();
    wrenchTransform_map.block<3,3>(3,3) = rotationR_map.transpose();
    
    Ar_map = A_map*wrenchTransform_map;
    
    return true;
}

