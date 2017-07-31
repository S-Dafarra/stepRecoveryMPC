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
,m_mass(30.0)
,m_dT(0)
,m_horizon(0)
,m_impact(0)
,m_exitCode(-6)
,m_rightFootStep(true)
{
    m_gamma0.zero();
    m_fLPrev.resize(6);
    m_fLPrev.zero();
    m_fRPrev.resize(6);
    m_fRPrev.zero();
    m_desiredGamma.zero();
    m_gammaWeight.resize(9);
    m_gammaWeight.zero();
    m_gammaWeightImpact.resize(9);
    m_gammaWeightImpact.zero();
    m_wrenchWeight.resize(12);
    m_wrenchWeight.zero();
    m_derivativeWrenchWeight.resize(12);
    m_derivativeWrenchWeight.zero();
    m_EvGamma.resize(9,9);
    m_EvGamma.zero();
    m_EvGammaSparsePtr = std::make_shared<iDynTree::SparseMatrix>();
    m_EvGammaSparsePtr->resize(9,9);
    m_FGamma.resize(9,12);
    m_FGamma.zero();
    m_FGammaSparsePtr = std::make_shared<iDynTree::SparseMatrix>();
    m_FGammaSparsePtr->resize(9,12);
    m_skewBuffer.resize(3,3);
    m_bias.zero();
    m_minusIdentityPtr = std::make_shared<iDynTree::SparseMatrix>();
    m_minusIdentityPtr->resize(9,9);
    m_wrenchAlSparsePtr = std::make_shared<iDynTree::SparseMatrix>();
    m_wrenchArSparsePtr = std::make_shared<iDynTree::SparseMatrix>();
    iDynTree::Triplets values;
    values.addDiagonalMatrix(0,0,-1,9);
    m_minusIdentityPtr->setFromTriplets(values);
    
    m_wrenchTransform.resize(6,6);
    m_wrenchTransform.zero();
    
    iDynTree::Position dummy;
    dummy.zero();
    m_wHl.setRotation(iDynTree::Rotation::Identity());
    m_wHl.setPosition(dummy);
    m_wHr.setRotation(iDynTree::Rotation::Identity());
    m_wHr.setPosition(dummy);
    
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
    m_modelConstraintsJacobian.resize(3*m_horizon - 1);
    m_costHessian.resize(21*m_horizon, 21*m_horizon);
    
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

bool MPCIpOptSolver::setImpactInstant(unsigned int impact, bool rightFootStepping)
{
    m_impact = impact;
    m_rightFootStep = rightFootStepping;
    return true;
}

bool MPCIpOptSolver::setGamma0(const iDynTree::VectorFixSize<9>& gamma0)
{
    m_gamma0 = gamma0;
    //std::cerr << "Gamma0: " << m_gamma0.toString() << std::endl;
    return true;
}

bool MPCIpOptSolver::setPreviousWrench(const iDynTree::VectorDynSize& previousLeftWrench, const iDynTree::VectorDynSize& previousRightWrench)
{
    if(previousLeftWrench.size() != 6){
        std::cerr << "The previous left wrench is expected to be 6 dimensional" << std::endl;
        return false;
    }
    if(previousRightWrench.size() != 6){
        std::cerr << "The previous right wrench is expected to be 6 dimensional" << std::endl;
        return false;
    }
    m_fLPrev = previousLeftWrench;
    m_fRPrev = previousRightWrench;
    
    return true;
}


bool MPCIpOptSolver::setWrenchConstraints(const iDynTree::MatrixDynSize& wrenchConstraintsMatrix, const iDynTree::VectorDynSize& constraintsBounds, const iDynTree::VectorDynSize& afterImpactConstraintsBounds)
{
    if(wrenchConstraintsMatrix.cols() != 6){
        std::cerr << "The wrenchConstraintsMatrix is expected to have 6 columns" << std::endl;
        return false;
    }
    
    if(wrenchConstraintsMatrix.rows() != constraintsBounds.size()){
        std::cerr << "Unconsistent dimension between wrenchConstraintsMatrix and constraintsBounds." << std::endl;
        return false;
    }
    
    if(afterImpactConstraintsBounds.size() != constraintsBounds.size()){
        std::cerr << "The dimension of constraints should not change across the impact." << std::endl;
        return false;
    }
    
    m_wrenchA = wrenchConstraintsMatrix;
    m_wrenchAl = wrenchConstraintsMatrix;
    m_wrenchAr = wrenchConstraintsMatrix;
    m_wrenchb  = constraintsBounds;
    m_wrenchbImpact = afterImpactConstraintsBounds;
    
//     std::cerr << "Constr: " << std::endl << wrenchConstraintsMatrix.toString() << std::endl;
//     std::cerr << "b: " << std::endl << constraintsBounds.toString() << std::endl;
//     std::cerr << "bImpact: " << std::endl << afterImpactConstraintsBounds.toString() << std::endl;
    
    
    return true;
}


bool MPCIpOptSolver::setFeetTransforms(const iDynTree::Transform& w_H_l, const iDynTree::Transform& w_H_r)
{
    m_wHl = w_H_l;
    m_wHr = w_H_r;
    
//     std::cerr << "Left Foot Transform: " << std::endl << m_wHl.toString() << std::endl;
//     std::cerr << "Right Foot Transform: " << std::endl << m_wHr.toString() << std::endl;
    
    if(m_wrenchAl.rows() != 0){
        if(!computeWrenchConstraints()){
            std::cerr << "Error while computing the wrench constraints." << std::endl;
            return false;
        }
    }
    return true;
}


bool MPCIpOptSolver::setDesiredCOMPosition(const iDynTree::Position& desiredCOM)
{
    Eigen::Map<const Eigen::VectorXd> desiredCOM_map(desiredCOM.data(), 3);
    iDynTree::toEigen(m_desiredGamma).head<3>() = desiredCOM_map;
    
    std::cerr << "desiredGamma: "<< m_desiredGamma.toString() << std::endl;
    
    return true;
}

bool MPCIpOptSolver::setGammaWeight(const iDynTree::VectorDynSize& gammaWeight)
{
    if(gammaWeight.size() != 9){
        std::cerr << "The gammaWeight vector is expected to be 9 dimensional" << std::endl;
        return false;
    }
    m_gammaWeight = gammaWeight;
    
  //  std::cerr << "gammaWeight: "<< m_gammaWeight.toString() << std::endl;

    return true;
}

bool MPCIpOptSolver::setPostImpactGammaWeight(const iDynTree::VectorDynSize& gammaImpactWeight)
{
    if(gammaImpactWeight.size() != 9){
        std::cerr << "The gammaWeight vector is expected to be 9 dimensional" << std::endl;
        return false;
    }
    
    m_gammaWeightImpact = gammaImpactWeight;
    
   // std::cerr << "gammaWeightImpact: "<< m_gammaWeightImpact.toString() << std::endl;

    return true;
}

bool MPCIpOptSolver::setWrenchsWeight(const iDynTree::VectorDynSize& wrenchWeight)
{
    if(wrenchWeight.size() != 12){
        std::cerr << "The wrenchWeight vector is expected to be 12 dimensional" << std::endl;
        return false;
    }
    m_wrenchWeight = wrenchWeight;
    
   // std::cerr << "wrenchWeight: "<< m_wrenchWeight.toString() << std::endl;

    return true;
}

bool MPCIpOptSolver::setWrenchDerivativeWeight(const iDynTree::VectorDynSize& derivativeWrenchWeight)
{
    if(derivativeWrenchWeight.size() != 12){
        std::cerr << "The derivativeWrenchWeight vector is expected to be 12 dimensional" << std::endl;
        return false;
    }
    m_derivativeWrenchWeight = derivativeWrenchWeight;
    
    //std::cerr << "wrenchWeight: "<< derivativeWrenchWeight.toString() << std::endl;
    
    
    return true;
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

    Eigen::Map <Eigen::VectorXd> fl_map(m_fLPrev.data(), m_fLPrev.size());
    Eigen::Map <Eigen::VectorXd> fr_map (m_fRPrev.data(), m_fRPrev.size());
    Eigen::Vector3d temp = fl_map.head(3) + fr_map.head(3);
    iDynTree::toEigen(m_skewBuffer) = m_dT * iDynTree::skew(temp);
    map_EvGamma.block<3,3>(6,0) = iDynTree::toEigen(m_skewBuffer);
    valuesEV.addSubMatrix(6,0,m_skewBuffer);
    
    m_EvGammaSparsePtr->setFromTriplets(valuesEV);
    
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
    
    m_FGammaSparsePtr->setFromTriplets(valuesF);
    
//     std::cerr << "Ev: "<< std::endl << m_EvGamma.toString() << std::endl;
//    std::cerr << "F: " << std::endl << m_FGamma.toString() << std::endl;
//    std::cerr << "FSparse: " << std::endl << m_FGammaSparsePtr->description(true) << std::endl;
    
    
    return true;
}

bool MPCIpOptSolver::computeModelBias()
{
    Eigen::Map <Eigen::VectorXd> bias_map (m_bias.data(), 9);
    
    bias_map[5] = -m_g*m_dT;
    Eigen::Map <Eigen::VectorXd> fl_map (m_fLPrev.data(), m_fLPrev.size());
    Eigen::Map <Eigen::VectorXd> fr_map (m_fRPrev.data(), m_fRPrev.size());
    Eigen::Map <Eigen::VectorXd> gamma0_map(m_gamma0.data(), 9);
    
    Eigen::Vector3d temp = fl_map.head<3>() + fr_map.head<3>();
    iDynTree::toEigen(m_skewBuffer) = m_dT * iDynTree::skew(temp);
    bias_map.tail<3>() = -m_dT*iDynTree::toEigen(m_skewBuffer)*gamma0_map.head(3);
    
//     std::cerr << "Bias: "<< std::endl << m_bias.toString() << std::endl;
    
    return true;
}

bool MPCIpOptSolver::computeModelConstraintsJacobian()
{
    MatrixBlock templateEv, templateF, templateId;
    
    templateEv.blockPtr = m_EvGammaSparsePtr;
    templateF.blockPtr = m_FGammaSparsePtr;
    templateId.blockPtr = m_minusIdentityPtr;
    
    m_modelConstraintsJacobian.resize(3*m_horizon - 1);
    
    unsigned int row = 0;
    unsigned int col = 0;
    
    templateId.rowOffset = row;
    templateId.colOffset = col;
    m_modelConstraintsJacobian[0] = templateId;
    col += 9;
    
    templateF.rowOffset = row;
    templateF.colOffset = col;
    m_modelConstraintsJacobian[1] = templateF;
    row = 9;
    
    int index = 2;
    for(int t = 1; t < m_horizon; ++t){
        col = 21*(t-1);
        
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
    }

    return true;
}

bool MPCIpOptSolver::computeWrenchConstraints()
{
    if(m_wrenchAl.rows() == 0){
        std::cerr << "First you have to load the wrenchConstraintsMatrix." << std::endl;
        return false;
    }
    
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
    
    m_wrenchAlSparsePtr->resize(m_wrenchAl.rows(), m_wrenchAl.cols());
    m_wrenchArSparsePtr->resize(m_wrenchAr.rows(), m_wrenchAr.cols());
    
    iDynTree::Triplets valuesL, valuesR;
    valuesL.setSubMatrix(0,0,m_wrenchAl); //maybe we can ask for the sparsity of this matrix
    m_wrenchAlSparsePtr->setFromTriplets(valuesL);
    
    valuesR.setSubMatrix(0,0,m_wrenchAr);
    m_wrenchArSparsePtr->setFromTriplets(valuesR);
    
//    std::cerr << "Al:" << std::endl << m_wrenchAl.toString() << std::endl;
//    std::cerr << "Ar:" << std::endl << m_wrenchAr.toString() << std::endl;
    
    return true;
}

bool MPCIpOptSolver::computeWrenchConstraintsJacobian()
{
    MatrixBlock templateLeft, templateRight;
    
    templateLeft.blockPtr = m_wrenchAlSparsePtr;
    templateRight.blockPtr =  m_wrenchArSparsePtr;
    
    m_wrenchConstraintJacobian.resize(m_horizon*2);
    
    unsigned int row = 0;
    unsigned int col = 9;
    unsigned int index = 0;
    
    unsigned int nConstraintsL = m_wrenchAl.rows();
    unsigned int nConstraintsR = m_wrenchAr.rows();
    
    for (int t=0; t<m_horizon; ++t){
        col = 9 + 21*t; //you have to pick the right variables inside the decision variables vector
        templateLeft.rowOffset = row;
        templateLeft.colOffset = col;
        m_wrenchConstraintJacobian[index] = templateLeft;
        index++;
        col += 6;
        row += nConstraintsL;
        
        templateRight.rowOffset = row;
        templateRight.colOffset = col;
        m_wrenchConstraintJacobian[index] = templateRight;
        index++; 
        
        row += nConstraintsR;
    }
    
    return true;
}

bool MPCIpOptSolver::computeSingleCostHessian()
{
    m_gammaWeightHessian.resize(9,9);
    m_gammaWeightImpactHessian.resize(9,9);
    m_wrenchWeightHessian.resize(12,12);
    m_derivativeWrenchWeightHessian.resize(12,12);
    m_negativeDerWrenchHessian.resize(12,12);
    
    iDynTree::Triplets valuesGamma, valuesGammaImpact;
    for(int i=0; i < m_gammaWeight.size(); ++i){
        valuesGamma.addDiagonalMatrix(i, i, m_gammaWeight(i), 1);
        valuesGammaImpact.addDiagonalMatrix(i, i, m_gammaWeightImpact(i), 1);
    }
    m_gammaWeightHessian.setFromTriplets(valuesGamma);
    m_gammaWeightImpactHessian.setFromTriplets(valuesGammaImpact);
    
    iDynTree::Triplets valuesWrench, valuesWrenchDerivative, valuesNegDerivative;
    for(int i=0; i < m_wrenchWeight.size(); ++i){
        valuesWrench.addDiagonalMatrix(i, i, m_wrenchWeight(i),1);
        valuesWrenchDerivative.addDiagonalMatrix(i, i, m_derivativeWrenchWeight(i), 1);
        valuesNegDerivative.addDiagonalMatrix(i, i, -m_derivativeWrenchWeight(i), 1);
    }
    m_wrenchWeightHessian.setFromTriplets(valuesWrench);
    m_derivativeWrenchWeightHessian.setFromTriplets(valuesWrenchDerivative);
    m_negativeDerWrenchHessian.setFromTriplets(valuesNegDerivative);
    
    return true;
}

bool MPCIpOptSolver::computeCostHessian()
{
    if(!computeSingleCostHessian())
        return false;
    
    iDynTree::Triplets values;
    
    for(int t=0; t<m_horizon; ++t){
        values.addSubMatrix(21*t, 21*t, m_gammaWeightHessian);
        if((t >= m_impact)&&(t < (m_horizon-1))){
            values.addSubMatrix(21*t, 21*t, m_gammaWeightImpactHessian);
        }
        values.addSubMatrix(9 + 21*t, 9 + 21*t, m_wrenchWeightHessian);
        if(t == 0){
            values.addSubMatrix(9, 9, m_derivativeWrenchWeightHessian);
        }
        else {
            values.addSubMatrix(9 + 21*t, 9+21*t, m_derivativeWrenchWeightHessian);
            values.addSubMatrix(9 + 21*t, 9+21*(t-1), m_negativeDerWrenchHessian);
            values.addSubMatrix(9 + 21*(t-1), 9 + 21*(t-1), m_derivativeWrenchWeightHessian);
            values.addSubMatrix(9 + 21*(t-1), 9 + 21*t, m_negativeDerWrenchHessian);
        }
    }
    values.addSubMatrix(21*(m_horizon-1), 21*(m_horizon-1), m_gammaWeightImpactHessian); //terminal cost
    m_costHessian.setFromTriplets(values);
    
    //std::cerr << "Cost hessian: " << std::endl << m_costHessian.description(true) << std::endl;
    
    return true;
}

bool MPCIpOptSolver::prepareProblem()
{
    if(!computeModelMatrices()){
        std::cerr << "Error while computing the model matrices." << std::endl;
        return false;
    }
    if(!computeModelBias()){
        std::cerr << "Error while computing the model bias." << std::endl;
        return false;
    }
    if(!computeModelConstraintsJacobian()){
        std::cerr << "Error while computing the model constraints jacobian." << std::endl;
        return false;
    }
    if(!computeWrenchConstraints()){
        std::cerr << "Error while computing the wrench constraints." << std::endl;
        return false;
    }
    if(!computeWrenchConstraintsJacobian()){
        std::cerr << "Error while computing the wrench constraints jacobian." << std::endl;
        return false;
    }
    if(!computeCostHessian()){
        std::cerr << "Error while computing the cost hessian." << std::endl;
        return false;
    }

    // initialize buffers for multipliers
    size_t n = 21*m_horizon;
    size_t m = (9+m_wrenchAl.rows()+m_wrenchAr.rows())*m_horizon;
    m_lowerBoundMultipliers.resize(n);
    m_lowerBoundMultipliers.zero();
    m_upperBoundMultipliers.resize(n);
    m_upperBoundMultipliers.zero();
    m_contraintMultipliers.resize(m);
    m_contraintMultipliers.zero();
    m_previousSolution.resize(n);
    m_previousSolution.zero();
    
    return true;
}


bool MPCIpOptSolver::updateProblem()
{
    if(!computeModelMatrices()){
        std::cerr << "Error while computing the model matrices." << std::endl;
        return false;
    }
    if(!computeModelBias()){
        std::cerr << "Error while computing the model bias." << std::endl;
        return false;
    }
    if(!computeWrenchConstraints()){
        std::cerr << "Error while computing the wrench constraints." << std::endl;
        return false;
    }
    
 //   printJacobian("Model jacobian:", m_modelConstraintsJacobian, 9*m_horizon, 21*m_horizon);
//    printJacobian("Wrench jacobian:", m_wrenchConstraintJacobian, (m_wrenchAl.rows()+m_wrenchAr.rows())*m_horizon, 21*m_horizon);
    
    
    return true;
}

int MPCIpOptSolver::getSolution(iDynTree::VectorDynSize& fL, iDynTree::VectorDynSize& fR, iDynTree::VectorDynSize& lastGamma)
{
    if(m_previousSolution.size() >= 21){
        Eigen::Map<Eigen::VectorXd> solution_map(m_previousSolution.data(),m_previousSolution.size());
        fL.resize(6);
        fR.resize(6);
        lastGamma.resize(9);
        iDynTree::toEigen(fL) = solution_map.segment<6>(9);
        iDynTree::toEigen(fR) = solution_map.segment<6>(9+6);
        iDynTree::toEigen(lastGamma) = solution_map.segment<9>(21*(m_horizon-1));
    }
    
   // std::cerr << "Solution left: " << fL.toString() << std::endl;
    //std::cerr << "Solution right: " << fR.toString() << std::endl;
    //std::cerr << "Solution lastGamma: " << lastGamma.toString() << std::endl;
   // std::cerr << "Full gamma: "<<m_previousSolution.toString() << std::endl;
    
    return m_exitCode;
}


bool MPCIpOptSolver::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style)
{
    n = 21*m_horizon;
    m = (9+m_wrenchAl.rows()+m_wrenchAr.rows())*m_horizon;
    nnz_jac_g = 0;
    for(int i=0; i<m_modelConstraintsJacobian.size(); ++i){
        nnz_jac_g += m_modelConstraintsJacobian[i].blockPtr->numberOfNonZeros();
    }
    for(int i=0; i<m_wrenchConstraintJacobian.size(); ++i){
        nnz_jac_g += m_wrenchConstraintJacobian[i].blockPtr->numberOfNonZeros();
    }
    nnz_h_lag = m_costHessian.numberOfNonZeros();
    
    index_style = Ipopt::TNLP::C_STYLE;
    return true;
}

bool MPCIpOptSolver::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u)
{
    Eigen::Map<Eigen::VectorXd> x_l_map(x_l,n);
    Eigen::Map<Eigen::VectorXd> x_u_map(x_u,n);
    
    x_l_map.setConstant(-2e+19);
    x_u_map.setConstant( 2e+19);
    
    for(int t = 0; t < std::min(m_horizon,m_impact); ++t){
        x_l_map.segment<6>(21*t + 9 + 6).setZero(); //assuming right foot impact
        x_u_map.segment<6>(21*t + 9 + 6).setZero();
    }
    
    for (Ipopt::Index c = 0; c < m; ++c) {
        if(c < (9*m_horizon)){
            g_l[c] = -1e-9*0;  //State constraints are equality
            g_u[c] = 1e-9*0;
        }
        else{
            g_l[c] = -2e+19;
            g_u[c] =  0;
        }
    }
    
    return true;
}

bool MPCIpOptSolver::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda)
{
    if (init_z) {
        Eigen::Map<Eigen::VectorXd> zlMap(z_L, n);
        Eigen::Map<Eigen::VectorXd> zuMap(z_U, n);
        zlMap = iDynTree::toEigen(m_lowerBoundMultipliers);
        zuMap = iDynTree::toEigen(m_upperBoundMultipliers);
    }
    if (init_lambda)
    {
        Eigen::Map<Eigen::VectorXd> lambdaMap(lambda, m);
        lambdaMap = iDynTree::toEigen(m_contraintMultipliers);
    }
    
    if(init_x){
        Eigen::Map<Eigen::VectorXd> x_map(x, n);
        x_map = iDynTree::toEigen(m_previousSolution);


//        iDynTree::iDynTreeEigenMatrixMap ev_map = iDynTree::toEigen(m_EvGamma);
//        iDynTree::iDynTreeEigenMatrixMap f_map = iDynTree::toEigen(m_FGamma);
//        Eigen::Map<Eigen::VectorXd> bias_map(m_bias.data(), 9);
//        Eigen::Map<Eigen::VectorXd> gamma0_map (m_gamma0.data(), m_gamma0.size());

//         x_map.head((m_horizon-1)*21) = prevSol_map.tail((m_horizon-1)*21);
//         x_map.segment<9>((m_horizon-1)*21) = ev_map*prevSol_map.segment<9>((m_horizon-1)*21) + f_map*prevSol_map.tail<12>() + bias_map;
//         x_map.tail<12>() = prevSol_map.tail<12>();
        
        
//         x_map.setZero();
//         for(int t=0; t < m_horizon; ++t){
//             if(t < m_impact){
//                 x_map(21*t + 9 + 2) = 15; //fz
//                 x_map(21*t + 9 + 6 + 2) = 15;  //fz
//             }
//             if(t == 0){
//                 x_map.head<9>() = ev_map*gamma0_map + f_map*x_map.segment<12>(9) + bias_map;
//             }
//             else{
//                 x_map.segment<9>(21*t) = ev_map*x_map.segment<9>(21*(t-1)) + f_map*x_map.segment<12>(9 + 21*t) + bias_map;
//             }
//         }
    }
//    std::cerr << "Init-------\n";
//
//    Eigen::Map<Eigen::VectorXd> zlMap(z_L, n);
//    std::cerr << zlMap.transpose() << std::endl<< std::endl;
//    Eigen::Map<Eigen::VectorXd> zuMap(z_U, n);
//    std::cerr << zuMap.transpose() << std::endl<< std::endl;
//    Eigen::Map<Eigen::VectorXd> lambdaMap(lambda, m);
//    std::cerr << lambdaMap.transpose() << std::endl<< std::endl;
//    Eigen::Map<Eigen::VectorXd> x_map(x, n);
//    std::cerr << x_map.transpose() << std::endl<< std::endl;

    return true;
}

bool MPCIpOptSolver::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value)
{
    Eigen::Map<const Eigen::VectorXd> x_map(x, n);
    Eigen::Map<const Eigen::VectorXd > gammaWeight_map(m_gammaWeight.data(),m_gammaWeight.size());
    Eigen::Map<const Eigen::VectorXd > gammaImpactWeight_map(m_gammaWeightImpact.data(),m_gammaWeightImpact.size());
    Eigen::Map<const Eigen::VectorXd > wrenchWeight_map(m_wrenchWeight.data(),m_wrenchWeight.size());
    Eigen::Map<const Eigen::VectorXd > derivativeWrenchWeight_map(m_derivativeWrenchWeight.data(),m_derivativeWrenchWeight.size());
    
    double gammaCost = 0;
    double wrenchCost = 0;
    double impactCost = 0;
    double wrenchDiffCost = 0;
    Eigen::Map <Eigen::VectorXd> fl_map(m_fLPrev.data(), m_fLPrev.size());
    Eigen::Map <Eigen::VectorXd> fr_map(m_fRPrev.data(), m_fRPrev.size());
        
    for(int t=0; t < m_horizon; ++t){
        gammaCost += 0.5*(x_map.segment<9>(21*t) - iDynTree::toEigen(m_desiredGamma)).transpose() * gammaWeight_map.asDiagonal() * (x_map.segment<9>(21*t) - iDynTree::toEigen(m_desiredGamma));
        
        if((t >= m_impact)&&(t < (m_horizon-1))){
            impactCost += 0.5*(x_map.segment<9>(21*t) - iDynTree::toEigen(m_desiredGamma)).transpose() * gammaImpactWeight_map.asDiagonal() * (x_map.segment<9>(21*t) - iDynTree::toEigen(m_desiredGamma));
        }
//         std::cerr << "Gamma: " << x_map.segment<9>(21*t).transpose() << std::endl;
//         std::cerr << "Gamma Des: " << m_desiredGamma.toString() << std::endl;
//         std::cerr << "Gamma error: "<< x_map.segment<9>(21*t).transpose() - iDynTree::toEigen(m_desiredGamma).transpose() << std::endl;
        
        wrenchCost += 0.5*x_map.segment<12>(9 + 21*t).transpose() * wrenchWeight_map.asDiagonal() * x_map.segment<12>(9 + 21*t);
        
        if(t == 0){
            wrenchDiffCost += 0.5*(x_map.segment<6>(9) - fl_map).transpose() * derivativeWrenchWeight_map.head<6>().asDiagonal() * (x_map.segment<6>(9) - fl_map);
            wrenchDiffCost += 0.5*(x_map.segment<6>(9+6) - fr_map).transpose() * derivativeWrenchWeight_map.tail<6>().asDiagonal() * (x_map.segment<6>(9+6) - fr_map);
        }
        else{
            wrenchDiffCost += 0.5*(x_map.segment<12>(9+21*t) - x_map.segment<12>(9+21*(t-1))).transpose() * derivativeWrenchWeight_map.asDiagonal() * (x_map.segment<12>(9+21*t) - x_map.segment<12>(9+21*(t-1)));
        }
    }

    double terminalCost = 0.5*(x_map.segment<9>(21*(m_horizon-1)) - iDynTree::toEigen(m_desiredGamma)).transpose() * gammaImpactWeight_map.asDiagonal() * (x_map.segment<9>(21*(m_horizon-1)) - iDynTree::toEigen(m_desiredGamma));

    obj_value = gammaCost + impactCost + wrenchCost + wrenchDiffCost + terminalCost;
    return true;
}

bool MPCIpOptSolver::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{
    Eigen::Map<const Eigen::VectorXd> x_map(x, n);
    Eigen::Map<Eigen::VectorXd> grad_map(grad_f, n);
    Eigen::Map<const Eigen::VectorXd > gammaWeight_map(m_gammaWeight.data(),m_gammaWeight.size());
    Eigen::Map<const Eigen::VectorXd > gammaImpactWeight_map(m_gammaWeightImpact.data(),m_gammaWeightImpact.size());
    Eigen::Map<const Eigen::VectorXd > wrenchWeight_map(m_wrenchWeight.data(),m_wrenchWeight.size());
    Eigen::Map<const Eigen::VectorXd > derivativeWrenchWeight_map(m_derivativeWrenchWeight.data(),m_derivativeWrenchWeight.size());
    Eigen::Map <Eigen::VectorXd> fl_map(m_fLPrev.data(), m_fLPrev.size());
    Eigen::Map <Eigen::VectorXd> fr_map(m_fRPrev.data(), m_fRPrev.size());
    
    for(int t=0; t<m_horizon; ++t){
        //Gamma
        grad_map.segment<9>(21*t) = gammaWeight_map.asDiagonal() * (x_map.segment<9>(21*t) - iDynTree::toEigen(m_desiredGamma));
        //Gamma Impact
        if((t >= m_impact)&&(t < (m_horizon-1))){
            grad_map.segment<9>(21*t) += gammaImpactWeight_map.asDiagonal() * (x_map.segment<9>(21*t) - iDynTree::toEigen(m_desiredGamma));
        }
        //Wrench
        grad_map.segment<12>(9 + 21*t) = wrenchWeight_map.asDiagonal() * x_map.segment<12>(9 + 21*t);
        //Derivative Wrench
        if(t == 0){
            grad_map.segment<6>(9) += derivativeWrenchWeight_map.head<6>().asDiagonal() * (x_map.segment<6>(9) - fl_map);
            grad_map.segment<6>(9+6) += derivativeWrenchWeight_map.tail<6>().asDiagonal() * (x_map.segment<6>(9+6) - fr_map);
        }
        else{
            grad_map.segment<12>(9 + 21*t) += derivativeWrenchWeight_map.asDiagonal() * (x_map.segment<12>(9+21*t) - x_map.segment<12>(9+21*(t-1)));
            grad_map.segment<12>(9 + 21*(t-1)) += -1*derivativeWrenchWeight_map.asDiagonal() * (x_map.segment<12>(9+21*t) - x_map.segment<12>(9+21*(t-1)));
        }
    }
    //terminal
    grad_map.segment<9>(21*(m_horizon-1)) += gammaImpactWeight_map.asDiagonal() * (x_map.segment<9>(21*(m_horizon-1)) - iDynTree::toEigen(m_desiredGamma));
    
    return true;
}

bool MPCIpOptSolver::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g)
{
    Eigen::Map<const Eigen::VectorXd> x_map(x, n);
    Eigen::Map<Eigen::VectorXd> g_map(g, m);
    iDynTree::iDynTreeEigenMatrixMap map_EvGamma = iDynTree::toEigen(m_EvGamma);
    iDynTree::iDynTreeEigenMatrixMap map_FGamma = iDynTree::toEigen(m_FGamma);
    Eigen::Map <Eigen::VectorXd> gamma0_map(m_gamma0.data(), 9);
    Eigen::Map <Eigen::VectorXd> bias_map (m_bias.data(), 9);
    iDynTree::iDynTreeEigenMatrixMap Al_map = iDynTree::toEigen(m_wrenchAl);
    iDynTree::iDynTreeEigenMatrixMap Ar_map = iDynTree::toEigen(m_wrenchAr);
    Eigen::Map <Eigen::VectorXd> b_map(m_wrenchb.data(),m_wrenchb.size());
    Eigen::Map <Eigen::VectorXd> bImpact_map(m_wrenchbImpact.data(),m_wrenchbImpact.size());
    
    
    unsigned int nConstraintsL = m_wrenchAl.rows();
    unsigned int nConstraintsR = m_wrenchAr.rows();
    unsigned int constraintOffset = 0;
    
    for(int t=0; t<m_horizon; ++t){
        if(t==0){
            g_map.head<9>() = -x_map.head<9>() + map_EvGamma*gamma0_map + map_FGamma*x_map.segment<12>(9) + bias_map;
        }
        else{
            g_map.segment<9>(9*t) = -x_map.segment<9>(21*t) + map_EvGamma*x_map.segment<9>(21*(t-1)) + map_FGamma*x_map.segment<12>(9 + 21*t) + bias_map;
        }
        
        constraintOffset = 9*m_horizon + (nConstraintsL+ nConstraintsR)*t;
        if(t < m_impact){
            if(m_rightFootStep){
                g_map.segment(constraintOffset, nConstraintsL) = Al_map*x_map.segment<6>(9 + 21*t) - bImpact_map;
                constraintOffset += nConstraintsL;
                g_map.segment(constraintOffset, nConstraintsR) = Ar_map*x_map.segment<6>(9 + 21*t + 6) - b_map;
            }
            else{
                g_map.segment(constraintOffset, nConstraintsL) = Al_map*x_map.segment<6>(9 + 21*t) - b_map;
                constraintOffset += nConstraintsL;
                g_map.segment(constraintOffset, nConstraintsR) = Ar_map*x_map.segment<6>(9 + 21*t + 6) - bImpact_map;
            }
        }
        else{
            g_map.segment(constraintOffset, nConstraintsL) = Al_map*x_map.segment<6>(9 + 21*t) - bImpact_map;
            constraintOffset += nConstraintsL;
            g_map.segment(constraintOffset, nConstraintsR) = Ar_map*x_map.segment<6>(9 + 21*t + 6) - bImpact_map;
        }
        //std::cerr << "Test frictionL: " << x_map.segment<2>(9 + 21*t).norm()/x_map(9 + 21*t + 2) << std::endl;
        //std::cerr << "Test frictionR: " << x_map.segment<2>(9 + 21*t + 6).norm()/x_map(9 + 21*t + 8) << std::endl;
    }
    //std::cerr << "Remainders constraints: " << g_map.transpose() << std::endl;
    return true;
}

bool MPCIpOptSolver::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{
    if(values == NULL){
        int index = 0;
        temp_iRow.resize(nele_jac);
        temp_jCol.resize(nele_jac);
        
        for(int block = 0; block < m_modelConstraintsJacobian.size(); ++block){
            for(iDynTree::SparseMatrix::Iterator val = m_modelConstraintsJacobian[block].blockPtr->begin(); val != m_modelConstraintsJacobian[block].blockPtr->end(); ++val){
                iRow[index] = val->row() + m_modelConstraintsJacobian[block].rowOffset;
                jCol[index] = val->column() + m_modelConstraintsJacobian[block].colOffset;
                temp_iRow(index) = iRow[index];
                temp_jCol(index) = jCol[index];
                index++;
            }
        }
        int rowOffset = m_horizon*9;
        for(int block=0; block < m_wrenchConstraintJacobian.size(); ++block){
            for(iDynTree::SparseMatrix::Iterator val = m_wrenchConstraintJacobian[block].blockPtr->begin(); val != m_wrenchConstraintJacobian[block].blockPtr->end(); ++val){
                iRow[index] = val->row() + m_wrenchConstraintJacobian[block].rowOffset + rowOffset;
                jCol[index] = val->column() + m_wrenchConstraintJacobian[block].colOffset;
                temp_iRow(index) = iRow[index];
                temp_jCol(index) = jCol[index];
                index++;
            }
        }
    }
    else{
        int index = 0;
        for(int block = 0; block < m_modelConstraintsJacobian.size(); ++block){
            for(iDynTree::SparseMatrix::Iterator val = m_modelConstraintsJacobian[block].blockPtr->begin(); val != m_modelConstraintsJacobian[block].blockPtr->end(); ++val){
                values[index] = val->value();
                index++;
            }
        }
        for(int block=0; block < m_wrenchConstraintJacobian.size(); ++block){
            for(iDynTree::SparseMatrix::Iterator val = m_wrenchConstraintJacobian[block].blockPtr->begin(); val != m_wrenchConstraintJacobian[block].blockPtr->end(); ++val){
                values[index] = val->value();
                index++;
            }
        }
        
       // printIpOptMatrix("Constraints Jacobian:", nele_jac, temp_iRow, temp_jCol, values, m, n);
        
    }
    return true;
}


bool MPCIpOptSolver::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{
    if(values == NULL){
        int index = 0;
        for(iDynTree::SparseMatrix::Iterator val = m_costHessian.begin(); val != m_costHessian.end(); ++val){
            iRow[index] = val->row();
            jCol[index] = val->column();
            index++;
        }
    }
    else{
        int index = 0;
        for(iDynTree::SparseMatrix::Iterator val = m_costHessian.begin(); val != m_costHessian.end(); ++val){
            values[index] = obj_factor*val->value();
            index++;
        }
    }
    return true;
}

void MPCIpOptSolver::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    if((status == Ipopt::SUCCESS)||status == Ipopt::STOP_AT_ACCEPTABLE_POINT){
        Eigen::Map< const Eigen::VectorXd > x_map (x, n);

        m_previousSolution.resize(n);
        iDynTree::toEigen(m_previousSolution) = x_map;

        Eigen::Map<const Eigen::VectorXd> zlMap(z_L, n);
        Eigen::Map<const Eigen::VectorXd> zuMap(z_U, n);
        iDynTree::toEigen(m_lowerBoundMultipliers) = zlMap;
        iDynTree::toEigen(m_upperBoundMultipliers) = zuMap;

        Eigen::Map<const Eigen::VectorXd> lambdaMap(lambda, m);
        iDynTree::toEigen(m_contraintMultipliers) = lambdaMap;


//        std::cerr << "Final solution" << std::endl;
//        std::cerr << zlMap.transpose() << std::endl<< std::endl;
//        std::cerr << zuMap.transpose() << std::endl<< std::endl;
//        std::cerr << lambdaMap.transpose() << std::endl<< std::endl;
//        std::cerr << x_map.transpose() << std::endl<< std::endl;

        if(status == Ipopt::SUCCESS){
            m_exitCode = 0;
        }
        else m_exitCode = 1;
    }
    else {
        switch(status){
            
            case Ipopt::LOCAL_INFEASIBILITY:
                m_exitCode = -1;
                break;
                
            case Ipopt::MAXITER_EXCEEDED:
            case Ipopt::CPUTIME_EXCEEDED:
            case Ipopt::STOP_AT_TINY_STEP: //PREMATURE STOP
                m_exitCode = -2;
                break;
                
            case Ipopt::INVALID_NUMBER_DETECTED: //WRONG DATA INSERTION
                m_exitCode = -3;
                break;
                
            case Ipopt::DIVERGING_ITERATES:
            case Ipopt::INTERNAL_ERROR:
            case Ipopt::ERROR_IN_STEP_COMPUTATION:
            case Ipopt::RESTORATION_FAILURE: //IPOPT INTERNAL PROBLEM
                m_exitCode = -4;
                break;
                
            case Ipopt::USER_REQUESTED_STOP: //USER REQUESTED STOP
                m_exitCode = -5;
                break;
                
            default:
                m_exitCode = -6;
                
        }
    }
}

bool MPCIpOptSolver::get_warm_start_iterate(Ipopt::IteratesVector& warm_start_iterate)
{
    std::cerr << "Hello inside warm start\n";
    return true;
}

void MPCIpOptSolver::printJacobian(const std::string header, const std::vector<MatrixBlock>& input, int rows, int cols)
{
    iDynTree::SparseMatrix printer;
    printer.resize(rows, cols);
    
    iDynTree::Triplets filler;
    for(int block = 0; block < input.size(); ++block){
        for(iDynTree::SparseMatrix::Iterator val = input[block].blockPtr->begin(); val != input[block].blockPtr->end(); ++val){
            filler.pushTriplet(iDynTree::Triplet(val->row()+input[block].rowOffset, val->column()+input[block].colOffset, val->value()));
        }
    }
    printer.setFromTriplets(filler);
    std::cerr << header << std::endl << printer.description(true) << std::endl;
}

void MPCIpOptSolver::printIpOptMatrix(const std::string header, Ipopt::Index nele_jac, const iDynTree::VectorDynSize& iRow, const iDynTree::VectorDynSize& jCol, Ipopt::Number* values, int rows, int cols)
{
    iDynTree::SparseMatrix printer;
    printer.resize(rows, cols);
    
    iDynTree::Triplets filler;
    for(int i=0; i < nele_jac; ++i){
        filler.pushTriplet(iDynTree::Triplet(iRow(i),jCol(i),values[i]));
    }
    printer.setFromTriplets(filler);
    std::cerr << header << std::endl << printer.description(true) << std::endl;
}

