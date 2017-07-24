
/******************************************************************************
 *                                                                            *
 * Copyright (C) 2017 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @file MPCIpOptSolverImplicit.cpp
 * @authors: Stefano Dafarra <stefano.dafarra@iit.it>
 */

#include <MPCIpOptSolverImplicit.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/EigenSparseHelpers.h>
#include <Eigen/Core>

MPCIpOptSolverImplicit::MPCIpOptSolverImplicit()
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

MPCIpOptSolverImplicit::~MPCIpOptSolverImplicit()
{
}

bool MPCIpOptSolverImplicit::setTimeSettings(double dT, unsigned int horizon)
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
    m_costHessian.resize(12*m_horizon, 12*m_horizon);
    m_iteratedModel.resize(m_horizon);
    m_controllabilityMatrix.resize(m_horizon);
    m_biasImplicit.resize(m_horizon);
    m_iteratedGamma0.resize(m_horizon);
    m_gammaBuffer.resize(m_horizon);
    m_summedControllability.resize(m_horizon);
    for(int t=0; t<m_horizon; ++t){
        m_biasImplicit[t].resize(9);
        m_iteratedGamma0[t].resize(9);
        m_gammaBuffer[t].resize(9);
    }
    
    return true;
}

bool MPCIpOptSolverImplicit::setRobotMass(const double mass)
{
    if(mass <= 0){
        std::cerr << "The mass is expected to be a positive number." << std::endl;
        return false;
    }
    m_mass = mass;
    return true;
}

bool MPCIpOptSolverImplicit::setImpactInstant(unsigned int impact, bool rightFootStepping)
{
    m_impact = impact;
    m_rightFootStep = rightFootStepping;
    return true;
}

bool MPCIpOptSolverImplicit::setGamma0(const iDynTree::VectorFixSize<9>& gamma0)
{
    m_gamma0 = gamma0;
    //std::cerr << "Gamma0: " << m_gamma0.toString() << std::endl;
    return true;
}

bool MPCIpOptSolverImplicit::setPreviousWrench(const iDynTree::VectorDynSize& previousLeftWrench, const iDynTree::VectorDynSize& previousRightWrench)
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


bool MPCIpOptSolverImplicit::setWrenchConstraints(const iDynTree::MatrixDynSize& wrenchConstraintsMatrix, const iDynTree::VectorDynSize& constraintsBounds, const iDynTree::VectorDynSize& afterImpactConstraintsBounds)
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


bool MPCIpOptSolverImplicit::setFeetTransforms(const iDynTree::Transform& w_H_l, const iDynTree::Transform& w_H_r)
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


bool MPCIpOptSolverImplicit::setDesiredCOMPosition(const iDynTree::Position& desiredCOM)
{
    Eigen::Map<const Eigen::VectorXd> desiredCOM_map(desiredCOM.data(), 3);
    iDynTree::toEigen(m_desiredGamma).head<3>() = desiredCOM_map;
    
    std::cerr << "desiredGamma: "<< m_desiredGamma.toString() << std::endl;
    
    return true;
}

bool MPCIpOptSolverImplicit::setGammaWeight(const iDynTree::VectorDynSize& gammaWeight)
{
    if(gammaWeight.size() != 9){
        std::cerr << "The gammaWeight vector is expected to be 9 dimensional" << std::endl;
        return false;
    }
    m_gammaWeight = gammaWeight;
    
    //  std::cerr << "gammaWeight: "<< m_gammaWeight.toString() << std::endl;
    
    return true;
}

bool MPCIpOptSolverImplicit::setPostImpactGammaWeight(const iDynTree::VectorDynSize& gammaImpactWeight)
{
    if(gammaImpactWeight.size() != 9){
        std::cerr << "The gammaWeight vector is expected to be 9 dimensional" << std::endl;
        return false;
    }
    
    m_gammaWeightImpact = gammaImpactWeight;
    
    // std::cerr << "gammaWeightImpact: "<< m_gammaWeightImpact.toString() << std::endl;
    
    return true;
}

bool MPCIpOptSolverImplicit::setWrenchsWeight(const iDynTree::VectorDynSize& wrenchWeight)
{
    if(wrenchWeight.size() != 12){
        std::cerr << "The wrenchWeight vector is expected to be 12 dimensional" << std::endl;
        return false;
    }
    m_wrenchWeight = wrenchWeight;
    
    // std::cerr << "wrenchWeight: "<< m_wrenchWeight.toString() << std::endl;
    
    return true;
}

bool MPCIpOptSolverImplicit::setWrenchDerivativeWeight(const iDynTree::VectorDynSize& derivativeWrenchWeight)
{
    if(derivativeWrenchWeight.size() != 12){
        std::cerr << "The derivativeWrenchWeight vector is expected to be 12 dimensional" << std::endl;
        return false;
    }
    m_derivativeWrenchWeight = derivativeWrenchWeight;
    
    //std::cerr << "wrenchWeight: "<< derivativeWrenchWeight.toString() << std::endl;
    
    
    return true;
}

bool MPCIpOptSolverImplicit::computeModelMatrices()
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
    
    m_EvGammaEigen = iDynTree::toEigen(*m_EvGammaSparsePtr);
    
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
    
    m_FGammaEigen = iDynTree::toEigen(*m_FGammaSparsePtr);
    
    
    //     std::cerr << "Ev: "<< std::endl << m_EvGamma.toString() << std::endl;
    //    std::cerr << "F: " << std::endl << m_FGamma.toString() << std::endl;
    //    std::cerr << "FSparse: " << std::endl << m_FGammaSparsePtr->description(true) << std::endl;
    
    
    return true;
}

bool MPCIpOptSolverImplicit::computeModelBias()
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

bool MPCIpOptSolverImplicit::computeWrenchConstraints()
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

bool MPCIpOptSolverImplicit::computeWrenchConstraintsJacobian()
{
    MatrixBlock templateLeft, templateRight;
    
    templateLeft.blockPtr = m_wrenchAlSparsePtr;
    templateRight.blockPtr =  m_wrenchArSparsePtr;
    
    m_wrenchConstraintJacobian.resize(m_horizon*2);
    
    unsigned int row = 0;
    unsigned int col = 0;
    unsigned int index = 0;
    
    unsigned int nConstraintsL = m_wrenchAl.rows();
    unsigned int nConstraintsR = m_wrenchAr.rows();
    
    for (int t=0; t<m_horizon; ++t){
        col = 12*t; //you have to pick the right variables inside the decision variables vector
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

bool MPCIpOptSolverImplicit::computeSingleCostHessian()
{
    m_gammaWeightHessian.resize(9,9);
    m_gammaWeightImpactHessian.resize(9,9);
    m_wrenchWeightHessian.resize(12,12);
    m_derivativeWrenchWeightHessian.resize(12,12);
    m_negativeDerWrenchHessian.resize(12,12);
    
    m_gammaWeightHessian = iDynTree::toEigen(m_gammaWeight).asDiagonal();
    
    m_gammaWeightImpactHessian = iDynTree::toEigen(m_gammaWeightImpact).asDiagonal();
    m_gammaWeightImpactHessian.makeCompressed();

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

iDynTree::SparseMatrix MPCIpOptSolverImplicit::convertSparse(Eigen::SparseMatrix<double>& input)
{
    iDynTree::SparseMatrix output;
    output.resize(input.rows(), input.cols());
    iDynTree::Triplets values;
    values.reserve(input.nonZeros());
    for(int k=0; k<input.outerSize(); ++k){
        for(Eigen::SparseMatrix<double>::InnerIterator it(input, k); it; ++it){
            values.pushTriplet(iDynTree::Triplet(it.row(),it.col(),it.value()));
        }
    }
    output.setFromTriplets(values);
    return output;
}


bool MPCIpOptSolverImplicit::computeCostHessian()
{
    m_summedControllability[m_horizon-1] = m_controllabilityMatrix[m_horizon - 1];
    for(int t=m_horizon-2; t >= 0; --t){
        m_summedControllability[t] = m_summedControllability[t+1] + m_controllabilityMatrix[t];
    }
    //m_summedControllability contains the gradient of sum_over_t (gamma(t)/df(i))
    //e.g. f(m_horizon-1) is contained only in gamma(m_horizon) with the matrix m_controllability[0] = F_gamma
    //     f(m_horizon-2) is contained in gamma(m_horizon) with the matrix Ev_gamma*F_Gamma and in gamma(m_horizon-1) with F_Gamma
    //     f(0) is present in every gamma, thus the gradient is the sum of all the elements of the controllability matrix
    
    Eigen::SparseMatrix<double> temp;
    
    iDynTree::Triplets values;
    values.reserve((9*9*3 + 12*5)*m_horizon); //HUUUGE
    
    for(int t=0; t<m_horizon; ++t){
        temp = m_summedControllability[t].transpose()*m_gammaWeightHessian*m_summedControllability[t]; //running cost on gamma
        values.addSubMatrix(12*t, 12*t, convertSparse(temp));
        
        if(t >= m_impact){
            temp = m_summedControllability[t].transpose()*m_gammaWeightImpactHessian*m_summedControllability[t]; //the forces after the impact have a full contribution over the post impact states
            values.addSubMatrix(12*t, 12*t, convertSparse(temp));
        }
        else{
            temp = m_controllabilityMatrix[t].transpose()*m_gammaWeightImpactHessian*m_controllabilityMatrix[t]; //terminal cost, every force has a contribution equal to the corresponding term in the controllability matrix.
            if(m_impact < m_horizon){ //if not, we have already considered the terminal cost
                    if(t < m_horizon-1){ //if not, we have already considered the terminal cost
                        for(int k = t+1; k < (m_horizon - (m_impact - t)); ++k){ //I have to sum the m_horizon - impact elements starting from t. But the first element has already been summed up in the terminal cost
                            temp += m_controllabilityMatrix[k].transpose()*m_gammaWeightImpactHessian*m_controllabilityMatrix[k]; 
                        }
                }
            }
            values.addSubMatrix(12*t, 12*t, convertSparse(temp));
        }
        values.addSubMatrix(12*t, 12*t, m_wrenchWeightHessian);
        if(t == 0){
            values.addSubMatrix(0, 0, m_derivativeWrenchWeightHessian);
        }
        else {
            values.addSubMatrix(12*t, 12*t, m_derivativeWrenchWeightHessian);
            values.addSubMatrix(12*t, 12*(t-1), m_negativeDerWrenchHessian);
            values.addSubMatrix(12*(t-1), 12*(t-1), m_derivativeWrenchWeightHessian);
            values.addSubMatrix(12*(t-1), 12*t, m_negativeDerWrenchHessian);
        }
    }
    m_costHessian.setFromTriplets(values);
    //std::cerr << "Cost hessian: " << std::endl << m_costHessian.description(true) << std::endl;
    
    return true;
}

bool MPCIpOptSolverImplicit::computeIteratedModelMatrices()
{
    m_iteratedModel[0] = m_EvGammaEigen;
    
    m_controllabilityMatrix[m_horizon-1] = m_FGammaEigen;
    
    m_biasImplicit[0] = iDynTree::toEigen(m_bias);
    
    Eigen::Map <Eigen::VectorXd> gamma0_map(m_gamma0.data(), 9);
    m_iteratedGamma0[0] = m_EvGammaEigen*gamma0_map;
    
    for(int t=1; t<m_horizon; ++t){
        m_iteratedModel[t] = m_iteratedModel[t-1]*m_EvGammaEigen;
        m_controllabilityMatrix[m_horizon-1 - t] = m_iteratedModel[t-1]*m_FGammaEigen; //[a^(n-1)b a^(n-2)b ... ab b]
        m_biasImplicit[t] = m_iteratedModel[t-1]*iDynTree::toEigen(m_bias);
        m_iteratedGamma0 [t] = m_iteratedModel[t]*gamma0_map;
    }
    
    return true;
}

bool MPCIpOptSolverImplicit::updateGamma(const Eigen::VectorXd& solution)
{
    for(int t=0; t<m_horizon; ++t){
        m_gammaBuffer[t] = m_iteratedGamma0[t];
        for(int k = 0; k <= t; ++k){
            m_gammaBuffer[t] += m_controllabilityMatrix[m_horizon -1 -t + k]*solution.segment<12>(12*k) + m_biasImplicit[k];
        }
    }
    return true;
}


bool MPCIpOptSolverImplicit::prepareProblem()
{
    if(!computeModelMatrices()){
        std::cerr << "Error while computing the model matrices." << std::endl;
        return false;
    }
    if(!computeModelBias()){
        std::cerr << "Error while computing the model bias." << std::endl;
        return false;
    }
    if(!computeIteratedModelMatrices()){
        std::cerr << "Error while the iterated model matrices." << std::endl;
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
    if(!computeSingleCostHessian()){
        std::cerr << "Error while computing the single cost hessian." << std::endl;
        return false;
    }
    if(!computeCostHessian()){
        std::cerr << "Error while computing the cost hessian." << std::endl;
        return false;
    }
    
    return true;
}


bool MPCIpOptSolverImplicit::updateProblem()
{
    if(!computeModelMatrices()){
        std::cerr << "Error while computing the model matrices." << std::endl;
        return false;
    }
    if(!computeModelBias()){
        std::cerr << "Error while computing the model bias." << std::endl;
        return false;
    }
    if(!computeIteratedModelMatrices()){
        std::cerr << "Error while the iterated model matrices." << std::endl;
        return false;
    }
    if(!computeWrenchConstraints()){
        std::cerr << "Error while computing the wrench constraints." << std::endl;
        return false;
    }
    if(!computeCostHessian()){
        std::cerr << "Error while computing the cost hessian." << std::endl;
        return false;
    }
    //   printJacobian("Model jacobian:", m_modelConstraintsJacobian, 9*m_horizon, 21*m_horizon);
    //    printJacobian("Wrench jacobian:", m_wrenchConstraintJacobian, (m_wrenchAl.rows()+m_wrenchAr.rows())*m_horizon, 21*m_horizon);
    
    
    return true;
}

int MPCIpOptSolverImplicit::getSolution(iDynTree::VectorDynSize& fL, iDynTree::VectorDynSize& fR, iDynTree::VectorDynSize& lastGamma)
{
    if(m_previousSolution.size() >= 12){
        Eigen::Map<Eigen::VectorXd> solution_map(m_previousSolution.data(),m_previousSolution.size());
        Eigen::Map <Eigen::VectorXd> gamma0_map(m_gamma0.data(), 9);
        
        fL.resize(6);
        fR.resize(6);
        lastGamma.resize(9);
        iDynTree::toEigen(fL) = solution_map.head<6>();
        iDynTree::toEigen(fR) = solution_map.segment<6>(6);
        iDynTree::toEigen(lastGamma) = m_iteratedGamma0[m_horizon-1];
        for(int t=0; t < m_horizon; ++t){
            iDynTree::toEigen(lastGamma) += m_controllabilityMatrix[t]*solution_map.segment<12>(12*t) + m_biasImplicit[t];
        }
    }
    
    // std::cerr << "Solution left: " << fL.toString() << std::endl;
    //std::cerr << "Solution right: " << fR.toString() << std::endl;
    //std::cerr << "Solution lastGamma: " << lastGamma.toString() << std::endl;
    // std::cerr << "Full gamma: "<<m_previousSolution.toString() << std::endl;
    
    return m_exitCode;
}


bool MPCIpOptSolverImplicit::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style)
{
    n = 12*m_horizon;
    m = (m_wrenchAl.rows()+m_wrenchAr.rows())*m_horizon; //TODO FIX
    nnz_jac_g = 0;
    for(int i=0; i<m_wrenchConstraintJacobian.size(); ++i){
        nnz_jac_g += m_wrenchConstraintJacobian[i].blockPtr->numberOfNonZeros(); //TODO FIX
    }
    nnz_h_lag = m_costHessian.numberOfNonZeros();
    
    index_style = Ipopt::TNLP::C_STYLE;
    
    return true;
}

bool MPCIpOptSolverImplicit::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u)
{
    for (Ipopt::Index i = 0; i < n; ++i) {
        x_l[i] = -2e+19;
        x_u[i] =  2e+19;
    }
    
    for (Ipopt::Index c = 0; c < m; ++c) {
        g_l[c] = -2e+19;
        g_u[c] =  0;
    }
    
    return true;
}

bool MPCIpOptSolverImplicit::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda)
{
    if(init_z) return false;
    if(init_lambda) return false;
    
    if(init_x){
        Eigen::Map<Eigen::VectorXd> x_map(x, n);
        
        x_map.setZero();
        for(int t=0; t < m_horizon; ++t){
            x_map(12*t + 2) = 20;
        }
        
        return true;
    }
    return false;
}

bool MPCIpOptSolverImplicit::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value)
{
    Eigen::Map<const Eigen::VectorXd> x_map(x, n);
    
    if(new_x){
        if(! updateGamma(x_map))
            return false;
    }
    
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
        gammaCost += 0.5*(m_gammaBuffer[t] - iDynTree::toEigen(m_desiredGamma)).transpose() * gammaWeight_map.asDiagonal() * (m_gammaBuffer[t] - iDynTree::toEigen(m_desiredGamma));
        
        if((t >= m_impact)&&(t < (m_horizon-1))){
            impactCost += 0.5*(m_gammaBuffer[t] - iDynTree::toEigen(m_desiredGamma)).transpose() * gammaImpactWeight_map.asDiagonal() * (m_gammaBuffer[t] - iDynTree::toEigen(m_desiredGamma));
        }
        //         std::cerr << "Gamma: " << x_map.segment<9>(21*t).transpose() << std::endl;
        //         std::cerr << "Gamma Des: " << m_desiredGamma.toString() << std::endl;
        //         std::cerr << "Gamma error: "<< x_map.segment<9>(21*t).transpose() - iDynTree::toEigen(m_desiredGamma).transpose() << std::endl;
        
        wrenchCost += 0.5*x_map.segment<12>(12*t).transpose() * wrenchWeight_map.asDiagonal() * x_map.segment<12>(12*t);
        
        if(t == 0){
            wrenchDiffCost += 0.5*(x_map.segment<6>(0) - fl_map).transpose() * derivativeWrenchWeight_map.head<6>().asDiagonal() * (x_map.segment<6>(0) - fl_map);
            wrenchDiffCost += 0.5*(x_map.segment<6>(6) - fr_map).transpose() * derivativeWrenchWeight_map.tail<6>().asDiagonal() * (x_map.segment<6>(6) - fr_map);
        }
        else{
            wrenchDiffCost += 0.5*(x_map.segment<12>(12*t) - x_map.segment<12>(12*(t-1))).transpose() * derivativeWrenchWeight_map.asDiagonal() * (x_map.segment<12>(12*t) - x_map.segment<12>(12*(t-1)));
        }
    }
    
    double terminalCost = 0.5*(m_gammaBuffer[m_horizon-1] - iDynTree::toEigen(m_desiredGamma)).transpose() * gammaImpactWeight_map.asDiagonal() * (m_gammaBuffer[m_horizon-1] - iDynTree::toEigen(m_desiredGamma));
    
    obj_value = gammaCost + impactCost + wrenchCost + wrenchDiffCost + terminalCost;
    return true;
}

bool MPCIpOptSolverImplicit::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{
    Eigen::Map<const Eigen::VectorXd> x_map(x, n);
    
    if(new_x){
        if(! updateGamma(x_map))
            return false;
    }
    
    Eigen::Map<Eigen::VectorXd> grad_map(grad_f, n);
    grad_map.setZero();
    
    Eigen::Map<const Eigen::VectorXd > gammaWeight_map(m_gammaWeight.data(),m_gammaWeight.size());
    Eigen::Map<const Eigen::VectorXd > gammaImpactWeight_map(m_gammaWeightImpact.data(),m_gammaWeightImpact.size());
    Eigen::Map<const Eigen::VectorXd > wrenchWeight_map(m_wrenchWeight.data(),m_wrenchWeight.size());
    Eigen::Map<const Eigen::VectorXd > derivativeWrenchWeight_map(m_derivativeWrenchWeight.data(),m_derivativeWrenchWeight.size());
    Eigen::Map <Eigen::VectorXd> fl_map(m_fLPrev.data(), m_fLPrev.size());
    Eigen::Map <Eigen::VectorXd> fr_map(m_fRPrev.data(), m_fRPrev.size());
    
    for(int t=0; t<m_horizon; ++t){
        for(int k = t; k < m_horizon; k++){
            grad_map.segment<12>(12*t) += m_controllabilityMatrix[k].transpose()*gammaWeight_map.asDiagonal()*(m_gammaBuffer[k] - iDynTree::toEigen(m_desiredGamma));
            if((k >= m_impact)&&(k < m_horizon - 1)){
                grad_map.segment<12>(12*t) += m_controllabilityMatrix[k].transpose()*gammaImpactWeight_map.asDiagonal()*(m_gammaBuffer[k] - iDynTree::toEigen(m_desiredGamma));
            }
        }
        grad_map.segment<12>(12*t) += m_controllabilityMatrix[t].transpose()*gammaImpactWeight_map.asDiagonal()*(m_gammaBuffer[m_horizon-1] - iDynTree::toEigen(m_desiredGamma));
        //Wrench
        grad_map.segment<12>(12*t) += wrenchWeight_map.asDiagonal() * x_map.segment<12>(12*t);
        //Derivative Wrench
        if(t == 0){
            grad_map.segment<6>(0) += derivativeWrenchWeight_map.head<6>().asDiagonal() * (x_map.segment<6>(0) - fl_map);
            grad_map.segment<6>(6) += derivativeWrenchWeight_map.tail<6>().asDiagonal() * (x_map.segment<6>(6) - fr_map);
        }
        else{
            grad_map.segment<12>(12*t) += derivativeWrenchWeight_map.asDiagonal() * (x_map.segment<12>(12*t) - x_map.segment<12>(12*(t-1)));
            grad_map.segment<12>(12*(t-1)) += -1*derivativeWrenchWeight_map.asDiagonal() * (x_map.segment<12>(12*t) - x_map.segment<12>(12*(t-1)));
        }
    }
    
    return true;
}

bool MPCIpOptSolverImplicit::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g)
{
    Eigen::Map<const Eigen::VectorXd> x_map(x, n);
    Eigen::Map<Eigen::VectorXd> g_map(g, m);
    iDynTree::iDynTreeEigenMatrixMap Al_map = iDynTree::toEigen(m_wrenchAl);
    iDynTree::iDynTreeEigenMatrixMap Ar_map = iDynTree::toEigen(m_wrenchAr);
    Eigen::Map <Eigen::VectorXd> b_map(m_wrenchb.data(),m_wrenchb.size());
    Eigen::Map <Eigen::VectorXd> bImpact_map(m_wrenchbImpact.data(),m_wrenchbImpact.size());
    
    
    unsigned int nConstraintsL = m_wrenchAl.rows();
    unsigned int nConstraintsR = m_wrenchAr.rows();
    unsigned int constraintOffset = 0;
    
    for(int t=0; t<m_horizon; ++t){
        constraintOffset = (nConstraintsL+ nConstraintsR)*t;
        if(t < m_impact){
            if(m_rightFootStep){
                g_map.segment(constraintOffset, nConstraintsL) = Al_map*x_map.segment<6>(12*t) - bImpact_map;
                constraintOffset += nConstraintsL;
                g_map.segment(constraintOffset, nConstraintsR) = Ar_map*x_map.segment<6>(12*t + 6) - b_map;
            }
            else{
                g_map.segment(constraintOffset, nConstraintsL) = Al_map*x_map.segment<6>(12*t) - b_map;
                constraintOffset += nConstraintsL;
                g_map.segment(constraintOffset, nConstraintsR) = Ar_map*x_map.segment<6>(12*t + 6) - bImpact_map;
            }
        }
        else{
            g_map.segment(constraintOffset, nConstraintsL) = Al_map*x_map.segment<6>(12*t) - bImpact_map;
            constraintOffset += nConstraintsL;
            g_map.segment(constraintOffset, nConstraintsR) = Ar_map*x_map.segment<6>(12*t + 6) - bImpact_map;
        }
        //std::cerr << "Test frictionL: " << x_map.segment<2>(9 + 21*t).norm()/x_map(9 + 21*t + 2) << std::endl;
        //std::cerr << "Test frictionR: " << x_map.segment<2>(9 + 21*t + 6).norm()/x_map(9 + 21*t + 8) << std::endl;
    }
    //std::cerr << "Remainders constraints: " << g_map.transpose() << std::endl;
    return true;
}

bool MPCIpOptSolverImplicit::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
{
    if(values == NULL){
        int index = 0;
        temp_iRow.resize(nele_jac);
        temp_jCol.resize(nele_jac);
        
        for(int block=0; block < m_wrenchConstraintJacobian.size(); ++block){
            for(iDynTree::SparseMatrix::Iterator val = m_wrenchConstraintJacobian[block].blockPtr->begin(); val != m_wrenchConstraintJacobian[block].blockPtr->end(); ++val){
                iRow[index] = val->row() + m_wrenchConstraintJacobian[block].rowOffset;
                jCol[index] = val->column() + m_wrenchConstraintJacobian[block].colOffset;
                temp_iRow(index) = iRow[index];
                temp_jCol(index) = jCol[index];
                index++;
            }
        }
    }
    else{
        int index = 0;
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


bool MPCIpOptSolverImplicit::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values)
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

void MPCIpOptSolverImplicit::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    if((status == Ipopt::SUCCESS)||status == Ipopt::STOP_AT_ACCEPTABLE_POINT){
        Eigen::Map< const Eigen::VectorXd > x_map (x, n);
        
        m_previousSolution.resize(n);
        iDynTree::toEigen(m_previousSolution) = x_map;
        
        
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

void MPCIpOptSolverImplicit::printJacobian(const std::string header, const std::vector<MatrixBlock>& input, int rows, int cols)
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

void MPCIpOptSolverImplicit::printIpOptMatrix(const std::string header, Ipopt::Index nele_jac, const iDynTree::VectorDynSize& iRow, const iDynTree::VectorDynSize& jCol, Ipopt::Number* values, int rows, int cols)
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

