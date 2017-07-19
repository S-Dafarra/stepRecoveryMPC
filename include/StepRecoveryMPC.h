
/******************************************************************************
 *                                                                            *
 * Copyright (C) 2017 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @file StepRecoveryMPC.h
 * @authors: Stefano Dafarra <stefano.dafarra@iit.it>
 */

#ifndef STEPRECOVERYMPC_H
#define STEPRECOVERYMPC_H

#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <iDynTree/Core/Transform.h>
#include <iDynTree/Core/SparseMatrix.h>
#include <yarp/os/Searchable.h>
#include "MPCIpOptSolver.h"
#include <IpIpoptApplication.hpp>

enum FootState{
    Standing,
    Swinging,
    Floating
    };
    
typedef struct{
    iDynTree::Transform leftTransform, rightTransform;
    iDynTree::VectorFixSize<9> gamma; //TO BE REMOVED
    unsigned int kImpact; //TO BE REMOVED
    int controllerState;
    double comZDes;
    double robotMass; //TO BE REMOVED
    unsigned int expectedControllerDim = 27;
} ControllerData;

class StepRecoveryMPC {
    
    Ipopt::SmartPtr< MPCIpOptSolver > solverPointer;
    Ipopt::SmartPtr< Ipopt::IpoptApplication > loader;
    std::vector< std::pair<double, double> > m_footDimensions;
    double m_dT, m_stepDuration;
    unsigned int m_horizon;
    
    bool m_configured, m_reOptimize;
    
    double m_frictionCoeff, m_torsionalFrictionCoeff;
    int m_edgesPyramid;
    double m_fzMin, m_fzMax;
    iDynTree::MatrixDynSize m_wrenchConstrMatrix;
    iDynTree::VectorDynSize m_wrenchConstrVector, m_afterImpactwrenchConstrVector;
    
    FootState m_stateL, m_stateR;
    ControllerData m_inputData;
    unsigned int m_kImpact;
    double m_robotMass;
    iDynTree::VectorFixSize<9> m_gamma0;
    iDynTree::VectorDynSize m_prevL, m_prevR;
    
    bool getVectorFromValue(yarp::os::Value& input, iDynTree::VectorDynSize& output);
    
    bool computeWrenchConstraints();
    bool getControllerData(const iDynTree::VectorDynSize& controllerData);
    bool getGamma();
    bool computeImpactInstant();
    bool setDesiredCoM();
    bool setPreviousWrench();
    
    int dryRun();
    
public:
    
    StepRecoveryMPC();
    ~StepRecoveryMPC();
    
    bool solve(const iDynTree::VectorDynSize& controllerData, iDynTree::VectorDynSize& fL, iDynTree::VectorDynSize& fR, iDynTree::VectorDynSize& lastGamma);
    
    bool configure(yarp::os::Searchable& mpcOptions);
    
    void setVerbosity(unsigned int level);
};

#endif
