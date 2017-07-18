
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

class StepRecoveryMPC {
    
    Ipopt::SmartPtr< MPCIpOptSolver > solverPointer;
    Ipopt::SmartPtr< Ipopt::IpoptApplication > loader;
    std::vector< std::pair<double, double> > m_footDimensions;
    bool m_configured;
    
    bool getVectorFromValue(yarp::os::Value& input, iDynTree::VectorDynSize& output);
    
    bool computeWrenchConstraints();
    
    int dryRun();
    int optimize();
    int reOptimize();
    
    bool getCoMPosition();
    bool setPreviousWrench();
    
    bool setLeftFootConstraints();
    bool setRightFootConstraints();
    
    bool updateConstraints();
    
public:
    
    StepRecoveryMPC();
    ~StepRecoveryMPC();
    
    bool solve();
    
    bool getControllerData();
    
    bool configure(yarp::os::Searchable& mpcOptions);
    
};

#endif
