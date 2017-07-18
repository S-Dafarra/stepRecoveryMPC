/******************************************************************************
 *                                                                            *
 * Copyright (C) 2017 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @file StepRecoveryMPC.cpp
 * @authors: Stefano Dafarra <stefano.dafarra@iit.it>
 */

#include "StepRecoveryMPC.h"
#include <yarp/os/all.h>
#include <yarp/dev/all.h>
#include <yarp/sig/all.h>

StepRecoveryMPC::StepRecoveryMPC()
:m_configured(false)
{
    solverPointer = new MPCIpOptSolver();
    loader = IpoptApplicationFactory();
    //TODO INSERT CONSTANT JACOBIAN AND HESSIAN
}

StepRecoveryMPC::~StepRecoveryMPC()
{
}

bool StepRecoveryMPC::getVectorFromValue(yarp::os::Value& input, iDynTree::VectorDynSize& output)
{
    if(input.isNull()){
        std::cerr << "Empty input value." << std::endl;
        return false;
    }
    if (!input.isList() || !input.asList()){
        std::cerr << "Unable to read the input list." << std::endl;
        return false;
    }
    yarp::os::Bottle *inputPtr = input.asList();
    
    output.resize(inputPtr->size());
    
    for(int i = 0; i < inputPtr->size(); ++i){
        if(!inputPtr->get(i).isDouble() && !inputPtr->get(i).isInt()){
            std::cerr << "The input is expected to be a double" << std::endl;
        }
        output(i) = inputPtr->get(i).asDouble();
    }
    return true;
}


bool StepRecoveryMPC::configure(yarp::os::Searchable& mpcOptions)
{
    double dT = mpcOptions.check("dT", yarp::os::Value(0.01)).asDouble();
    int horizon = mpcOptions.check("horizon", yarp::os::Value(25)).asInt();
    
    yarp::os::Value feetDimensions = mpcOptions.find("foot_size");
    if (feetDimensions.isNull() || !feetDimensions.isList())
    {
        std::cerr << "Unable to find feet dimensions in the configuration file." <<std::endl;
        return false;
    }
    
    yarp::os::Bottle *feetDimensionsPointer = feetDimensions.asList();
    if (!feetDimensionsPointer || feetDimensionsPointer->size() != 2)
    {
        std::cerr << "Error while reading the feet dimensions. Wrong number of elements." <<std::endl;
        return false;
    }
    
    yarp::os::Value& xLimits = feetDimensionsPointer->get(0);
    if (xLimits.isNull() || !xLimits.isList())
    {
        std::cerr << "Error while reading the X limits." << std::endl;
        return false;
    }
    
    yarp::os::Bottle *xLimitsPtr = xLimits.asList();
    if (!xLimitsPtr || xLimitsPtr->size() != 2)
    {
        std::cerr << "Error while reading the X limits. Wrong dimensions." << std::endl;
        return false;
    }
    double limit1 = xLimitsPtr->get(0).asDouble();
    double limit2 = xLimitsPtr->get(1).asDouble();
    m_footDimensions.resize(2);
    m_footDimensions[0].first = std::min(limit1, limit2);
    m_footDimensions[0].second = std::max(limit1, limit2);
    
    yarp::os::Value& yLimits = feetDimensionsPointer->get(1);
    if (yLimits.isNull() || !yLimits.isList())
    {
        std::cerr << "Error while reading the Y limits." << std::endl;
        
        return false;
    }
    
    yarp::os::Bottle *yLimitsPtr = yLimits.asList();
    if (!yLimitsPtr || yLimitsPtr->size() != 2)
    {
        std::cerr << "Error while reading the Y limits. Wrong dimensions." << std::endl;
        return false;
    }
    limit1 = yLimitsPtr->get(0).asDouble();
    limit2 = yLimitsPtr->get(1).asDouble();
    m_footDimensions[1].first = std::min(limit1, limit2);
    m_footDimensions[1].second = std::max(limit1, limit2);
    
    yarp::os::Value tempValue;
    iDynTree::VectorDynSize tempVector;
    tempValue = mpcOptions.find("CoM_Weight");
    
    //TODO READ DI TUTTI I VALORI
    
    m_configured = true;
    return true;
}

