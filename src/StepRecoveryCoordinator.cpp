/******************************************************************************
 *                                                                            *
 * Copyright (C) 2017 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @file StepRecoveryCoordinator.cpp
 * @authors: Stefano Dafarra <stefano.dafarra@iit.it>
 */

#include "StepRecoveryCoordinator.h"
#include <iDynTree/yarp/YARPConversions.h>
#include <iDynTree/yarp/YARPEigenConversions.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <Eigen/Core>
#include <yarp/os/Network.h>
#include <yarp/os/LogStream.h>

int main(int argc, char * argv[])
{
    // YARP setting
    yarp::os::Network yarp;
    if (!yarp::os::Network::checkNetwork(5.0))
    {
        yError() << " YARP server not available!";
        return EXIT_FAILURE;
    }
    
    // Configure ResourceFinder
    yarp::os::ResourceFinder &rf = yarp::os::ResourceFinder::getResourceFinderSingleton();
    rf.setDefaultContext("StepRecoveryMPC");
    rf.setDefaultConfigFile("mpc.ini");
    rf.configure(argc, argv);
    
    // Configure the module
    StepRecoveryCoordinator module;
    
    return module.runModule(rf);
}


StepRecoveryCoordinator::StepRecoveryCoordinator(){
    
    inputBuffer.resize(EXPECTED_INPUT_DIM);
    fL.resize(6);
    fL.zero();
    fL_prev.resize(6);
    fL_prev.zero();
    fR.resize(6);
    fR.zero();
    fR_prev.resize(6);
    fR_prev.zero();
    lastGamma.resize(9);
    lastGamma.zero();
}

StepRecoveryCoordinator::~StepRecoveryCoordinator()
{
}

bool StepRecoveryCoordinator::configure(yarp::os::ResourceFinder& rf)
{
    std::string name = rf.check("name", yarp::os::Value("stepRecoveryMPC")).asString();
    std::string robot = rf.check("robot", yarp::os::Value("icubSim")).asString();
    int verbosityLevel = rf.check("solver_verbosity", yarp::os::Value(0)).asInt();
    
    controllerMPC.setVerbosity(verbosityLevel);
    
    yarp::os::Bottle &MPCSection = rf.findGroup("MPC");
    
    if(!controllerMPC.configure(MPCSection))
        return false;
    
    if(!controllerPort.open("/"+name+":in")){
        std::cerr << "Failed to open the input port." << std::endl;
        return false;
    }
    controllerPort.setInputMode(true);
    
    if(!outputPort.open("/"+name+":out")){
        std::cerr << "Failed to open the output port." << std::endl;
        return false;
    }
    outputPort.setOutputMode(true);
    
    return true;
}

bool StepRecoveryCoordinator::updateModule()
{
    yarp::sig::Vector *inputData = controllerPort.read(true);
    
    if(inputData == NULL)
        return false;

    if(!iDynTree::toiDynTree(*inputData, inputBuffer)){
        std::cerr << "Wrong dimension on the input data" << std::endl;
        return false;
    }
    
    bool solved = controllerMPC.solve(inputBuffer, fL, fR, lastGamma);
    
    yarp::sig::Vector& outputBuffer = outputPort.prepare();
    outputBuffer.resize(21);
    
    if(!solved){
        std::cerr << "Failed to find a solution" << std::endl;
        iDynTree::toEigen(outputBuffer) << iDynTree::toEigen(fL_prev), iDynTree::toEigen(fR_prev), iDynTree::toEigen(lastGamma);
    }
    else{
        iDynTree::toEigen(outputBuffer) << iDynTree::toEigen(fL), iDynTree::toEigen(fR), iDynTree::toEigen(lastGamma);
        fL_prev = fL;
        fR_prev = fR;
    }
    outputPort.write();
    
    return true;
}

double StepRecoveryCoordinator::getPeriod()
{
    return 10.0;
}

bool StepRecoveryCoordinator::close()
{
    controllerPort.interrupt();
    controllerPort.close();
    outputPort.interrupt();
    outputPort.close();
    return true;
}
