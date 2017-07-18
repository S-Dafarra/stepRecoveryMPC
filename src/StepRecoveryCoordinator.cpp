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
}

StepRecoveryCoordinator::~StepRecoveryCoordinator()
{
}

bool StepRecoveryCoordinator::configure(yarp::os::ResourceFinder& rf)
{
}

bool StepRecoveryCoordinator::updateModule()
{
}

double StepRecoveryCoordinator::getPeriod()
{
}

bool StepRecoveryCoordinator::close()
{
}

