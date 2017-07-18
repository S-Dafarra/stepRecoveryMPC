
/******************************************************************************
 *                                                                            *
 * Copyright (C) 2017 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @file StepRecoveryCoordinator.h
 * @authors: Stefano Dafarra <stefano.dafarra@iit.it>
 */

#ifndef STEPRECOVERYCOORDINATOR_H
#define STEPRECOVERYCOORDINATOR_H

#include <yarp/os/RFModule.h>
#include "MPCIpOptSolver.h"

class StepRecoveryCoordinator : public yarp::os::RFModule {
    
public:
    StepRecoveryCoordinator();
    ~StepRecoveryCoordinator();
    // RFModule methods
    virtual bool configure(yarp::os::ResourceFinder &rf);
    virtual double getPeriod();
    virtual bool updateModule();
    virtual bool close();
    
};



#endif
