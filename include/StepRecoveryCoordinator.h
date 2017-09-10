
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
#include "StepRecoveryMPC.h"
#include <yarp/os/all.h>
#include <yarp/sig/all.h>



class StepRecoveryCoordinator : public yarp::os::RFModule {
    yarp::os::BufferedPort< yarp::sig::Vector > controllerPort;
    yarp::os::BufferedPort< yarp::sig::Vector > outputPort;
    StepRecoveryMPC controllerMPC;
    
    double m_period;
    
    iDynTree::VectorDynSize inputBuffer, fL, fL_prev, fR, fR_prev, newGamma, lastGamma;
    
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
