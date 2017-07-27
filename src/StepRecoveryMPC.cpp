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
#include <iDynTree/Core/EigenHelpers.h>
#include <Eigen/Core>
#include <cmath>
#include <ctime>

StepRecoveryMPC::StepRecoveryMPC()
:m_configured(false)
,m_reOptimize(false)
,m_dT(0)
,m_horizon(0)
,m_stepDuration(0)
,m_stateL(Standing)
,m_stateR(Standing)
,m_kImpact(0)
,m_robotMass(30)
{
    m_prevL.resize(6);
    m_prevL.zero();
    m_prevR.resize(6);
    m_prevR.zero();
    
    solverPointer = new MPCIpOptSolver();
    loader = IpoptApplicationFactory();
    loader->Options()->SetStringValue("jac_c_constant",   "yes");
    loader->Options()->SetStringValue("jac_d_constant",   "yes");
    loader->Options()->SetStringValue("hessian_constant", "yes");
//    loader->Options()->SetStringValue("print_timing_statistics", "yes");
    //loader->Options()->SetStringValue("derivative_test", "second-order");
    loader->Options()->SetIntegerValue("print_level",0);
    
    loader->Options()->SetStringValue("mehrotra_algorithm", "yes");
   //loader->Options()->SetStringValue("mu_strategy", "monotone");
    
    loader->Options()->SetNumericValue("bound_push", 1e-10);
    loader->Options()->SetNumericValue("bound_frac", 1e-10);
    
    loader->Options()->SetNumericValue("warm_start_bound_push", 1e-10);
    loader->Options()->SetNumericValue("warm_start_bound_frac", 1e-10);
    
    //loader->Options()->SetStringValue("bound_mult_init_method", "mu-based");
    
    //loader->Options()->SetIntegerValue("max_iter",6000);
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
    std::string solvername = mpcOptions.check("solver_name", yarp::os::Value("mumps")).asString();
    loader->Options()->SetStringValue("linear_solver", solvername);
    
    loader->Options()->SetNumericValue("acceptable_tol", mpcOptions.check("accettable_tolerance", yarp::os::Value(1e-6)).asDouble());
    loader->Options()->SetIntegerValue("acceptable_iter", mpcOptions.check("acceptable_iterations", yarp::os::Value(15)).asInt());

    loader->Options()->SetNumericValue("nlp_scaling_max_gradient", mpcOptions.check("nlp_scaling_max_gradient", yarp::os::Value(100.0)).asDouble());

    loader->Options()->SetNumericValue("nlp_scaling_min_value", mpcOptions.check("nlp_scaling_min_value", yarp::os::Value(1e-8)).asDouble());
    
    m_dT = mpcOptions.check("dT", yarp::os::Value(0.01)).asDouble();
    m_horizon = mpcOptions.check("horizon", yarp::os::Value(25)).asInt();
    
    if(! solverPointer->setTimeSettings(m_dT, m_horizon))
        return false;
    
    m_stepDuration = mpcOptions.check("step_duration", yarp::os::Value(0.6)).asDouble();
    m_frictionCoeff = mpcOptions.check("friction_coefficient", yarp::os::Value(0.3333)).asDouble();
    m_torsionalFrictionCoeff = mpcOptions.check("torsional_friction_coeff", yarp::os::Value(0.0133)).asDouble();
    m_edgesPyramid = mpcOptions.check("edges_friction_pyramid", yarp::os::Value(4)).asInt();
    m_fzMin = mpcOptions.check("normalForce_min", yarp::os::Value(10.0)).asDouble();
    m_fzMax = mpcOptions.check("normailForce_max", yarp::os::Value(800.0)).asDouble();
    
    
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
    iDynTree::VectorDynSize pos, vel, ang, gamma;
    
    tempValue = mpcOptions.find("CoM_Weight");
    if(! this->getVectorFromValue(tempValue, pos)){
        std::cerr << "Initialization failed while reading CoM_Weight vector." << std::endl;
        return false;
    }
    
    tempValue = mpcOptions.find("CoMVelocity_Weight");
    if(! this->getVectorFromValue(tempValue, vel)){
        std::cerr << "Initialization failed while reading CoMVelocity_Weight vector." << std::endl;
        return false;
    }
    
    tempValue = mpcOptions.find("AngMom_Weight");
    if(! this->getVectorFromValue(tempValue, ang)){
        std::cerr << "Initialization failed while reading AngMom_Weight vector." << std::endl;
        return false;
    }
    
    gamma.resize(pos.size()+vel.size()+ang.size());
    iDynTree::toEigen(gamma).head(pos.size()) = iDynTree::toEigen(pos);
    iDynTree::toEigen(gamma).segment(pos.size(), vel.size()) = iDynTree::toEigen(vel);
    iDynTree::toEigen(gamma).tail(ang.size()) = iDynTree::toEigen(ang);
    
    if (! solverPointer->setGammaWeight(gamma))
        return false;
    
    tempValue = mpcOptions.find("CoM_ImpactWeight");
    if(! this->getVectorFromValue(tempValue, pos)){
        std::cerr << "Initialization failed while reading CoM_ImpactWeight vector." << std::endl;
        return false;
    }
    
    tempValue = mpcOptions.find("CoMVelocity_ImpactWeight");
    if(! this->getVectorFromValue(tempValue, vel)){
        std::cerr << "Initialization failed while reading CoMVelocity_ImpactWeight vector." << std::endl;
        return false;
    }
    
    tempValue = mpcOptions.find("AngMom_ImpactWeight");
    if(! this->getVectorFromValue(tempValue, ang)){
        std::cerr << "Initialization failed while reading AngMom_ImpactWeight vector." << std::endl;
        return false;
    }
    
    gamma.resize(pos.size()+vel.size()+ang.size());
    iDynTree::toEigen(gamma).head(pos.size()) = iDynTree::toEigen(pos);
    iDynTree::toEigen(gamma).segment(pos.size(), vel.size()) = iDynTree::toEigen(vel);
    iDynTree::toEigen(gamma).tail(ang.size()) = iDynTree::toEigen(ang);
    
    if(! solverPointer->setPostImpactGammaWeight(gamma))
        return false;
    
    iDynTree::VectorDynSize left, right, wrenches;
    tempValue = mpcOptions.find("LeftWrench_Weight");
    if(! this->getVectorFromValue(tempValue, left)){
        std::cerr << "Initialization failed while reading LeftWrench_Weight vector." << std::endl;
        return false;
    }
    
    tempValue = mpcOptions.find("RightWrench_Weight");
    if(! this->getVectorFromValue(tempValue, right)){
        std::cerr << "Initialization failed while reading RightWrench_Weight vector." << std::endl;
        return false;
    }

    wrenches.resize(left.size()+right.size());
    iDynTree::toEigen(wrenches).head(left.size()) = iDynTree::toEigen(left);
    iDynTree::toEigen(wrenches).tail(right.size()) = iDynTree::toEigen(right);
    
    if(! solverPointer->setWrenchsWeight(wrenches))
        return false;
    
    tempValue = mpcOptions.find("LeftWrench_DiffWeight");
    if(! this->getVectorFromValue(tempValue, left)){
        std::cerr << "Initialization failed while reading LeftWrench_DiffWeight vector." << std::endl;
        return false;
    }
    
    tempValue = mpcOptions.find("RightWrench_DiffWeight");
    if(! this->getVectorFromValue(tempValue, right)){
        std::cerr << "Initialization failed while reading RightWrench_DiffWeight vector." << std::endl;
        return false;
    }
    
    wrenches.resize(left.size()+right.size());
    iDynTree::toEigen(wrenches).head(left.size()) = iDynTree::toEigen(left);
    iDynTree::toEigen(wrenches).tail(right.size()) = iDynTree::toEigen(right);
    
    if(! solverPointer->setWrenchDerivativeWeight(wrenches))
        return false;
    
    if(!computeWrenchConstraints())
        return false;
    
    
    if(! solverPointer->prepareProblem()){
        std::cerr << "Error while preparing the optimization problem." << std::endl;
        return false;
    }
    
    Ipopt::ApplicationReturnStatus status = loader->Initialize();
    
    if(status != Ipopt::Solve_Succeeded){
        std::cerr<<"Error during IPOPT solver initialization"<< std::endl;
        return false;
    }
    
    m_configured = true;
    
    if(!dryRun()){
        std::cerr << "Initial test failed." << std::endl;
        return false;
    }

    if(!dryRun()){
        std::cerr << "Initial test #2 failed." << std::endl;
        return false;
    }
    
    if(!dryRun()){
        std::cerr << "Initial test #3 failed." << std::endl;
        return false;
    }
    return true;
}

bool StepRecoveryMPC::computeWrenchConstraints()
{
    if(m_frictionCoeff <= 0){
        std::cerr << "The friction coefficient is expected to be positive." <<std::endl;
        return false;
    }
    if(m_torsionalFrictionCoeff <= 0){
        std::cerr << "The torsional friction coefficient is expected to be positive." <<std::endl;
        return false;
    }
    
    if(m_edgesPyramid < 3){
        std::cerr << "The number of edges of the friction pyramid should be greater than 4" <<std::endl;
        return false;
    }
    
    if(m_fzMin <= 0){
        std::cerr << "The minimum normal force is expected to be positive." <<std::endl;
        return false;
    }
    
    if(m_fzMax <= 0){
        std::cerr << "The maximum normal force is expected to be positive." <<std::endl;
        return false;
    }
    
    if(m_fzMax < m_fzMin){
        std::cerr << "The maximum normal force is expected to be greater than the minimum." <<std::endl;
        return false;
    }
    
    m_wrenchConstrMatrix.resize(m_edgesPyramid+8, 6); //There are m_edgesPyramid constraints for friction pyramid, 4 for the CoP, 2 for the torsional friction, 2 for the bounds on the normal force.
    m_wrenchConstrMatrix.zero();
    
    m_wrenchConstrVector.resize(m_wrenchConstrMatrix.rows());
    m_wrenchConstrVector.zero();
    m_afterImpactwrenchConstrVector.resize(m_wrenchConstrMatrix.rows());
    m_afterImpactwrenchConstrVector.zero();
    
    double angle = 0;
    double sectorAngle = 2*M_PI/m_edgesPyramid;
    double angularCoeff, offset;
    int inequalityFactor = 1;
    std::pair <double, double> point1, point2;
    for(int edge = 0; edge < m_edgesPyramid; ++edge){
        point1.first = std::cos(angle);
        point1.second = std::sin(angle);
        
        if(edge == m_edgesPyramid-1){
            point2.first = std::cos(0);
            point2.second = std::sin(0);
        }
        else{
            point2.first = std::cos(angle+sectorAngle);
            point2.second = std::sin(angle+sectorAngle);
        }
        
        angularCoeff = (point2.second - point1.second)/(point2.first - point1.first);
        offset = point1.second - angularCoeff*point1.first;
        
        inequalityFactor = 1;
        if((angle > M_PI)|| ((angle+sectorAngle) > M_PI))
            inequalityFactor = -1;
        
        m_wrenchConstrMatrix(edge, 0) = -inequalityFactor*angularCoeff;
        m_wrenchConstrMatrix(edge, 1) = inequalityFactor;
        m_wrenchConstrMatrix(edge, 2) = -inequalityFactor*offset*m_frictionCoeff;
        
        angle += sectorAngle;
    }
    unsigned int row = m_edgesPyramid;
    iDynTree::toEigen(m_wrenchConstrMatrix).row(row) << 0, 0, -m_torsionalFrictionCoeff, 0, 0, 1;
    row++;
    iDynTree::toEigen(m_wrenchConstrMatrix).row(row) << 0, 0, -m_torsionalFrictionCoeff, 0, 0, -1;
    row++;
    iDynTree::toEigen(m_wrenchConstrMatrix).row(row) << 0, 0, m_footDimensions[0].first, 0, 1, 0;
    row++;
    iDynTree::toEigen(m_wrenchConstrMatrix).row(row) << 0, 0, -m_footDimensions[0].second, 0, -1, 0;
    row++;
    iDynTree::toEigen(m_wrenchConstrMatrix).row(row) << 0, 0, m_footDimensions[1].first, -1, 0, 0;
    row++;
    iDynTree::toEigen(m_wrenchConstrMatrix).row(row) << 0, 0, -m_footDimensions[1].second, 1, 0, 0;
    row++;
    iDynTree::toEigen(m_wrenchConstrMatrix).row(row) << 0, 0, -1, 0, 0, 0;
    row++;
    iDynTree::toEigen(m_wrenchConstrMatrix).row(row) << 0, 0, 1, 0, 0, 0;
    
    m_wrenchConstrVector(m_wrenchConstrVector.size() - 2) = 0*1e-6;
    m_wrenchConstrVector(m_wrenchConstrVector.size() - 1) = 0*1e-6;
    
    m_afterImpactwrenchConstrVector(m_afterImpactwrenchConstrVector.size() - 2) = -m_fzMin;
    m_afterImpactwrenchConstrVector(m_afterImpactwrenchConstrVector.size() - 1) = m_fzMax;
    
    solverPointer->setWrenchConstraints(m_wrenchConstrMatrix, m_wrenchConstrVector, m_afterImpactwrenchConstrVector);
    
    return true;
}

bool StepRecoveryMPC::getControllerData(const iDynTree::VectorDynSize& controllerData)
{
    if(controllerData.size() != m_inputData.expectedControllerDim){
        std::cerr << "The dimension of the data coming from the controller doens't match the expected dimension. Input = "<< controllerData.size() << ". Expected = " << m_inputData.expectedControllerDim << std::endl;
        return false;
    }
    Eigen::Map<const Eigen::VectorXd> controllerData_map(controllerData.data(), m_inputData.expectedControllerDim);
    
    iDynTree::Position pos;
    iDynTree::Vector4 quat;
    iDynTree::Rotation rot;
    unsigned int row = 0;
    
    iDynTree::toEigen(pos) = controllerData_map.head<3>();
    row += 3;
    iDynTree::toEigen(quat) = controllerData_map.segment<4>(row);
    row +=4;
    rot.fromQuaternion(quat);
    m_inputData.leftTransform.setPosition(pos);
    m_inputData.leftTransform.setRotation(rot);
    
    
    iDynTree::toEigen(pos) = controllerData_map.segment<3>(row);
    row += 3;
    iDynTree::toEigen(quat) = controllerData_map.segment<4>(row);
    row += 4;
    rot.fromQuaternion(quat);
    m_inputData.rightTransform.setPosition(pos);
    m_inputData.rightTransform.setRotation(rot);
    
    iDynTree::toEigen(m_inputData.gamma) = controllerData_map.segment<9>(row);
    row += 9;
    
    m_inputData.kImpact = std::round(controllerData(row));
    row++;
    
    m_inputData.comZDes = controllerData(row);
    row++;
    
    m_inputData.robotMass = controllerData(row);
    row++;
    
    m_inputData.controllerState = std::round(controllerData(row));
    row++;
    
    if(m_inputData.controllerState == 14){
        m_stateL = Standing;
        m_stateR = Swinging;
    }
    
    if(m_inputData.controllerState == 15){
        m_stateL = Standing;
        m_stateR = Standing;
    }
    
    if(m_inputData.controllerState == 16){
        m_stateL = Standing;
        m_stateR = Floating;
    }

    return true;
}

bool StepRecoveryMPC::setFeetTransform()
{
    return solverPointer->setFeetTransforms(m_inputData.leftTransform, m_inputData.rightTransform);
}


bool StepRecoveryMPC::getGamma()
{
    m_gamma0 = m_inputData.gamma;
    m_robotMass = m_inputData.robotMass;
    return solverPointer->setGamma0(m_gamma0) && solverPointer->setRobotMass(m_robotMass);
}

bool StepRecoveryMPC::computeImpactInstant()
{
    m_kImpact = m_inputData.kImpact;
    return solverPointer->setImpactInstant(m_kImpact, true); //TODO BETTER HANDLING OF THE SWING FOOT, MAYBE USING THE FootState
}

bool StepRecoveryMPC::setDesiredCoM()
{
    iDynTree::Position desiredCoM;
    Eigen::Map <Eigen::VectorXd> desiredCoM_map(desiredCoM.data(), 3);
    if(m_stateL == Standing){
        if(m_stateR == Standing){
            desiredCoM = m_inputData.leftTransform.getPosition() + m_inputData.rightTransform.getPosition();
            desiredCoM_map = desiredCoM_map/2;
            desiredCoM(2) = m_inputData.comZDes;
        }
        else{
            if(m_stateR == Floating){
                desiredCoM = m_inputData.leftTransform.getPosition();
                desiredCoM(2) = m_inputData.comZDes;
            }
            else{
                desiredCoM = m_inputData.leftTransform.getPosition() + m_inputData.rightTransform.getPosition();
                desiredCoM_map = desiredCoM_map/2;
                desiredCoM(2) = m_inputData.comZDes;
            }
        }
    }
    else{
        if(m_stateR == Standing){
            if(m_stateL == Swinging){
                desiredCoM = m_inputData.leftTransform.getPosition() + m_inputData.rightTransform.getPosition();
                desiredCoM_map = desiredCoM_map/2;
                desiredCoM(2) = m_inputData.comZDes;
            }
            else{
                desiredCoM = m_inputData.rightTransform.getPosition();
                desiredCoM(2) = m_inputData.comZDes;
            }
        }
        else{
            desiredCoM = m_inputData.rightTransform.getPosition();
            desiredCoM(2) = m_inputData.comZDes;
        }
    }
    
    return solverPointer->setDesiredCOMPosition(desiredCoM);
}

bool StepRecoveryMPC::setPreviousWrench()
{
    return solverPointer->setPreviousWrench(m_prevL, m_prevR);
}


bool StepRecoveryMPC::solve(const iDynTree::VectorDynSize& controllerData, iDynTree::VectorDynSize& fL, iDynTree::VectorDynSize& fR, iDynTree::VectorDynSize& lastGamma)
{
    if(!m_configured){
        std::cerr <<"First you have to call the configure method." << std::endl;
        return false;
    }
    
    if(!getControllerData(controllerData)){
        std::cerr << "Error while reading the data from the controller." << std::endl;
        return false;
    }
    
    if(!setFeetTransform()){
        std::cerr << "Error while setting the feet transformations." << std::endl;
        return false;
    }
    
    if(!getGamma()){
        std::cerr << "Error while setting the feedback." << std::endl;
        return false;
    }
    
    if(!computeImpactInstant()){
        std::cerr << "Error while setting the impact instant." << std::endl;
        return false;
    }
    
    if(!setDesiredCoM()){
        std::cerr << "Error while setting the desired CoM." << std::endl;
        return false;
    }
    
    if(!setPreviousWrench()){
        std::cerr << "Error while setting the previous wrench." << std::endl;
        return false;
    }

    if(! solverPointer->updateProblem()){
        std::cerr << "Error while updating the optimization problem." << std::endl;
        return false;
    }
    
    int exitCode;
    
    clock_t begin, end;

    if(m_reOptimize){
        loader->Options()->SetStringValue("warm_start_init_point", "yes");
        loader->Options()->SetNumericValue("warm_start_bound_frac", 1e-6);
        loader->Options()->SetNumericValue("warm_start_bound_push", 1e-6);
        loader->Options()->SetNumericValue("warm_start_mult_bound_push", 1e-6);
        loader->Options()->SetNumericValue("warm_start_slack_bound_frac", 1e-6);
        loader->Options()->SetNumericValue("warm_start_slack_bound_push", 1e-6);
    }

//    if(m_reOptimize){
//        //std::cerr << "ReOptimizeTNLP!!" << std::endl;
//        begin = clock();
//        loader->ReOptimizeTNLP(solverPointer);
//        end = clock();
//        exitCode = solverPointer->getSolution(fL, fR, lastGamma);
//        if(exitCode < 0){
//            std::cerr << "Optimization problem failed!" << std::endl;
//            return false;
//        }
//        m_prevL = fL;
//        m_prevR = fR;
//    }
//    else{
        begin = clock();
        loader->OptimizeTNLP(solverPointer);
        end = clock();
        exitCode = solverPointer->getSolution(fL, fR, lastGamma);
        if(exitCode < 0){
            std::cerr << "Optimization problem failed!" << std::endl;
            return false;
        }
        m_prevL = fL;
        m_prevR = fR;
        m_reOptimize = true;
//    }
    std::cerr << "Solved in: " << double(end - begin) / CLOCKS_PER_SEC << "sec." << std::endl;
    return true;
}

int StepRecoveryMPC::dryRun()
{
    Eigen::Vector3d pL, pR;
    iDynTree::Vector4 quatL, quatR;
    pL << 0, 0, 0;
    quatL = iDynTree::Rotation::Identity().asQuaternion();
    pR << 0.05, -0.3, 0;
    quatR = iDynTree::Rotation::RotZ(M_PI/6).asQuaternion();
    
    Eigen::VectorXd gamma0(9);
    gamma0.setZero();
    gamma0(2) = 0.5;
    
    double k_impact = 2;
    
    double comZDes = 0.5;
    
    double mass = 30.0;
    
    double state = 14;
    
    iDynTree::VectorDynSize dummyController;
    dummyController.resize(27);
    iDynTree::toEigen(dummyController) << pL, iDynTree::toEigen(quatL), pR, iDynTree::toEigen(quatR), gamma0, k_impact, comZDes, mass, state;
    m_prevL(2) = 150;
    m_prevR(2) = 150;
    
    iDynTree::VectorDynSize dummyVector;
    return solve(dummyController, m_prevL, m_prevR, dummyVector);
}


void StepRecoveryMPC::setVerbosity(unsigned int level)
{
    loader->Options()->SetIntegerValue("print_level",level);
}

