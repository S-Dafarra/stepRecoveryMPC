/******************************************************************************
 *                                                                            *
 * Copyright (C) 2017 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @file MPCIpOptSolver.h
 * @authors: Stefano Dafarra <stefano.dafarra@iit.it>
 */

#ifndef MPCIPOPTSOLVER_H
#define MPCIPOPTSOLVER_H

#define HAVE_STDDEF_H
#define HAVE_CSTDDEF
#include <IpTNLP.hpp>
#undef HAVE_STDDEF_H
#undef HAVE_CSTDDEF //workaroud for missing libraries

#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <iDynTree/Core/Transform.h>
#include <iDynTree/Core/SparseMatrix.h>
#include <memory>

typedef struct {
    std::shared_ptr<iDynTree::SparseMatrix> blockPtr;
    unsigned int rowOffset;
    unsigned int colOffset;
} MatrixBlock;

class MPCIpOptSolver : public Ipopt::TNLP {
    
    double m_dT;
    unsigned int m_horizon;
    double m_g;
    unsigned int m_impact;
    bool m_rightFootStep;
    double m_mass;
    
    iDynTree::VectorFixSize<9> m_gamma0;
    iDynTree::VectorDynSize m_fLPrev, m_fRPrev;
    
    iDynTree::MatrixDynSize m_wrenchA, m_wrenchAl, m_wrenchAr;
    std::shared_ptr<iDynTree::SparseMatrix> m_wrenchAlSparsePtr, m_wrenchArSparsePtr;
    iDynTree::VectorDynSize m_wrenchb, m_wrenchbImpact;
    
    iDynTree::Transform m_wHl, m_wHr;
    
    iDynTree::VectorFixSize<9> m_desiredGamma;
    
    iDynTree::VectorDynSize m_gammaWeight, m_gammaWeightImpact;
    iDynTree::VectorDynSize m_wrenchWeight, m_derivativeWrenchWeight;
    
    iDynTree::MatrixDynSize m_EvGamma;
    std::shared_ptr<iDynTree::SparseMatrix> m_EvGammaSparsePtr;
    iDynTree::MatrixDynSize m_FGamma;
    std::shared_ptr<iDynTree::SparseMatrix> m_FGammaSparsePtr;
    iDynTree::VectorFixSize<9> m_bias;
    
    std::vector<MatrixBlock> m_modelConstraintsJacobian;
    std::vector<MatrixBlock> m_wrenchConstraintJacobian;
    
    iDynTree::MatrixDynSize m_skewBuffer;
    std::shared_ptr<iDynTree::SparseMatrix> m_minusIdentityPtr;
    iDynTree::MatrixDynSize m_wrenchTransform;
    iDynTree::SparseMatrix m_costHessian, m_gammaWeightHessian, m_gammaWeightImpactHessian, m_wrenchWeightHessian, m_derivativeWrenchWeightHessian, m_negativeDerWrenchHessian;
    
    iDynTree::VectorDynSize m_previousSolution;
    int m_exitCode;
    
    bool computeModelMatrices();
    bool computeModelBias();
    bool computeModelConstraintsJacobian();
    bool computeWrenchConstraints();
    bool computeWrenchConstraintsJacobian();
    bool computeSingleCostHessian();
    bool computeCostHessian();
    
    
    void printJacobian(std::string header, std::vector<MatrixBlock>& input, int rows, int cols);
    
public:
    MPCIpOptSolver();
    
    ~MPCIpOptSolver();
    
    bool setTimeSettings(double dT, unsigned int horizon);
    
    bool setRobotMass(const double mass);
    
    bool setImpactInstant(unsigned int impact, bool rightFootStepping = true);
    
    bool setGamma0(const iDynTree::VectorFixSize<9>& gamma0);
    
    bool setPreviousWrench(const iDynTree::VectorDynSize& previousLeftWrench, const iDynTree::VectorDynSize& previousRightWrench);
    
    bool setWrenchConstraints(const iDynTree::MatrixDynSize& wrenchConstraintsMatrix, const iDynTree::VectorDynSize& constraintsBounds, const iDynTree::VectorDynSize& afterImpactConstraintsBounds);
    
    bool setFeetTransforms(const iDynTree::Transform& w_H_l, const iDynTree::Transform& w_H_r);
    
    bool setDesiredCOMPosition(const iDynTree::Position& desiredCOM);
    
    bool setGammaWeight(const iDynTree::VectorDynSize& gammaWeight);
    
    bool setPostImpactGammaWeight(const iDynTree::VectorDynSize& gammaImpactWeight);
    
    bool setWrenchsWeight(const iDynTree::VectorDynSize& wrenchWeight);
    
    bool setWrenchDerivativeWeight(const iDynTree::VectorDynSize& derivativeWrenchWeight);
    
    bool prepareProblem();
    
    bool updateProblem();
    
    int getSolution(iDynTree::VectorDynSize& fL, iDynTree::VectorDynSize& fR, iDynTree::VectorDynSize& lastGamma);
    
    // IPOPT methods redefinition
    
    bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                      Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style);
    
    bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                         Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u);
    
    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                            bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                            Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda);
    
    bool eval_f(Ipopt::Index n, const Ipopt::Number* x,
                bool new_x, Ipopt::Number& obj_value);
    
    bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                     Ipopt::Number* grad_f);
    
    bool eval_g(Ipopt::Index n, const Ipopt::Number* x,
                bool new_x, Ipopt::Index m, Ipopt::Number* g);
    
    bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                    Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
                    Ipopt::Index *jCol, Ipopt::Number* values);
    
    bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                Ipopt::Index* jCol, Ipopt::Number* values);
    
    void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                           const Ipopt::Number* x, const Ipopt::Number* z_L,
                           const Ipopt::Number* z_U, Ipopt::Index m, const Ipopt::Number* g,
                           const Ipopt::Number* lambda, Ipopt::Number obj_value,
                           const Ipopt::IpoptData* ip_data,
                           Ipopt::IpoptCalculatedQuantities* ip_cq);

};
#endif
