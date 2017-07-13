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

class MatrixBlock{
    iDynTree::VectorDynSize m_values;
    std::vector<unsigned int> m_nzRowIndeces;
    std::vector<unsigned int> m_nzColIndeces;
    std::vector<unsigned int> m_nzRowIndecesOffset;
    std::vector<unsigned int> m_nzColIndecesOffset;
    double m_tol;
    std::pair< unsigned int, unsigned int> m_offsets;
    
    void addOffset();
    
public:
    MatrixBlock();
    ~MatrixBlock();
    
    bool setTolerance(double tol);
    bool setBlock(const iDynTree::MatrixDynSize& block);
    bool setBlock(const iDynTree::MatrixDynSize& block,
                  const std::vector<unsigned int>& nonZerosRowIndeces,
                  const std::vector<unsigned int>& nonZerosColumnIndeces);
    void setOffsets(unsigned int rowOffset, unsigned int colOffset);
    iDynTree::VectorDynSize getValues();
    std::vector<unsigned int> getRowIndeces(bool plusOffset = true);
    std::vector<unsigned int> getColsIndeces(bool plusOffset = true);
    
};

class MPCIpOptSolver : public Ipopt::TNLP {
    
    double m_dT;
    unsigned int m_horizon;
    double m_g;
    
public:
    MPCIpOptSolver();
    ~MPCIpOptSolver();
    
    bool setTimeSettings(double dT, unsigned int horizon);
    
    bool setImpactInstant(unsigned int impact);
    
    bool setGamma0(iDynTree::VectorFixSize<9>& gamma0);
    
    bool setPreviousWrench(iDynTree::VectorFixSize<6>& previousLeftWrench, iDynTree::VectorFixSize<6>& previousRightWrench);
    
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
