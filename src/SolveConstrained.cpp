#include "SolveConstrained.hpp"
#include <iostream>

void solveConstrainedSymmetric(const Eigen::SparseMatrix<double>& A, const Eigen::MatrixXd& B,
                               const std::vector<int>& constr, const Eigen::MatrixXd& constrValues,
                               Eigen::MatrixXd& X)
{
    assert(constrValues.rows() == constr.size());
    assert(constrValues.cols() == B.cols());
    
    // compute new indizes
    const int n = A.cols();
    std::vector<int> idMap(n, -1);
    
    int cnt = n;
    for(int j = constr.size() - 1; j >= 0; --j)
    {
        if(idMap[constr[j]] != -1) std::cout << "err0" << std::endl;
        
        if(cnt == 0) std::cout << "err1" << std::endl;
        
        idMap[constr[j]] = --cnt;
    }
    
    const int nc = n - cnt;
    const int ni = n - nc;
    
    for(int& i : idMap)
        if(i == -1) i = --cnt;
    
    if(cnt != 0) std::cout << "err2" << std::endl;
    assert(cnt == 0);
    
    // build matrices AII and AIB
    std::vector<Eigen::Triplet<double>> tripII, tripIB;
    
    for(int i = 0; i < A.cols(); ++i)
    {
        const int i2 = idMap[i];
        
        if(i2 < ni)
        {
            for(int j = A.outerIndexPtr()[i]; j < A.outerIndexPtr()[i + 1]; ++j)
            {
                const int j2 = idMap[A.innerIndexPtr()[j]];
                if(j2 < ni) tripII.emplace_back(j2, i2, A.valuePtr()[j]);
            }
        } else
        {
            for(int j = A.outerIndexPtr()[i]; j < A.outerIndexPtr()[i + 1]; ++j)
            {
                const int j2 = idMap[A.innerIndexPtr()[j]];
                if(j2 < ni) tripIB.emplace_back(j2, i2 - ni, A.valuePtr()[j]);
            }
        }
    }
    
    // select inner values of rhs
    Eigen::MatrixXd BI(ni, B.cols());
    for(int j = 0; j < n; ++j)
    {
        if(idMap[j] < ni)
        {
            BI.row(idMap[j]) = B.row(j);
        }
    }
    
    // build, factorize and solve system
    Eigen::SparseMatrix<double> AII(ni, ni), AIB(ni, nc);
     std::cout << "set triplets....";
    AII.setFromTriplets(tripII.begin(), tripII.end());   std::cout << "done" << std::endl;
    AIB.setFromTriplets(tripIB.begin(), tripIB.end());
    
   // Eigen::SparseLU<Eigen::SparseMatrix<double>> chol;
    
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol;
       
        std::cout << "factorize....";
    chol.analyzePattern(AII);
    chol.factorize(AII);
       std::cout << "done" << std::endl;
    Eigen::MatrixXd XI = chol.solve(BI - AIB * constrValues);
    
    // construct final solution
    X.resize(n, B.cols());
    
    for(int j = 0; j < n; ++j)
    {
        if(idMap[j] < ni)
        {
            X.row(j) = XI.row(idMap[j]);
        } else
        {
            X.row(j) = constrValues.row(idMap[j] - ni);
        }
    }
}


