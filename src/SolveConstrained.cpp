#include "SolveConstrained.hpp"

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
    for(int j = constr.size() - 1; j >= 0; --j) idMap[constr[j]] = --cnt;
    
    const int nc = n - cnt;
    const int ni = n - nc;
    
    for(int& i : idMap)
        if(i == -1) i = --cnt;
    
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
    AII.setFromTriplets(tripII.begin(), tripII.end());
    AIB.setFromTriplets(tripIB.begin(), tripIB.end());
    
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol;
    
    chol.analyzePattern(AII);
    chol.factorize(AII);
    
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


