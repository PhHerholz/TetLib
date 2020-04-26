#include <iostream>
#include "IndexedTetMesh.hpp"
#include "VectorOperations.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unordered_map>
#include <fstream>
#include <array>

IndexedTetMesh::IndexedTetMesh()
{
    
}


IndexedTetMesh::~IndexedTetMesh()
{
    
}

void IndexedTetMesh::write(std::string fname)
{
    std::ofstream file(fname);
    file << vertices.size() << " " << tets.size() << "\n";
    
    for(auto& v : vertices)
        file << v << "\n";
   
    for(auto& t : tets)
        file << t << "\n";
      
    file.close();
}

void IndexedTetMesh::read(std::string fname)
{
    std::ifstream file(fname);
    int nv, nt;
    
    file >> nv;
    file >> nt;
    
    vertices.resize(nv);
    tets.resize(nt);
    
    for(int i = 0; i < nv; ++i)
    {
        file >> vertices[i][0];
        file >> vertices[i][1];
        file >> vertices[i][2];
    }

    for(int i = 0; i < nt; ++i)
    {
        file >> tets[i][0];
        file >> tets[i][1];
        file >> tets[i][2];
        file >> tets[i][3];
    }
    
    file.close();
}

// vert-vert neighbours
void
IndexedTetMesh::buildVertexNeighbours()
{
    const int NV = vertices.size();
   
    vertexNeighbours.resize(NV);
    
    for(auto& t : tets)
    {
        for(int i = 0; i < 4; ++i)
        {
            for(int j = 0; j < 4; ++j)
            {
                if(i != j)
                {
                    vertexNeighbours[t[i]].push_back(t[j]);
                }
            }
        }
    }
    
    for(auto& nbh : vertexNeighbours)
    {
        sort(nbh.begin(), nbh.end());
        nbh.erase(unique(nbh.begin(), nbh.end()) , nbh.end());
    }
}

std::vector<int>
IndexedTetMesh::sortRing(const std::vector<std::array<int, 2>>& ring)
{
    const int nr = ring.size();
   
    for(int i = 0; i < nr; ++i)
    {
        flag[ring[i][1]] = i;
    }
    
    int first = 0;
    
    for(int i = 0; i < nr; ++i)
    {
        if(flag[ring[i][0]] == -1)
        {
            first = i;
        }
    }
    
    for(int i = 0; i < nr; ++i)
    {
        flag[ring[i][1]] = -1;
    }
    
    for(int i = 0; i < nr; ++i)
    {
        flag[ring[i][0]] = i;
    }
    
    std::vector<int> ret;
    ret.reserve(nr);
    
    int i = first;
    int j;
    do
    {
        ret.push_back(ring[i][0]);
        j = ring[i][1];
        i = flag[ring[i][1]];
    } while(i != first && i != -1);
    
    ret.push_back(j);
    
    for(int i = 0; i < nr; ++i)
    {
        flag[ring[i][0]] = -1;
        flag[ring[i][1]] = -1;
    }
    
    return ret;
}


std::vector<char>
IndexedTetMesh::boundaryVertexFlag()
{
    std::vector<char> ret(vertices.size(), 0);
    
    if(tetNeighbours.size() != vertices.size())
        buildTetNeighbours();
        
    for(int i = 0; i < tets.size(); ++i)
    {
        for(int j = 0; j < 4; ++j)
        {
            if(tetNeighbours[i][j] == -1)
            {
                for(int k = 1; k < 4; ++k)
                    ret[ tets[i][(j + k) % 4] ] = 1;
            }
        }
    }
    
    return ret;
}


void
IndexedTetMesh::buildEdgeNeighours()
{
    using namespace std;
    flag = vector<int>(vertices.size(), -1);
    unordered_map<array<unsigned int, 2>, vector<array<int, 2>>> edgeNeighbours0;
    
    const int flipEdge[4][4][2]
    {
        {{-1, -1},
        {2, 3},
        {3, 1},
        {1, 2}},
        
        {{3, 2},
        {-1, -1},
        {0, 3},
        {2, 0}},
        
        {{1, 3},
        {3, 0},
        {-1, -1},
        {0, 1}},
        
        {{2, 1},
        {0, 2},
        {1, 0},
        {-1, -1}}
    };
    
    int cnt = 0;
    for(auto& t : tets)
    {
        for(int i = 0; i < 4; ++i)
        {
            for(int j = 0; j < 4; ++j)
            {
                if(t[i] < t[j])
                {
                    edgeNeighbours0[array<unsigned int, 2>{{t[i], t[j]}}].push_back(array<int, 2>{{(int)t[flipEdge[i][j][0]], (int)t[flipEdge[i][j][1]]}});
                }
            }
        }
        
        ++cnt;
    }
    
    for(const auto& x : edgeNeighbours0)
    {
        edgeCellNeighbours[x.first] = sortRing(x.second);
    }
}


void
IndexedTetMesh::buildTetNeighbours()
{
    using namespace std;
    
    typedef unsigned long long int Index;
        
    const int NT = tets.size();
    const Index nv = vertices.size();
    const Index  nv2 = nv * nv;
    
    auto faceToIndex = [=](const array<unsigned int, 4>& cell, const int facet)
    {
        unsigned f[3]{cell[(facet+1)%4], cell[(facet+2)%4], cell[(facet+3)%4]};
        
        if (f[0] > f[2])
           swap(f[0], f[2]);

        if (f[0] > f[1])
           swap(f[0], f[1]);

        if (f[1] > f[2])
           swap(f[1], f[2]);
        
        return f[0] + nv * f[1] + nv2 * f[2];
    };
    
    
    unordered_map<Index, std::array<int,2>> faceCell;
    
    tetNeighbours.clear();
    tetNeighbours.resize(NT, std::array<int, 4>{-1,-1,-1,-1});
    
    for(int i = 0; i < NT; ++i)
    {
        for(int j = 0; j < 4; ++j)
        {
            auto id = faceToIndex(tets[i], j);
            auto it = faceCell.find(id);
            
            if(it != faceCell.end())
            {
                tetNeighbours[i][j] = it->second[0];
                tetNeighbours[it->second[0]][it->second[1]] = i;
                faceCell.erase(it);
            } else
            {
                faceCell[id] = std::array<int, 2>{i, j};
            }
        }
    }
}



void IndexedTetMesh::circumcenter(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, Eigen::Vector3d& cc)
{
    const double l[3]{
        (b - c).squaredNorm(),
        (a - c).squaredNorm(),
        (a - b).squaredNorm()
    };
    
    const double ba[3]{l[0] * (l[1] + l[2] - l[0]), l[1] * (l[2] + l[0] - l[1]), l[2] * (l[0] + l[1] - l[2])};
    const double sum = ba[0] + ba[1] + ba[2];
    
    cc = (ba[0] / sum) * a + (ba[1] / sum) * b + (ba[2] / sum) * c;
}

void IndexedTetMesh::circumcenter(const Eigen::Matrix<double, 4, 3>& t, Eigen::Vector3d& c)
{
    Eigen::Matrix3d A;
    Eigen::Vector3d b;
    
    const double n0 = t.row(0).squaredNorm();
    
    for(int k = 0; k < 3; ++k)
    {
        A.row(k) = t.row(k + 1) - t.row(0);
        b(k) = t.row(k + 1).squaredNorm() - n0;
    }
    
    c = 0.5 * A.fullPivHouseholderQr().solve(b);
}

double IndexedTetMesh::volume(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d)
{
    Eigen::Matrix3d A;
    A.col(0) = b - a;
    A.col(1) = c - a;
    A.col(2) = d - a;
    
    return A.determinant() / 6.;
}

void IndexedTetMesh::dualLaplace(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T, Eigen::SparseMatrix<double>& L, Eigen::SparseMatrix<double>& M)
{
    const size_t nt = T.rows();
    const size_t nv = V.rows();
    
    const int turn[4][4]
    {
        {-1, 2, 3, 1},
        {3, -1, 0, 2},
        {1, 3, -1, 0},
        {2, 0, 1, -1}
    };

    auto getTet = [&](const int i, Eigen::Matrix<double, 4, 3>& t)
    {
        for(int k = 0; k < 4; ++k)
        {
            t.row(k) = V.row(T(i, k));
        }
    };
    
    std::vector<Eigen::Triplet<double>> tripL, tripM;
    
    Eigen::Vector3d cc;
    Eigen::Matrix<double, 4, 3> t;
    
    for(int k = 0; k < nt; ++k)
    {
        getTet(k, t);
        circumcenter(t, cc);
        
        for(int i = 0; i < 4; ++i)
        {
            for(int j = 0; j < 4; ++j)
            {
                if(i != j)
                {
                    Eigen::Vector3d cf;
                    circumcenter(t.row(i), t.row(j), t.row(turn[i][j]), cf);
                    
                    const Eigen::Vector3d ce = 0.5 * (t.row(i) + t.row(j));
                    
                    const double vol = volume(t.row(i), ce, cf, cc);
                    const double wij = 6. * vol / (t.row(i) - t.row(j)).squaredNorm();
                    
                    tripL.emplace_back(T(k, i), T(k, j), wij);
                    tripL.emplace_back(T(k, j), T(k, i), wij);
                    
                    tripL.emplace_back(T(k, i), T(k, i), -wij);
                    tripL.emplace_back(T(k, j), T(k, j), -wij);
                    
                    tripM.emplace_back(T(k, i), T(k, i), vol);
                    tripM.emplace_back(T(k, j), T(k, j), vol);
                }
            }
        }
    }
    
    L.resize(nv, nv);
    M.resize(nv, nv);
    
    L.setFromTriplets(tripL.begin(), tripL.end());
    M.setFromTriplets(tripM.begin(), tripM.end());
}


void IndexedTetMesh::dualLaplace(Eigen::SparseMatrix<double>& L, Eigen::SparseMatrix<double>& M) {

	int nv = vertices.size();
	int nt = tets.size();

	std::cout << "... calc laplace for " << nv << " vertices and " << nt << " tets" << std::endl;

	Eigen::MatrixXd V(nv, 3);
	Eigen::MatrixXi T(nt, 4); 

	for (int i=0; i<nv; ++i) {
		V(i,0) = vertices[i][0];
		V(i,1) = vertices[i][1];
		V(i,2) = vertices[i][2];
	}

	for (int i=0; i<nt; ++i) {
		T(i, 0) = tets[i][0];
		T(i, 1) = tets[i][1];
		T(i, 2) = tets[i][2];
		T(i, 3) = tets[i][3];
	}

	dualLaplace(V,T,L,M);
}
