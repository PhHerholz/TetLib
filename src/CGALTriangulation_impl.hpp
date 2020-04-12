#include <CGAL/Kernel/global_functions.h>
#include <Eigen/Sparse>
#include <fstream>
#include <unordered_map>


template<class TKernel>
CGALTriangulation<TKernel>::CGALTriangulation()
{
}

template<class TKernel>
void
CGALTriangulation<TKernel>::setIndizes()
{
    int cnt = 0;
    for(auto it = mesh.finite_vertices_begin(); it != mesh.finite_vertices_end(); ++it)
    {
        it->info() = cnt++;
    }
    
    mesh.infinite_vertex()->info() = -1;
    
    cnt = 0;
    for(auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it)
    {
        if(mesh.is_infinite(it)) it->info() = -1;
        else it->info() = cnt++;
    }
}

template<class TKernel>
CGALTriangulation<TKernel>::CGALTriangulation(const std::vector<typename TKernel::Point_3>& pts)
{
    Delaunay del;
    int cnt = 0;
   
    del.infinite_vertex()->info() = -1;
    
    for(auto& p : pts)
    {
        del.insert(p)->info() = cnt++;
    }
    
    cnt = 0;
    for(auto it = del.cells_begin(); it != del.cells_end(); ++it)
    {
        if(del.is_infinite(it)) it->info() = -1;
        else it->info() = cnt++;
    }
    
    mesh = del;
}

template<class TKernel>
CGALTriangulation<TKernel>::CGALTriangulation(Triangulation&& _mesh)
: mesh(_mesh)
{
    
}

template<class TKernel>
CGALTriangulation<TKernel>::~CGALTriangulation()
{
    
}


template<class TKernel>
IndexedTetMesh
CGALTriangulation<TKernel>::toIndexed() const
{
    IndexedTetMesh ret;
    
    for(auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it)
        if(it->info() != -1)
            ret.vertices.push_back(std::array<double, 3>{it->point().x(), it->point().y(), it->point().z()});
                                   
    for(auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it)
        if(it->info() != -1)
            ret.tets.push_back(std::array<unsigned int, 4>{it->vertex(0)->info(), it->vertex(1)->info(), it->vertex(2)->info(), it->vertex(3)->info()});
    
    return ret;
}

template<class TKernel>
void
CGALTriangulation<TKernel>::write(std::string fname)
{
    std::ofstream file(fname);
    
    int nv = mesh.number_of_vertices();
    int nt = 0;
    
    // should not be necessary as there should be only one infinite vertex
    for(auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it)
        if(it->info() != -1) ++nt;
    
    file << nv << " " << nt << "\n";
    
    for(auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it)
        if(it->info() != -1) file << it->point().x() << " " << it->point().y() << " " << it->point().z() << "\n";
    
    for(auto it = mesh.finite_cells_begin(); it != mesh.finite_cells_end(); ++it)
        if(it->info() != -1 && !mesh.is_infinite(it))
        {
            file << it->vertex(0)->info() << " " <<   it->vertex(1)->info() <<  " " <<  it->vertex(2)->info() << " " <<  it->vertex(3)->info() <<"\n";
        }
    
    file.close();
}

template<class TKernel>
void
CGALTriangulation<TKernel>::read(std::string fname)
{
    IndexedTetMesh indexed;
    indexed.read(fname);
    indexed.convert(*this);
}


template<class TKernel>
std::vector<int>
CGALTriangulation<TKernel>::surfaceVertices() const
{
    using namespace std;
    vector<typename Triangulation::Vertex_handle> adj;
    std::vector<int> ret;
    mesh.finite_adjacent_vertices(mesh.infinite_vertex(), back_inserter(adj));
    
    for(auto h : adj) ret.push_back(h->info());
    
    return ret;
}

template<class TKernel>
double
CGALTriangulation<TKernel>::meanEdgeLengthSquared()
{
    // does not account for boundary
    
    double ret = .0;
    int cnt = 0;
    
    for(auto h : mesh.finite_vertex_handles())
    {
        assert(h->info() != -1);
        std::vector<typename Triangulation::Vertex_handle> adj;
        mesh.finite_adjacent_vertices(h, back_inserter(adj));
        
        for(auto hi : adj)
            if(!mesh.is_infinite(hi))
            {
                assert(hi->info() != -1);
                ret += sqrt((h->point() - hi->point()).squared_length());
                ++cnt;
            }
    }
    
    ret /= cnt;
    
    return ret * ret;
}

template<class TKernel>
void
CGALTriangulation<TKernel>::setPoints(const Eigen::MatrixXd& V)
{
    if(V.rows() == mesh.number_of_vertices())
    {
        for(auto h : mesh.finite_vertex_handles())
        {
            if(h->info() == -1 || h->info() >= V.rows()) std::cout << "err" << std::endl;
            
            h->point() = Point(V(h->info(), 0), V(h->info(), 1), V(h->info(), 2));
        }
        
    } else std::cout << "V.rows() does not match" << std::endl;
}

template<class TKernel>
std::vector<double>
CGALTriangulation<TKernel>::tetrahedraVolumes()
{
    std::vector<double> volumes(mesh.number_of_finite_cells());
    
    for(auto it = mesh.finite_cells_begin(); it != mesh.finite_cells_end(); ++it)
    {
        if(it->info() >= volumes.size()) std::cout << "error 2" << std::endl;
        volumes[it->info()] = mesh.tetrahedron(it).volume();
    }
        
    return volumes;
}

template<class TKernel>
void
CGALTriangulation<TKernel>::getPoints(Eigen::MatrixXd& V)
{
    V.resize(mesh.number_of_vertices(), 3);
    
    for(auto h : mesh.finite_vertex_handles())
    {
        for(int k = 0; k < 3; ++k)
        {
            V(h->info(), k) = h->point()[k];
        }
    }
}


template<class TKernel>
int
CGALTriangulation<TKernel>::centerVertex()
{
    double mean[3]{0,0,0};
    for(auto h : mesh.finite_vertex_handles())
    {
        mean[0] += h->point().x();
        mean[1] += h->point().y();
        mean[2] += h->point().z();
    }
            
    Point c(mean[0] / mesh.number_of_vertices(), mean[1] / mesh.number_of_vertices(), mean[2] / mesh.number_of_vertices());
 
    int minId = -1;
    double minDist = std::numeric_limits<double>::max();
    
    for(auto h : mesh.finite_vertex_handles())
    {
        double d = (h->point() - c).squared_length();
        if(d < minDist)
        {
            minDist = d;
            minId = h->info();
        }
    }
    
    return minId;
}

struct Edge
{
    int i, j;
    
    Edge(const int i_, const int j_)
    : i(std::min(i_, j_)), j(std::max(i_, j_))
    {
    }
    
    bool operator==(const Edge& e) const
    {
        return i == e.i && j == e.j;
    }
    
    bool operator<(const Edge& e) const
    {
        return (i < e.i) || (i == e.i && j < e.j);
    }
};
    
namespace std
{
    template<>
    struct hash<Edge> {
    public:
        size_t operator()(const Edge& e) const
        {
                return std::hash<int>()(e.i) ^ (std::hash<int>()(e.j) << 1);
        }
    };
}



template<class TKernel>
void
CGALTriangulation<TKernel>::marchingTets(const Eigen::VectorXd& x, Eigen::MatrixXd& V, Eigen::MatrixXi& F, const double val)
{
    assert(x.rows() == mesh.number_of_vertices());
    
    using namespace std;
    
    const int mtEdges[6][2]
    {
        {0, 1},
        {0, 2},
        {0, 3},
        {1, 2},
        {2, 3},
        {3, 1}
    };

    const int mtVertices[16][4]
    {
        {-1, -1, -1, -1},
        {0, 2, 1, -1},
        {0, 3, 5, -1},
        {1, 3, 5, 2},
        {1, 4, 3, -1},
        {4, 3, 0, 2},
        {1, 4, 5, 0},
        {5, 2, 4, -1},
        {2, 5, 4, -1},
        {0, 5, 4, 1},
        {2, 0, 3, 4},
        {4, 1, 3, -1},
        {2, 5, 3, 1},
        {3, 0, 5, -1},
        {2, 0, 1, -1},
        {-1, -1, -1, -1}
    };

    const int tetFaces[4][3]
    {
        {1, 2, 3},
        {0, 3, 2},
        {0, 1, 3},
        {0, 2, 1}
    };
    
    auto linearInterpolate = [&](const Point& p0, const Point& p1, const double v0, const double v1)
    {
        const double t = v0 / (v0 - v1);
        return p0 + t * (p1 - p0);
    };
    
    std::vector<Point> pts;
    std::vector<std::array<int, 3>> tris;
    
    std::unordered_map<Edge, int> vertexMap;
    
    std::array<int, 4> vids;
    
    for(auto h : mesh.finite_cell_handles())
    {
        const int t[4]{h->vertex(0)->info(), h->vertex(1)->info(), h->vertex(2)->info(), h->vertex(3)->info()};
        const double xh[4]{x(t[0]) - val, x(t[1]) - val, x(t[2]) - val, x(t[3]) - val};
        
        unsigned char code = 0;
          
        for(int k = 0; k < 4; ++k)
            if(xh[k] > 0)
                code |= (1 << k);
        
        if(code && code != 15)
        {
            int k = 0;
            for(; k < 4; ++k)
            {
                int eid = mtVertices[code][k];
                if(eid == -1)
                {
                    k = 0;
                    break;
                }
                
                const int id0 = mtEdges[eid][0];
                const int id1 = mtEdges[eid][1];
                
                int vid0 = t[id0];
                int vid1 = t[id1];
                
                auto it = vertexMap.find(Edge(vid0, vid1));
                
                if(it == vertexMap.end())
                {
                    vids[k] = pts.size();
                    vertexMap[Edge(vid0, vid1)] = pts.size();
                    pts.push_back(linearInterpolate(h->vertex(id0)->point(), h->vertex(id1)->point(), xh[id0], xh[id1]));
                } else
                {
                    vids[k] = it->second;
                }
            }
            
            tris.push_back({vids[1], vids[0], vids[2]});

            if(k)
            {
                tris.push_back({vids[2], vids[0], vids[3]});
            }
        }
    }
    
    V.resize(pts.size(), 3);
    F.resize(tris.size(), 3);
    
    for(int i = 0; i < pts.size(); ++i)
    {
        V(i, 0) = pts[i][0];
        V(i, 1) = pts[i][1];
        V(i, 2) = pts[i][2];
    }
    
     for(int i = 0; i < tris.size(); ++i)
     {
         F(i, 0) = tris[i][0];
         F(i, 1) = tris[i][1];
         F(i, 2) = tris[i][2];
     }
}
  
template<class TKernel>
std::vector<int>
CGALTriangulation<TKernel>::surfaceMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    std::vector<int> idMap(mesh.number_of_vertices(), -1);
    std::vector<int> ids;
    std::vector<std::array<int, 3>> tris;
    std::vector<std::array<double, 3>> pts;
    
    int cnt = 0;
    
    for(auto h : mesh.finite_cell_handles())
    {
        for(int i = 0; i < 4; ++i)
        {
            const int vid = mesh.mirror_vertex(h, i)->info();
            std::array<int, 3> f;
            
            if(vid == -1)
            {
                for(int k = 0; k < 3; ++k)
                {
                    auto vh = h->vertex(Triangulation::vertex_triple_index(i, k));
                    int id = vh->info();
                
                    if(idMap[id] == -1)
                    {
                        pts.push_back(std::array<double, 3>{vh->point()[0], vh->point()[1], vh->point()[2]});
                        ids.push_back(id);
                        idMap[id] = cnt++;
                    }
                    
                    f[k] = idMap[id];
                }
                
                tris.push_back(f);
            }
        }
    }
    
    V.resize(cnt, 3);
    
    for(int i = 0; i < cnt; ++i)
    {
        for(int k = 0; k < 3; ++k)
            V(i, k) = pts[i][k];
    }
  
    F.resize(tris.size(), 3);
    
    for(int i = 0; i < tris.size(); ++i)
    {
        F(i, 0) = tris[i][0];
        F(i, 1) = tris[i][2];
        F(i, 2) = tris[i][1];
    }
    
    return ids;
}
    
template<class TKernel>
std::vector<std::vector<int>>
CGALTriangulation<TKernel>::cutMesh(const std::array<double, 4>& plane, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    std::vector<char> sideFlag(mesh.number_of_vertices(), 3);
    std::vector<char> boundaryFlag(mesh.number_of_vertices(), 0);
    
    for(int i : surfaceVertices()) boundaryFlag[i] = 1;
    
    typedef typename TKernel::Plane_3 Plane;
    
    Plane p(plane[0], plane[1], plane[2], plane[3]);
    
    for(auto h : mesh.finite_vertex_handles())
    {
        sideFlag[h->info()] = 1 + p.has_on_positive_side(h->point());
    }
    
    std::vector<Point> tris;
    std::vector<int> ids;
	std::vector<int> cellids;
    
    for(auto h : mesh.finite_cell_handles())
    {
        for(int i = 0; i < 4; ++i)
        {
            const int vid = mesh.mirror_vertex(h, i)->info();
            int face[3]
            {
                Triangulation::vertex_triple_index(i, 0),
                Triangulation::vertex_triple_index(i, 1),
                Triangulation::vertex_triple_index(i, 2)
            };
            
            if(vid == -1 || sideFlag[vid] == 1)
            {
                if((sideFlag[h->vertex(face[0])->info()] |  sideFlag[h->vertex(face[1])->info()] |  sideFlag[h->vertex(face[2])->info()] ) == 2)
                {
                    for(int k = 0; k < 3; ++k)
                    {
                        tris.push_back(h->vertex(face[k])->point());
                        ids.push_back(h->vertex(face[k])->info());
                    }
					cellids.push_back(h->info());
                }
            }
        }
    }
    
    
    V.resize(tris.size(), 3);
    F.resize(tris.size() / 3, 3);
    
    for(int i = 0; i < tris.size(); ++i)
    {
        V(i, 0) = tris[i][0];
        V(i, 1) = tris[i][1];
        V(i, 2) = tris[i][2];
    }
    
    for(int i = 0; i < tris.size() / 3; ++i)
    {
        F(i, 0) = 3 * i;
        F(i, 1) = 3 * i + 2;
        F(i, 2) = 3 * i + 1;
    }
	std::vector<std::vector<int>> returnvec;
	returnvec.push_back(ids);
	returnvec.push_back(cellids);
    
    return returnvec;
}

template<class TKernel>
void
CGALTriangulation<TKernel>::writeDual(const std::string fname)
{
    std::vector<std::vector<std::vector<Point>>> dual(mesh.number_of_vertices());
    
    std::vector<char> flag(mesh.number_of_vertices(), 1);
    auto brd = surfaceVertices();
    for(int i : brd) flag[i] = 0;
    
    for(auto h : mesh.finite_vertex_handles())
    {
        if(flag[h->info()])
        {
            std::vector<typename Triangulation::Edge> adj;
            mesh.incident_edges(h, back_inserter(adj) );
        
            std::vector<std::vector<Point>> vorCell;
            
            for(auto& e : adj)
            {
                auto start = mesh.incident_cells(e);
                auto it = start;
                
                std::vector<Point> vorFace;
                
                do
                {
                    vorFace.push_back(CGAL::circumcenter(mesh.tetrahedron(it)));
                    
                } while(++it != start);
            
                vorCell.push_back(vorFace);
            }
        
            
            dual[h->info()] = vorCell;
        }
    }
    
    std::ofstream file("../dual");
    for(auto& c : dual) file << c.size() << " ";
    file << "\n";
    
    
    for(auto& c : dual)
        for(auto& f : c)
        {
            for(auto& p : f)
            {
                file << p.x() << " " << p.y() << " " << p.z() << " ";
            }
            
            file << "\n";
        }
    
    file.close();
}

template<class TKernel>
void
CGALTriangulation<TKernel>::DECLaplacian(Eigen::SparseMatrix<double>& L, Eigen::SparseMatrix<double>* M)
{
    std::vector<Eigen::Triplet<double>> triplets;
    const int nv = mesh.number_of_vertices();
    
    // turn off some costly sanity tests
    bool dbg = false;
    
    std::vector<typename TKernel::Vector_3> vecs(mesh.number_of_vertices());
    
    if(M)
    {
        M->resize(nv, nv);
        M->resizeNonZeros(nv);
        
        for(int i = 0; i < nv; ++i)
        {
            M->outerIndexPtr()[i] = i;
            M->innerIndexPtr()[i] = i;
            M->valuePtr()[i] = .0;
        }
        
        M->outerIndexPtr()[nv] = nv;
    }
    
    for(auto h : mesh.finite_cell_handles())
    {
        auto tet = mesh.tetrahedron(h);
        double vol = tet.volume();
        
        for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
            if( i != j )
            {
                const int k = Triangulation::next_around_edge(i, j);
                const int l = Triangulation::next_around_edge(j, i);
            
                auto a2 = tet[i] - tet[l];
                auto b = tet[j] - tet[l];
                auto c = tet[k] - tet[l];
                
                auto b2 = tet[i] - tet[k];
                auto c2 = tet[j] - tet[i];
                auto a =  tet[k] - tet[j];
                
                 
                const double n = (a * b2) * (b2 * c2) * (c * a2) + (b2 * c2) * (c2 * a) * (a2 * b)
                + (c2 * a) * (a * b2) * (b * c) + (a * b2) * (b2 * c2) * (c2 * a);
              
                const double d = 192 * vol * 0.25 * CGAL::cross_product(b2, a).squared_length();
                const double fac = std::abs(d) < 1e-24 ? .0 : 1. / d;
                const double val = -n * (a * b2) * fac;
                
     
                const int r = h->vertex(i)->info();
                const int s = h->vertex(j)->info();
            
                if(M)
                {
                    M->valuePtr()[r] += val * (c2 * c2) / 6.;
                    M->valuePtr()[s] += val * (c2 * c2) / 6.;
                }
    
                triplets.emplace_back(r, r, -val);
                triplets.emplace_back(s, s, -val);
                    
                triplets.emplace_back(r, s, val);
                triplets.emplace_back(s, r, val);
                
                

                if(dbg)
                {
                    auto cc = CGAL::circumcenter(tet);
                    auto ccf = CGAL::circumcenter(tet[i], tet[j], tet[k]);
                    auto cce = CGAL::circumcenter(tet[i], tet[j]);
                    
                    auto edge = tet[j] - tet[i];
            
                    auto nrml = 0.5 * CGAL::cross_product(ccf - cce, cc - cce);
                    double val2 = (nrml * edge) / edge.squared_length();
 
                    if(std::abs(val - val2) > 1e-10) std::cout << "error: " << val << " " << val2 << std::endl;
                  
                    vecs[r] += nrml;
                    vecs[s] -= nrml;
                }
            }
    }

    
    
   
    L.resize(nv, nv);
    L.setFromTriplets(triplets.begin(), triplets.end());

    
    if(dbg)
    {
        Eigen::MatrixXd V(nv, 3);
        
        for(auto h : mesh.finite_vertex_handles())
        {
            V(h->info(), 0) = h->point().x();
            V(h->info(), 1) = h->point().y();
            V(h->info(), 2) = h->point().z();
        }
        
        Eigen::MatrixXd LV = L * V;
        
        for(int i : surfaceVertices())
        {
            LV.row(i).setZero();
        }
        
        std::cout << "linear precision " << LV.norm() << std::endl;
    }
}


template<class TKernel>
void
CGALTriangulation<TKernel>::FEMLaplacian(Eigen::SparseMatrix<double>& L)
{
    Eigen::SparseMatrix<double> G, M;
    gradientOperator(G);
    massMatrix3(M);
    
    L = G.transpose() * M * G;
}

template<class TKernel>
void
CGALTriangulation<TKernel>::gradientOperator(Eigen::SparseMatrix<double>& G)
{
    std::vector<Eigen::Triplet<double>> triplets;
    const int nv = mesh.number_of_vertices();
    
    int cnt = 0;
    
    for(auto h : mesh.finite_cell_handles())
    {
        const int cid = h->info();

        if(cid != -1)
        {
            const double vol = mesh.tetrahedron(h).volume();
            
            for(int j = 0; j < 4; ++j)
            {
                auto t = mesh.triangle(h, j);
                auto va = -CGAL::cross_product(t[1] - t[0], t[2] - t[0]) / 2 / (3 * vol);
                
                const int vid = h->vertex(j)->info();
                assert(vid != -1);
                
                for(int k = 0; k < 3; ++k)
                {
                    triplets.emplace_back(3 * cid + k, vid, va[k]);
                }
            }
            
            ++cnt;
        }
    }
    
    G.resize(3 * cnt, nv);
    G.setFromTriplets(triplets.begin(), triplets.end());
}

template<class TKernel>
void
CGALTriangulation<TKernel>::massMatrix(Eigen::SparseMatrix<double>& M)
{
    std::vector<Eigen::Triplet<double>> triplets;
    const int nv = mesh.number_of_vertices () ;
    
    int cnt = 0;
    
    for(auto h : mesh.finite_cell_handles())
    {
        if(h->info() != -1)
        {
            const double vol = 0.25 * mesh.tetrahedron(h).volume();
            
            for(int i = 0; i < 4; ++i)
            {
                triplets.emplace_back(h->vertex(i)->info(), h->vertex(i)->info(), vol);
            }
     
            ++cnt;
        }
    }
    
    M.resize(nv, nv);
    M.setFromTriplets(triplets.begin(), triplets.end());
}

template<class TKernel>
void
CGALTriangulation<TKernel>::massMatrix3(Eigen::SparseMatrix<double>& M)
{
    std::vector<Eigen::Triplet<double>> triplets;
    const int nv = mesh.number_of_vertices () ;
    
    int cnt = 0;
    
    for(auto h : mesh.finite_cell_handles())
    {
        const int cid = h->info();
        if(cid != -1)
        {
            const double vol = mesh.tetrahedron(h).volume();
            
            for(int k = 0; k < 3; ++k)
                triplets.emplace_back(3 * cid + k, 3 * cid + k, vol);
            
            ++cnt;
        }
    }
    
    M.resize(3 * cnt, 3 * cnt);
    M.setFromTriplets(triplets.begin(), triplets.end());
}


template<class TKernel>
void
CGALTriangulation<TKernel>::umbrellaLaplacian(Eigen::SparseMatrix<double>& L)
{
    std::vector<Eigen::Triplet<double>> triplets;
    
    int cnt = 0;
    for(auto h : mesh.finite_vertex_handles())
    {
        const int cid = h->info();
        if(cid != -1)
        {
            std::vector<typename Triangulation::Vertex_handle> adj;
            mesh.finite_adjacent_vertices(h, back_inserter(adj));
           
            int cnti = 0;
            for(auto hi : adj)
            {
                if(hi->info() != -1)
                {
                    triplets.emplace_back(cid, hi->info(), 1.);
                    ++cnti;
                }
            }
            
            triplets.emplace_back(cid, cid, -cnti);
            
            ++cnt;
        }
    }
    
    L.resize(cnt, cnt);
    L.setFromTriplets(triplets.begin(), triplets.end());
}


template<class TKernel>
void
CGALTriangulation<TKernel>::calcMinAngleAllCells(Eigen::VectorXd &A) {

	A.resize(mesh.number_of_finite_cells());
	double minangle, angle;

    for(auto h : mesh.finite_cell_handles())
    {
		minangle = std::numeric_limits<double>::max();
        for(int i = 0; i < 3; ++i)
        {
            int face[3]
            {
                Triangulation::vertex_triple_index(i, 0),
                Triangulation::vertex_triple_index(i, 1),
                Triangulation::vertex_triple_index(i, 2)
            };
			//std::cout << "--" << std::endl;
			for(int j=0; j<3; ++j){
				angle = CGAL::approximate_dihedral_angle(h->vertex(i)->point(), h->vertex(face[j])->point(), h->vertex(face[(j+1)%3])->point(), h->vertex(face[(j+2)%3])->point());
				if (angle < minangle) minangle=angle;
				// approx angle of edge i, face[j]	
				//std::cout << "(" << i << ", " << face[j] << ")" << std::endl;
			}
		}


        const int cid = h->info();
		A(cid) = minangle;

	}
	//std::cout << std::endl;
}

template<class TKernel>
void
CGALTriangulation<TKernel>::calcVolumeAllCells(Eigen::VectorXd &V) {

    int N = mesh.number_of_finite_cells();
	V.resize(N);
	std::cout << "N: " << N << std::endl;

    for(auto h : mesh.finite_cell_handles())
    {
        const int cid = h->info();
        const double vol = mesh.tetrahedron(h).volume();

		V(cid) = vol;

	}
	std::cout << std::endl;
}


template<class TKernel>
void
CGALTriangulation<TKernel>::calcAMIPSAllCells(Eigen::VectorXd &E) {

    int N = mesh.number_of_finite_cells();
	E.resize(N);

	Eigen::Matrix3d tet, unittet, J;

	// unit tetrahedron in columns
	unittet <<	1,     0.5,					  0.5,
				0, sqrt(3)/2,			sqrt(3)/6,
				0,       0,		sqrt(2) / sqrt(3);

    for(auto h : mesh.finite_cell_handles())
    {
        const int cid = h->info();
		
		for(int col=0; col<3; ++col) {
			tet(0,col) = h->vertex(col+1)->point().x() - h->vertex(0)->point().x();
			tet(1,col) = h->vertex(col+1)->point().y() - h->vertex(0)->point().y();
			tet(2,col) = h->vertex(col+1)->point().z() - h->vertex(0)->point().z();
		}

		J = tet * unittet.inverse();
		E(cid) = (J.transpose()*J).trace() / J.determinant();
	}
}



// TODO: implement
// WARNING: this method resets the cell indices.
template<class TKernel>
void
CGALTriangulation<TKernel>::performRandomFlips(int num_flips, int try_its,  double edge_prob) {
	
	unsigned long j;
	srand( (unsigned)time(NULL) );
	std::random_device rd{}; // use to seed the rng
    std::mt19937 rng{rd()}; // rng
	std::bernoulli_distribution distribution(edge_prob);

	int flipped = 0;

	std::cout << "try to perform " << num_flips << " flips (edge prob " << edge_prob << ")" << std::endl;

	for (int t=0; t<try_its; ++t) {
		// don't try for more than try_its iterations
		if (flipped >= num_flips) {
			break;	
		}
		//std::cout << t+1 << "/" << try_its << std::endl;
		if (distribution(rng)) {
		//if(rand() % 2) {
            // try to flip an edge
			int strt = rand() % mesh.number_of_finite_edges();
			bool flipped_edge = false;
			int e = 0;
			for (auto a: mesh.finite_edges()) {
				if(e >= strt) {
					if(mesh.flip(a)) {
						//std::cout << "Performed Edgeflip (0)" << std::endl;
						++flipped;
						flipped_edge = true;
						break;
					}
				}
				++e;
			}

			if(flipped_edge) {
				continue;	
			} else {
				e = 0;
				for (auto a: mesh.finite_edges()) {
					if(e < strt) {
						if(mesh.flip(a)) {
							//std::cout << "Performed Edgeflip (1)" << std::endl;
							++flipped;
							break;
						}
					}
					++e;
				}
			
			}
		
		} else {
			// try to flip a facet	
			int strt = rand() % mesh.number_of_finite_facets();
			bool flipped_facet = false;
			int e = 0;
			for (auto a: mesh.finite_facets()) {
				if(e >= strt) {
					if(mesh.flip(a)) {
						//std::cout << "Performed Facetflip (0)" << std::endl;
						++flipped;
						flipped_facet = true;
						break;
					}
				}
				++e;
			}

			if(flipped_facet) {
				continue;	
			} else {
				e = 0;
				for (auto a: mesh.finite_facets()) {
					if(e < strt) {
						if(mesh.flip(a)) {
							//std::cout << "Performed Facetflip (1)" << std::endl;
							++flipped;
							break;
						}
					}
					++e;
				}
			
			}
		}
	}
	if (flipped < num_flips) {
		std::cout << "... only managed to perform " << flipped << "flips" << std::endl;
	}
	// reset vertex ids (old ones are no longer valid)
	setIndizes();
}
