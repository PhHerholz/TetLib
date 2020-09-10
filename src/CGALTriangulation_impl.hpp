#include <CGAL/Kernel/global_functions.h>
#include <CGAL/intersections.h>
#include <Eigen/Sparse>
#include <fstream>
#include <unordered_map>

#include <igl/grad.h>


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
            ret.tets.push_back(std::array<unsigned int, 4>{(unsigned int) it->vertex(0)->info(), (unsigned int) it->vertex(1)->info(), (unsigned int) it->vertex(2)->info(), (unsigned int) it->vertex(3)->info()});
    
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
std::vector<int>
CGALTriangulation<TKernel>::surfaceVerticesSlow() const
{
    std::vector<int> ret;
	for(auto h : mesh.finite_vertex_handles()) {
		std::vector<CGALTriangulation<Kernel>::Vertex_handle> adj;
		mesh.adjacent_vertices(h, std::back_inserter(adj));
		bool adjtoinf=false;
		for (auto h: adj) {
			if (h->info() == -1) {
				adjtoinf=true;
			}	
		}
		if (adjtoinf) {
			ret.push_back(h->info());
		}
	}
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
    bool dbg = true;
    
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
        
        for(int i : surfaceVerticesSlow())
        {
            LV.row(i).setZero();
        }
        
		std::cout << "dec laplacian " << std::endl;
        std::cout << "linear precision " << LV.norm() << std::endl;
    }
}



template<class TKernel>
void
CGALTriangulation<TKernel>::setLFromW(Eigen::SparseMatrix<double>& L, Eigen::VectorXd& w, std::vector<edge> edges)
{
	std::cout << " generate triplets" << std::endl;
	std::cout << "edges.size() = " << edges.size() << std::endl;
	std::cout << "w.size()     = " << w.size() << std::endl;
    std::vector<Eigen::Triplet<double>> triplets;
	for (int i=0; i < w.size(); ++i){
		triplets.emplace_back(edges[i][0], edges[i][1], w[i]);	
		triplets.emplace_back(edges[i][1], edges[i][0], w[i]);	
		triplets.emplace_back(edges[i][0], edges[i][0], -w[i]);	
		triplets.emplace_back(edges[i][1], edges[i][1], -w[i]);	
	}

	std::cout << " set L from w " << std::endl;
	L.setZero();
    L.setFromTriplets(triplets.begin(), triplets.end());
}
	

template<class TKernel>
void
CGALTriangulation<TKernel>::initAWMatrices(Eigen::SparseMatrix<double>& L, Eigen::VectorXd& w, Eigen::SparseMatrix<double>& A, std::unordered_map<edge, double>& edgeindexmap, std::vector<edge>& edges, std::vector<int> &constrVertices, std::vector<int> ignoreIndices, bool fixBoundaryEdges)
{
	// 1. construct the w vector 
	edgeindexmap.clear();	// to find index for a given edge
	edges.clear();			// to find an edge for a given index

	int ne = mesh.number_of_finite_edges();
	int nv = mesh.number_of_vertices();

	//std::cout << "ne " << ne << " nv " << nv << std::endl;
    w.resize(ne); // weight matrix w
	w.setZero();

	std::cout << L.outerSize() << std::endl;

	//		a) init with non-zero entries of L and build hashmap at the same time	
	int w_index=0;
	for (int k=0; k<L.outerSize(); ++k)
	{
		//std::cout << "k: " << k << std::endl;
		for (Eigen::SparseMatrix<double>::InnerIterator it(L,k); it; ++it)
		{
			if (it.row() != it.col()) {
				edge e{it.row(), it.col()};
				if (edgeindexmap.find(e) == edgeindexmap.end()) {
					w(w_index) = L.coeff(it.row(), it.col());
					edgeindexmap[e] = w_index;
					edges.push_back(e);
					w_index++;
				}
			}
		}
	}

	//		b) check that hashmap has exactly ne

	// 2. Construct A Matrix of constraints
    std::vector<Eigen::Triplet<double>> triplets;
	//		a) iterate over vertices
	
	int constrcnt=0;
	for (auto vh: mesh.finite_vertex_handles()){
		bool ignorevertex=false;
		for (int igni: ignoreIndices) {
			if (vh->info() == igni) {
				// dont add constraints for specific vertices (the surface)
				ignorevertex=true;	
			}	
		}
		if (!ignorevertex){
			constrVertices.push_back(vh->info());
			std::vector<Vertex_handle> adjacent_vertices;
			mesh.finite_adjacent_vertices(vh, std::back_inserter(adjacent_vertices));
			//	b) iterate over adjacent edges
			//		fill in constraint values
			for (Vertex_handle nh: adjacent_vertices) {
				edge e{vh->info(), nh->info()};
				int edgeindex = edgeindexmap[e];

				triplets.emplace_back(3*(constrcnt),   edgeindex, nh->point().x() - vh->point().x());
				triplets.emplace_back(3*(constrcnt)+1, edgeindex, nh->point().y() - vh->point().y());
				triplets.emplace_back(3*(constrcnt)+2, edgeindex, nh->point().z() - vh->point().z());
			}
			++constrcnt;
		} 
	}

	int beccnt = 0;
	if (fixBoundaryEdges) {
		for (edge e: edges) {
			bool boundaryEdge = false;
			for (int igni: ignoreIndices) {
				if (e[0] == igni || e[1] == igni) {
					boundaryEdge = true;	
					break;
				}	
			}
			if (boundaryEdge) {
				int edgeindex = edgeindexmap[e];
				triplets.emplace_back(3*constrcnt+beccnt, edgeindex, 1);
				++beccnt;
			}
		}	
	}

	std::cout << constrcnt << " constraints, " << nv << " vertices." << std::endl;
	A.resize(3*constrcnt + beccnt, ne); // constraint matrix A
	A.setZero();
    A.setFromTriplets(triplets.begin(), triplets.end());
}


template<class TKernel>
double
CGALTriangulation<TKernel>::calcLaplaceTarget(Eigen::VectorXd w, int targetstyle)
{
	if (targetstyle == 0) {
		// expstyle	
		double targetvalue = 0;
		for (int i=0; i < w.size(); ++i) {
			if (w(i) < 0) targetvalue += std::exp(-w(i));
		}
		return targetvalue;
	
	} else if (targetstyle == 1) {
		double targetvalue = 0;
		for (int i=0; i < w.size(); ++i) {
			if (w(i) < 0) targetvalue -= w(i);
		}
		return targetvalue;
	} else {
		//Default: - minval
		return - w.minCoeff();
	}
}

template<class TKernel>
Eigen::VectorXd
CGALTriangulation<TKernel>::calcLaplaceGradient(Eigen::VectorXd w, int targetstyle)
{

	if (targetstyle == 0) {
		// expstyle
		Eigen::VectorXd g(w.size());
		g.setZero();
		for (int i=0; i < w.size(); ++i) {
			if (w(i) < 0) {
				g(i) = - std::exp(-w(i));	
			} 
		}
		return g;	
	} else if (targetstyle == 1) {
		// expstyle
		Eigen::VectorXd g(w.size());
		g.setZero();
		for (int i=0; i < w.size(); ++i) {
			if (w(i) < 0) {
				g(i) = -1;	
			} 
		}
		return g;	
	} else {
		// DEFAULT: minvalue
		// find min indices
		std::vector<int> minIndices;
		double minval  = std::numeric_limits<double>::max();
		// TODO: epsilon correct, needed?
		double eps = 1e-8;
		for (int i=0; i<w.size();++i) {
			if (fabs(w(i) - minval) < eps) {
				// equal case
				minIndices.push_back(i);	
			} else if (w(i) < minval){
				minIndices.clear();
				minIndices.push_back(i);
				minval = w(i);
			}
		}

		// min value gradient
		Eigen::VectorXd g(w.size());
		g.setZero();
		for (int min_ind: minIndices) {
			g(min_ind) = 1.;	
		}
		return -g;	
	}
}

template<class TKernel>
void
CGALTriangulation<TKernel>::DECLaplacianOptimized(Eigen::SparseMatrix<double>& L, double alpha_init, int maxits, int targetstyle, std::vector<int> ignoreIndices, bool fixBoundaryEdges,  std::string logpath)
{
	bool debug = false;

	// 1. INIT A AND w from L
	Eigen::VectorXd w;
	Eigen::SparseMatrix<double> A;
	std::unordered_map<edge, double> edgeindexmap;
	std::vector<edge> edges;

	if (ignoreIndices.size() == 0) {
		ignoreIndices = surfaceVerticesSlow();
	}
	std::vector<int> constrIndices;
	initAWMatrices(L, w, A, edgeindexmap, edges, constrIndices, ignoreIndices, fixBoundaryEdges);
	std::cout << "initialized A and w" << std::endl;

	std::cout << "A.rows : " << A.rows() << ", constrIndices.size(): " << constrIndices.size() << std::endl;

	if (debug) {
		std::cout << "check linear precision of init values:" << std::endl;
		Eigen::VectorXd lpvals = A * w; // weight matrix w
		std::cout << lpvals.norm() << std::endl;
		std::vector<int> problemVertices;

		if (lpvals.norm() >1e-10) {
			std::cout << "No linear precision. highest value: " << std::endl;
			std::cout << lpvals.maxCoeff()<< std::endl;
			std::cout << "dem: ";
			int cntr = 0;
			for (int i=0; i < lpvals.size(); ++i) {
				if (lpvals(i)> 1e-5) {
					problemVertices.push_back(constrIndices[i/3]);
					++cntr;
				}
			}
			std::cout << std::endl;
		}
	}

	// 2. CALC THE GRADIENT
	Eigen::VectorXd gradient;
	// 3. PROJECT THE GRADIENT
	std::cout << "generate solver and start compute step" << std::endl;
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	std::cout << "..generate AAt Matrix..." << std::endl;
	Eigen::SparseMatrix<double> AAt = (A * A.transpose());
	std::cout << "...compute step" << std::endl;
    solver.compute(AAt);
	std::cout << "Finished compute, success?? : " << (solver.info()==Eigen::Success) << std::endl;
	if (! solver.info()==Eigen::Success) {
		std::cout << "...abort" << std::endl;	
		return;
	}
	double targetvalue = calcLaplaceTarget(w, targetstyle); 
	double minvalue    = w.minCoeff();
	//std::cout << "finished compute step, start descent" << std::endl;
	std::cout << "Init targetval: " << targetvalue << std::endl;


	// backtracking line search parameters
	double c   = 0.5;
	double tau = 0.5;
	int lrupdate = 100;
	double alpha = alpha_init; //stepsize;

	bool logtraining = false;
	std::ofstream tlogfile;
	if (logpath.length() > 0) {
		logtraining = true;
		tlogfile.open(logpath);
		// general logging:
		tlogfile << "targetstyle,lrupdate,fixBoudaryEdges,initAWnorm" << std::endl;
		tlogfile << targetstyle << "," <<  lrupdate << "," << fixBoundaryEdges << "," << (A*w).norm() << std::endl;

		tlogfile << "iteration,targetval,minval,stepsize" << std::endl;
		tlogfile << "0," << targetvalue << "," << minvalue << "," << alpha <<  std::endl;
	}

	for (int s_ind=0; s_ind<maxits; ++s_ind) {
		// calc gradient and project it
		gradient = calcLaplaceGradient(w, targetstyle);
		Eigen::VectorXd lambda = solver.solve(2*A*gradient);
		Eigen::VectorXd gradient_projected = gradient - .5*A.transpose()*lambda;

		std::cout << "(A * g_p).norm(): " << (A * gradient_projected).norm() << std::endl;

		if (targetstyle == 0) {
			alpha = alpha_init;
			// Backtracking stepsize
			double m = gradient_projected.norm();
			int j=0;
			while ( calcLaplaceTarget(w - alpha * gradient_projected, targetstyle) > targetvalue - alpha * (c*m) ) {
				alpha = tau * alpha;
			}
		} else {
			if ( (s_ind+1) % lrupdate == 0) {
				alpha = 1./ int((1+s_ind) / lrupdate);
			}
		}

		// update 
		w = w - alpha * gradient_projected;
		double oldtargetvalue = targetvalue;
		targetvalue = calcLaplaceTarget(w, targetstyle);
		minvalue    = w.minCoeff();

		std::cout << "it " << s_ind << ", targetval: " << targetvalue << std::endl;
		std::cout << "         (minval= " << minvalue  << std::endl;
		std::cout << "         (alpha = " << alpha     << ")" << std::endl;

		if (logtraining ) {
			tlogfile << s_ind << "," << targetvalue << "," << minvalue << "," << alpha <<  std::endl;
		}

		if (fabs(targetvalue - oldtargetvalue) < 1e-14) {
			//std::cout << "Alpha < 1e-16, -> break" << std::endl;	
			std::cout << "...seems converged" << std::endl;
			break;
		}
	}	

	if (logtraining) {
		Eigen::VectorXd lpvals = A * w; // weight matrix w
		double endAWnorm = lpvals.norm();

		tlogfile << "endAWnorm" << std::endl;		
		tlogfile << endAWnorm << std::endl;;

		tlogfile.close();	
	}

	// reload the optimized matrix and insert it in the result matrix L
	Eigen::SparseMatrix<double> L_optimized;
	L_optimized.resize(L.rows(), L.cols());
	setLFromW(L_optimized, w, edges);

	if (debug && fixBoundaryEdges) {
		Eigen::SparseMatrix<double> Ldiff = L_optimized - L;
		double ignoresum=0;
		for (int i: ignoreIndices) {
			ignoresum += Ldiff.row(i).norm();	
		}
		std::cout << "Ignoresum: " << ignoresum << std::endl;
	}

	L = L_optimized;
}

template<class TKernel>
void
CGALTriangulation<TKernel>::DECLaplacianRegular(CGALTriangulation<TKernel>::Regular reg, Eigen::SparseMatrix<double>& L, Eigen::SparseMatrix<double>* M)
{
    std::vector<Eigen::Triplet<double>> triplets;
    const int nv = reg.number_of_vertices();
    
    // turn off some costly sanity tests
    bool dbg = true;
	bool dbg2 = false;
    bool dbg3 = false;
    
    std::vector<typename TKernel::Vector_3> vecs(reg.number_of_vertices());
    
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
    
    for(auto h : reg.finite_cell_handles())
    {
        auto tet = reg.tetrahedron(h);
        double vol = tet.volume();
        
        for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
            if( i != j )
            {
                const int k = Triangulation::next_around_edge(i, j);
                const int l = Triangulation::next_around_edge(j, i);

				typedef typename TKernel::Plane_3 Plane;
				typedef typename TKernel::Line_3  Line;
				typedef typename TKernel::Segment_3  Segment;

				auto tet_dual  = reg.dual(h);
				auto face_dual_res = reg.dual(h, l); // facet (tet[i], tet[j], tet[k]);
				Segment face_dual;
				assign(face_dual, face_dual_res);

				// auto ccf = CGAL::circumcenter(tet[i], tet[j], tet[k]);

				Line  fd_line   = face_dual.supporting_line();
				Plane tri_plane(tet[i], tet[j], tet[k]);
				auto ccf_res = intersection(tri_plane, fd_line);
				auto ccf = boost::get<Point>(&*ccf_res);

				// auto cce = CGAL::circumcenter(tet[i], tet[j]);
				auto edge = tet[j] - tet[i];
				Line edge_line(tet[i], edge);
				Plane dual_plane(tet_dual, edge);
				//auto cce = intersection(dual_plane, edge_line);
				
				// access point directly since intersection cannot be a line in this case
				auto cce_res = intersection(dual_plane, edge_line);
				auto cce = boost::get<Point>(&*cce_res);

				auto edge_normalized = edge / sqrt(edge.squared_length());
				Point cce_const = CGAL::ORIGIN + (edge_normalized * (tet_dual - tet[i])) * edge_normalized;

				//std::cout << "Calced" << std::endl;

				// ----------------------------------------------
				// construct ccf hopefully correctly?
				double tri_area = 0.5 * sqrt(CGAL::cross_product(tet[j] - tet[i], tet[k] - tet[i]).squared_length());

				Point ccf_const = tet[i]; 
				// contrib of point k
				auto e_k = (tet[j] - tet[i]) / sqrt((tet[j]-tet[i]).squared_length());
				auto normal_k   = tet[k] - (tet[i] + e_k * ((tet[k]- tet[i]) * e_k));
				normal_k = normal_k / sqrt(normal_k.squared_length()) * sqrt((tet[j] - tet[i]).squared_length());
				ccf_const += (((tet[k] - tet[i]).squared_length() + h->vertex(i)->point().weight() - h->vertex(k)->point().weight()) * normal_k) / (4. * tri_area) ; 
				//contrib of point j
				auto e_j = (tet[k] - tet[i]) / sqrt((tet[k]-tet[i]).squared_length());
				auto normal_j   = tet[j] - (tet[i] + e_j * ((tet[j]- tet[i]) * e_j));
				normal_j = normal_j / sqrt(normal_j.squared_length()) * sqrt((tet[k] - tet[i]).squared_length());
				ccf_const += (((tet[j] - tet[i]).squared_length() + h->vertex(i)->point().weight() - h->vertex(j)->point().weight()) * normal_j) / (4. * tri_area); 

				//auto nrml_const = 0.5 * CGAL::cross_product(ccf_const - cce_const, tet_dual - cce_const); 
				auto nrml_const = 0.5 * CGAL::cross_product(ccf_const - *cce, tet_dual - *cce);
				double val = (nrml_const * edge) / edge.squared_length();
				// ----------------------------------------------
                const int r = h->vertex(i)->info();
                const int s = h->vertex(j)->info();
            
                if(M)
                {
					// TODO: eddge.squared_length the correct thing here?
                    M->valuePtr()[r] += val * (edge.squared_length()) / 6.;
                    M->valuePtr()[s] += val * (edge.squared_length()) / 6.;
                }

				if (dbg2) {
					// compare with non-weighte circumcenter calls (only makes sense if weights are all 0)
                    auto cc_nr  = CGAL::circumcenter(tet);
                    auto ccf_nr = CGAL::circumcenter(tet[i], tet[j], tet[k]);
                    auto cce_nr = CGAL::circumcenter(tet[i], tet[j]);
            
                    auto nrml_nr = 0.5 * CGAL::cross_product(ccf_nr - cce_nr, cc_nr - cce_nr);
                    double val2 = (nrml_nr * edge) / edge.squared_length();

					Point cce_const = tet[i] + 0.5 * (edge.squared_length() + h->vertex(i)->point().weight() - h->vertex(j)->point().weight()) * (edge / sqrt(edge.squared_length()));

					if (std::abs(val - val2) > 1e-13) {

						std::cout << std::endl << "Val  : " << val << std::endl;
						std::cout << "Val2 : " << val2 << std::endl;
						std::cout << "CC  - CC_nr  = " << tet_dual - cc_nr << std::endl;
						std::cout << "CCF - CCF_nr = " << *ccf - ccf_nr << std::endl;
						std::cout << "CCE - CCE_nr = " << *cce - cce_nr << std::endl;

						std::cout << std::endl;
						std::cout << "CCF      : " << *ccf      << std::endl;
						std::cout << "CCF_nr   : " << ccf_nr    << std::endl;
						std::cout << "ccf_const: " << ccf_const << std::endl;

						std::cout << std::endl;
						std::cout << "CCE       : " << *cce      << std::endl;
						std::cout << "CCE_nr    : " << cce_nr    << std::endl;
						std::cout << "CCE_const : " << cce_const << std::endl;

					}

                }

				if (dbg3) {
					// compare to phillipps dec impl (only makes sense if weights are all 0)
					
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
					const double val3 = -n * (a * b2) * fac;

					
					if (std::abs(val - val3) > 1e-13) {
						std::cout << "ERROR: " << val << "!= " << val3 << std::endl;
					}
				}

				if (r > nv) {
					std::cout << "r " << r << " OUT!!!!!!" << std::endl;	
				}
				if (s > nv) {
					std::cout << "s " << s << " OUT!!!!!!" << std::endl;	
				}

                triplets.emplace_back(r, r, -val);
                triplets.emplace_back(s, s, -val);
                    
                triplets.emplace_back(r, s, val);
                triplets.emplace_back(s, r, val);

				/*
				if(std::abs(val - val2) > 1e-10) std::cout << "error: " << val << " " << val2 << std::endl;
			  
				vecs[r] += nrml;
				vecs[s] -= nrml;
				*/

            }
    }
    
   
    L.resize(nv, nv);
    L.setFromTriplets(triplets.begin(), triplets.end());
    
    if(dbg)
    {
		// check linear precision
        Eigen::MatrixXd V(nv, 3);
        
        for(auto h : reg.finite_vertex_handles())
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
CGALTriangulation<TKernel>::DECLaplacianRegular(C3t3 c3t3, Eigen::SparseMatrix<double>& L, Eigen::SparseMatrix<double>* M)
{
    std::vector<Eigen::Triplet<double>> triplets;
    const int nv = c3t3.triangulation().number_of_vertices();
    
    // turn off some costly sanity tests
    bool dbg = true;
	bool dbg2 = false;
    bool dbg3 = false;
    
    std::vector<typename TKernel::Vector_3> vecs(nv);
    
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
    
    //for(auto h : reg.finite_cell_handles())
    for(auto h =  c3t3.cells_in_complex_begin(); h !=  c3t3.cells_in_complex_end(); ++h)
    {
        auto tet = c3t3.triangulation().tetrahedron(h);
        double vol = tet.volume();
        
        for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
            if( i != j )
            {
                const int k = Triangulation::next_around_edge(i, j);
                const int l = Triangulation::next_around_edge(j, i);

				typedef typename TKernel::Plane_3 Plane;
				typedef typename TKernel::Line_3  Line;
				typedef typename TKernel::Segment_3  Segment;

				auto tet_dual  = c3t3.triangulation().dual(h);
				auto face_dual_res = c3t3.triangulation().dual(h, l); // facet (tet[i], tet[j], tet[k]);
				Segment face_dual;
				assign(face_dual, face_dual_res);

				// auto ccf = CGAL::circumcenter(tet[i], tet[j], tet[k]);

				Line  fd_line   = face_dual.supporting_line();
				Plane tri_plane(tet[i], tet[j], tet[k]);
				auto ccf_res = intersection(tri_plane, fd_line);
				auto ccf = boost::get<Point>(&*ccf_res);

				// auto cce = CGAL::circumcenter(tet[i], tet[j]);
				auto edge = tet[j] - tet[i];
				Line edge_line(tet[i], edge);
				Plane dual_plane(tet_dual, edge);
				//auto cce = intersection(dual_plane, edge_line);
				
				// access point directly since intersection cannot be a line in this case
				auto cce_res = intersection(dual_plane, edge_line);
				auto cce = boost::get<Point>(&*cce_res);

				auto edge_normalized = edge / sqrt(edge.squared_length());
				Point cce_const = CGAL::ORIGIN + (edge_normalized * (tet_dual - tet[i])) * edge_normalized;

				//std::cout << "Calced" << std::endl;

				// ----------------------------------------------
				// construct ccf hopefully correctly?
				double tri_area = 0.5 * sqrt(CGAL::cross_product(tet[j] - tet[i], tet[k] - tet[i]).squared_length());

				Point ccf_const = tet[i]; 
				// contrib of point k
				auto e_k = (tet[j] - tet[i]) / sqrt((tet[j]-tet[i]).squared_length());
				auto normal_k   = tet[k] - (tet[i] + e_k * ((tet[k]- tet[i]) * e_k));
				normal_k = normal_k / sqrt(normal_k.squared_length()) * sqrt((tet[j] - tet[i]).squared_length());
				ccf_const += (((tet[k] - tet[i]).squared_length() + h->vertex(i)->point().weight() - h->vertex(k)->point().weight()) * normal_k) / (4. * tri_area) ; 
				//contrib of point j
				auto e_j = (tet[k] - tet[i]) / sqrt((tet[k]-tet[i]).squared_length());
				auto normal_j   = tet[j] - (tet[i] + e_j * ((tet[j]- tet[i]) * e_j));
				normal_j = normal_j / sqrt(normal_j.squared_length()) * sqrt((tet[k] - tet[i]).squared_length());
				ccf_const += (((tet[j] - tet[i]).squared_length() + h->vertex(i)->point().weight() - h->vertex(j)->point().weight()) * normal_j) / (4. * tri_area); 

				//auto nrml_const = 0.5 * CGAL::cross_product(ccf_const - cce_const, tet_dual - cce_const); 
				auto nrml_const = 0.5 * CGAL::cross_product(ccf_const - *cce, tet_dual - *cce);
				double val = (nrml_const * edge) / edge.squared_length();
				// ----------------------------------------------
                const int r = h->vertex(i)->info();
                const int s = h->vertex(j)->info();
            
                if(M)
                {
					// TODO: eddge.squared_length the correct thing here?
                    M->valuePtr()[r] += val * (edge.squared_length()) / 6.;
                    M->valuePtr()[s] += val * (edge.squared_length()) / 6.;
                }

				if (dbg2) {
					// compare with non-weighte circumcenter calls (only makes sense if weights are all 0)
                    auto cc_nr  = CGAL::circumcenter(tet);
                    auto ccf_nr = CGAL::circumcenter(tet[i], tet[j], tet[k]);
                    auto cce_nr = CGAL::circumcenter(tet[i], tet[j]);
            
                    auto nrml_nr = 0.5 * CGAL::cross_product(ccf_nr - cce_nr, cc_nr - cce_nr);
                    double val2 = (nrml_nr * edge) / edge.squared_length();

					Point cce_const = tet[i] + 0.5 * (edge.squared_length() + h->vertex(i)->point().weight() - h->vertex(j)->point().weight()) * (edge / sqrt(edge.squared_length()));

					if (std::abs(val - val2) > 1e-13) {

						std::cout << std::endl << "Val  : " << val << std::endl;
						std::cout << "Val2 : " << val2 << std::endl;
						std::cout << "CC  - CC_nr  = " << tet_dual - cc_nr << std::endl;
						std::cout << "CCF - CCF_nr = " << *ccf - ccf_nr << std::endl;
						std::cout << "CCE - CCE_nr = " << *cce - cce_nr << std::endl;

						std::cout << std::endl;
						std::cout << "CCF      : " << *ccf      << std::endl;
						std::cout << "CCF_nr   : " << ccf_nr    << std::endl;
						std::cout << "ccf_const: " << ccf_const << std::endl;

						std::cout << std::endl;
						std::cout << "CCE       : " << *cce      << std::endl;
						std::cout << "CCE_nr    : " << cce_nr    << std::endl;
						std::cout << "CCE_const : " << cce_const << std::endl;

					}

                }

				if (dbg3) {
					// compare to phillipps dec impl (only makes sense if weights are all 0)
					
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
					const double val3 = -n * (a * b2) * fac;

					
					if (std::abs(val - val3) > 1e-13) {
						std::cout << "ERROR: " << val << "!= " << val3 << std::endl;
					}
				}

				if (r > nv) {
					std::cout << "r " << r << " OUT!!!!!!" << std::endl;	
				}
				if (s > nv) {
					std::cout << "s " << s << " OUT!!!!!!" << std::endl;	
				}

                triplets.emplace_back(r, r, -val);
                triplets.emplace_back(s, s, -val);
                    
                triplets.emplace_back(r, s, val);
                triplets.emplace_back(s, r, val);

				/*
				if(std::abs(val - val2) > 1e-10) std::cout << "error: " << val << " " << val2 << std::endl;
			  
				vecs[r] += nrml;
				vecs[s] -= nrml;
				*/

            }
    }
    
   
    L.resize(nv, nv);
    L.setFromTriplets(triplets.begin(), triplets.end());
    
    if(dbg)
    {
		// check linear precision
        Eigen::MatrixXd V(nv, 3);
        
        for(auto h : c3t3.triangulation().finite_vertex_handles())
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
std::vector<char>
CGALTriangulation<TKernel>::surfaceVertexFlag()
{
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  std::vector<Vertex_handle> out;
  mesh.finite_adjacent_vertices(mesh.infinite_vertex(), std::back_inserter(out));
  std::vector<char> flag(mesh.number_of_vertices(), 0);
  for(auto h : out) flag[h->info()] = 1;
  return flag;
}


template<class TKernel>
void
CGALTriangulation<TKernel>::DECLaplacianMixed(Eigen::SparseMatrix<double>& L, Eigen::SparseMatrix<double>* M)
{
    std::vector<Eigen::Triplet<double>> triplets;
    const int nv = mesh.number_of_vertices();
    auto surfFlag = surfaceVertexFlag();
    
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
    
    
    std::vector<Point> pts(mesh.number_of_vertices());
    
    for(auto& h : mesh.finite_vertex_handles())
        pts[h->info()] = h->point();
    
    
    auto ccFace = [&](const int i, const int j, const int k)
    {
        char sm = surfFlag[i] + surfFlag[j] + surfFlag[k];
        return sm == 0 ? CGAL::circumcenter(pts[i], pts[j], pts[k]) : CGAL::centroid(pts[i], pts[j], pts[k]);
    };
    
    auto ccTet = [&](int tid[4])
    {
        char sm = surfFlag[tid[0]] + surfFlag[tid[1]] + surfFlag[tid[2]] + surfFlag[tid[3]];
        return sm ? CGAL::centroid(pts[tid[0]], pts[tid[1]], pts[tid[2]], pts[tid[3]]) : CGAL::circumcenter(pts[tid[0]], pts[tid[1]], pts[tid[2]], pts[tid[3]]);;
    };
    
    for(auto h : mesh.finite_cell_handles())
    {
        auto tet = mesh.tetrahedron(h);
       
        int tid[4]{
            h->vertex(0)->info(),
            h->vertex(1)->info(),
            h->vertex(2)->info(),
            h->vertex(3)->info()
        };
        
        for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
                if( i != j )
                {
                    const int k = Triangulation::next_around_edge(i, j);
                    
                    auto cc = ccTet(tid);
                    auto ccf = ccFace(tid[i], tid[j], tid[k]);
                    auto cce = CGAL::circumcenter(tet[i], tet[j]);
                    
                    auto edge = tet[j] - tet[i];
                    
                    auto nrml = 0.5 * CGAL::cross_product(ccf - cce, cc - cce);
                    double val = (nrml * edge) / edge.squared_length();
                    
                    
                    const int r = h->vertex(i)->info();
                    const int s = h->vertex(j)->info();
                    
                    if(M)
                    {
                        M->valuePtr()[r] += val * (edge * edge) / 6.;
                        M->valuePtr()[s] += val * (edge * edge) / 6.;
                    }
                    
                    triplets.emplace_back(r, r, -val);
                    triplets.emplace_back(s, s, -val);
                    
                    triplets.emplace_back(r, s, val);
                    triplets.emplace_back(s, r, val);
                }
    }
        
    L.resize(nv, nv);
    L.setFromTriplets(triplets.begin(), triplets.end());
	
	bool dbg = true;
    if(dbg)
    {
		// check linear precision
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
CGALTriangulation<TKernel>::calcMinDECEdgeContributionAllCells(Eigen::VectorXd &A) {


	Eigen::SparseMatrix<double> L;
	DECLaplacianMixed(L);

	A.resize(mesh.number_of_finite_cells());
	double minval, val;

    for(auto h : mesh.finite_cell_handles())
    {
		minval = std::numeric_limits<double>::max();
		//iterate over edges
		for(int i=0; i<3; i++){
			for (int j=i+1; j<4; j++) {
				// edge i-j
				val = L.coeff(h->vertex(i)->info(), h->vertex(j)->info());
				if (val < minval) minval = val;
			}	
		}
		A(h->info()) = minval;
	}
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
CGALTriangulation<TKernel>::calcContainsCircumcenterFlagAllCells(Eigen::VectorXd &C){
    int N = mesh.number_of_finite_cells();
	C.resize(N);
	C.setZero();

    for(auto h : mesh.finite_cell_handles()){
        auto tet = mesh.tetrahedron(h);
		Point cc = CGAL::circumcenter(tet);
        const int cid = h->info();
	
		bool ccc = true; // contains circumcenter
        for(int i = 0; i < 4; ++i)
		{
			if (!ccc) break;
            for(int j = i+1; j < 4; ++j)
            {
				if (!ccc) break;
				// face(i,j,k)
                const int k = Triangulation::next_around_edge(i, j);
                const int l = Triangulation::next_around_edge(j, i);
				if (!tet.has_on_bounded_side(cc)){
					ccc = false;
				}
			}	
		}
		C(cid) = (ccc)?1:0;
	}
}


template<class TKernel>
void
CGALTriangulation<TKernel>::calcIsDelaunayFlagAllCells(Eigen::VectorXd &D){

    int N = mesh.number_of_finite_cells();
	D.resize(N);
	D.setZero();

    for(auto h : mesh.finite_cell_handles()){
        auto tet = mesh.tetrahedron(h);
		Point cc = CGAL::circumcenter(tet);
		double radsq = CGAL::squared_distance(cc, tet[0]);
        const int cid = h->info();
	
		bool del = true; // is delaunay tet flag
		for (auto vh: mesh.finite_vertex_handles()) {
			if (h->has_vertex(vh)) break; // point is part of the tet
			if (CGAL::squared_distance(cc, vh->point()) < radsq) {
				del = false;	
				break;
			}
		}
		D(cid) = (del)?1:0;
	}
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


template<class TKernel>
void  
CGALTriangulation<TKernel>::calcDistToPointAllVertices(Eigen::VectorXd &D, Point p)
{
	D.resize(mesh.number_of_vertices());
    for(auto vh : mesh.finite_vertex_handles()) {	
		D(vh->info()) = sqrt(CGAL::squared_distance(vh->point(), p));
	}
}

template<class TKernel>
void
CGALTriangulation<TKernel>::calcCentroidAllCells(Eigen::MatrixXd &cellCentroids)
{
	cellCentroids.resize(mesh.number_of_finite_cells(), 3);
	for (auto ch: mesh.finite_cell_handles()) {
		auto tet = mesh.tetrahedron(ch);
		// calc and save centroid of the tet
		Point centroid = CGAL::centroid(tet.vertex(0), tet.vertex(1), tet.vertex(2), tet.vertex(3));
		cellCentroids(ch->info(), 0) = centroid.x();
		cellCentroids(ch->info(), 1) = centroid.y();
		cellCentroids(ch->info(), 2) = centroid.z();
	}
}

template<class TKernel>
void  
CGALTriangulation<TKernel>::calcHeatGradientAllCells(Eigen::MatrixXd h, Eigen::MatrixXd &heatGradField)
{
	heatGradField.resize(mesh.number_of_finite_cells(), 3);

	Eigen::MatrixXd F(1, 4);
	for (int i=0; i<4; ++i) F(0,i) = i;

	for (auto ch: mesh.finite_cell_handles()) {
		// calc and save gradient of the tet
		Eigen::MatrixXd V(4, 3);
		Eigen::MatrixXd U(4, 1);
		Eigen::SparseMatrix<double> G(1,3);
		for (int i = 0; i < 4; ++i) {
			V(i,0) = ch->vertex(i)->point().x();
			V(i,1) = ch->vertex(i)->point().y();
			V(i,2) = ch->vertex(i)->point().z();
		}
		// gradient operator call
		igl::grad(V,F,G);

		// apply gradient operator to values
		for (int i = 0; i < 4; ++i) U(i, 0) = h(ch->vertex(i)->info());
		Eigen::MatrixXd GU = Eigen::Map<const Eigen::MatrixXd>((G*U).eval().data(),F.rows(),3);

		for (int j=0; j<3; ++j) {
			heatGradField(ch->info(), j) = GU(0,j);
		}
	}
}

template<class TKernel>
void  
CGALTriangulation<TKernel>::evaluateHeatGradByDirectionPerCell(Eigen::MatrixXd &heatGradField, Eigen::MatrixXd &hGMetricRes)
{
	for (auto ch: mesh.finite_cell_handles()) {
		auto tet = mesh.tetrahedron(ch);

		// centroid as the point in the tet in which the gradient directions are evaluated
		Point centroid = CGAL::centroid(tet.vertex(0), tet.vertex(1), tet.vertex(2), tet.vertex(3));
		double pointnorm = sqrt(CGAL::squared_distance(CGAL::ORIGIN, centroid));
		
		// analytic gradient direction
		Eigen::MatrixXd gradDir(1,3);
		gradDir << centroid.x() / pointnorm, centroid.y() / pointnorm, centroid.z() / pointnorm;

		// heat gradient direction
		Eigen::MatrixXd heatGrad = heatGradField.row(ch->info());

		hGMetricRes(ch->info()) = gradDir * (heatGrad / heatGrad.norm());
	}
}

template<class TKernel>
void
CGALTriangulation<TKernel>::solveConstrainedLS(Eigen::VectorXd &grad, double &b, Eigen::MatrixXd X, Eigen::MatrixXd Y, Eigen::MatrixXd x0, double y0, double stab)
{
	using namespace Eigen;
	//std::cout << "X.rows(): " << X.rows() << std::endl;
	int N = X.rows();
	assert(Y.rows() == N);
	// shapes: X (N, 3), Y (N, 1), x0 (1, 3), y0 (1, 1) 
	
	grad.resize(3);
	grad.setZero();
	
	//std::cout << "p" << std::endl;
	// concat the one dimension to X and x0
	//MatrixXd ons(N, 1); ons.setOnes(N, 1);
	MatrixXd X_ (N, 4);  
	MatrixXd x0_ (1, 4); 
	// X_ << X, ons;
	// // x0_ << x0, 1.;
	//
	for (int i=0; i<N; ++i) {
		for (int j=0; j<3; ++j) {
			X_(i,j) = X(i,j);
		}
		X_(i, 3) = 1;
	}
	for (int j=0; j<3; ++j) {
		x0_(0, j) = x0(j);
	}
	x0_(0, 3) = 1.;

	/*
	std::cout << "X_ (" << X_.rows() << ", " << X_.cols() << "): " << std::endl;
	std::cout << X_ <<  std::endl;
	std::cout << "Y   : " << Y  <<  std::endl;
	std::cout << "x0_ : " << x0_ << std::endl;
	std::cout << "y0  : " << y0  << std::endl;
	*/

	MatrixXd XXtinv = (X_.transpose() * X_ + stab * MatrixXd::Identity(4, 4) ).inverse();

	double lbda_up = (x0_ * XXtinv * X_.transpose() * Y)(0,0) - y0;
	double lbda_do = (x0_ * XXtinv * x0_.transpose()).value();
	double lbda = lbda_up / lbda_do; 
	/*
	std::cout << "lbda_up: " << lbda_up << std::endl;
	std::cout << "lbda_do: " << lbda_do << std::endl;
	*/

	MatrixXd w = XXtinv * (X_.transpose() * Y  - x0_.transpose() * lbda); 

	/*
	std::cout << "XXtinv (" << XXtinv.rows() << ", " << XXtinv.cols() << ") : " << std::endl;
	std::cout << XXtinv << std::endl;
	std::cout << "lbda : " << lbda << std::endl;
	std::cout << "w : " << w << std::endl;
	*/
   
	for (int i=0; i < 3; ++i) {
		grad(i) = w(i, 0);
	}
	b = w(3,0);

	/*
	std::cout << "grad: " << grad << std::endl;
	std::cout << "b : " << b << std::endl;
	std::cout << "y0  " << y0 << std::endl;
	std::cout << "x0_*w" << x0_ * w << std::endl;
	*/
}

template<class TKernel>
void  
CGALTriangulation<TKernel>::calcHeatGradientAllVertices(Eigen::MatrixXd h, Eigen::MatrixXd &heatGradField)
{
	using namespace Eigen;
	heatGradField.resize(h.rows(), 3);

	for (auto vh : mesh.finite_vertex_handles()) {
		std::vector<typename Triangulation::Vertex_handle> adj;
		mesh.finite_adjacent_vertices(vh, back_inserter(adj));

		// fill the required matrices:
		MatrixXd X(adj.size(), 3), Y(adj.size(), 1), x0(1, 3); 

		// point itself:
		x0(0, 0) = vh->point().x();
		x0(0, 1) = vh->point().y();
		x0(0, 2) = vh->point().z();
		double y0 = h(vh->info());

		// the neighbors:
		int ncntr = 0;
		for (auto nh: adj) {
			X(ncntr, 0) = nh->point().x();
			X(ncntr, 1) = nh->point().y();
			X(ncntr, 2) = nh->point().z();
			Y(ncntr, 0) = h(nh->info());

			++ncntr;
		}

		/*
		std::cout << "X: " << std::endl;
		std::cout << X << std::endl;
		std::cout << "Y: " << std::endl;
		std::cout << Y << std::endl;
		std::cout << "x0: " << std::endl;
		std::cout << x0 << std::endl;
		std::cout << "y0: " << std::endl;
		std::cout << y0 << std::endl;
		std::cout << " 4 " << std::endl;
		*/

		// solve 
		VectorXd w(3);
		double b;

		solveConstrainedLS(w, b, X, Y, x0, y0, 0.);
		
		// check constraint and error
		//std::cout << "Check: w.transpose() * x0 + b = " << (w.transpose() * x0.transpose())(0, 0)   + b << std::endl;
		//std::cout << "                           y0 = " << y0 << std::endl;

		//std::cout << "Check: w.transpose() * X + b = " << (w.transpose() *

		// write gradient (w) to the gradfield:
		for (int i=0; i<3; ++i) {
			heatGradField(vh->info(), i) = w(i);
		}
	}
}


// WARNING: this method resets the cell indices.
template<class TKernel>
void
CGALTriangulation<TKernel>::performRandomFlips(int num_flips, int try_its,  double edge_prob) {
	
	unsigned long j;
	srand( (unsigned)time(NULL) );
//	std::random_device rd{}; // use to seed the rng
//    std::mt19937 rng{rd()}; // rng
//	std::bernoulli_distribution distribution(edge_prob);

	int flipped = 0;

	std::cout << "try to perform " << num_flips << " flips (edge prob " << edge_prob << ")" << std::endl;

	for (int t=0; t<try_its; ++t) {
		// don't try for more than try_its iterations
		if (flipped >= num_flips) {
			break;	
		}
		//std::cout << t+1 << "/" << try_its << std::endl;
	//	if (distribution(rng)) {
		if(rand() % 2) {
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
	// reset cell ids (old ones are no longer valid)

    int cnt = 0;
    for(auto h : mesh.finite_cell_handles()) {
        h->info() = cnt++;
    }

}


template<class TKernel>
typename CGALTriangulation<TKernel>::Regular
CGALTriangulation<TKernel>::generateRandomRegular(double variance){

	std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0, variance};
	
	std::vector< std::pair<WPoint,unsigned> > points;
    for (auto vh : mesh.finite_vertex_handles()) {
		double randomweight = fabs(d(gen)); 
		points.push_back( std::make_pair(WPoint(vh->point(),randomweight),vh->info()) );
    }
	std::cout << std::endl;


	//Regular reg;
	CGALTriangulation<TKernel>::Regular reg;
	reg.insert(points.begin(), points.end());
	std::cout << "Created regular triangulation. Valid?=" << reg.is_valid() << std::endl;

	//handle infos: (finite vertex infos have been set when inserting them)
    reg.infinite_vertex()->info() = -1;
	int cnt = 0;
	for(auto it = reg.cells_begin(); it != reg.cells_end(); ++it)
	{
		if(reg.is_infinite(it)) it->info() = -1;
		else it->info() = cnt++;
	}

	return reg;
}


template<class TKernel>
typename CGALTriangulation<TKernel>::Regular
CGALTriangulation<TKernel>::generateRegularFromWeightsfile(std::string weightsfilepath){

	std::vector<double> pointweights;
	std::ifstream weightsfile;
	weightsfile.open(weightsfilepath); 
	std::string line;
	std::getline(weightsfile, line); // skip the header line
	while(std::getline(weightsfile, line)) {
		//std::cout << line << std::endl;	
		pointweights.push_back(std::stod(line));
	}
	weightsfile.close();
	
	std::vector< std::pair<WPoint,unsigned> > points;
	int wind=0;
    for (auto vh : mesh.finite_vertex_handles()) {
		points.push_back( std::make_pair(WPoint(vh->point(),pointweights[wind]),vh->info()) );
		++wind;
    }
	std::cout << std::endl;
	std::cout << "Created file list. pointweights.size(): " << pointweights.size() << ", points.size() " << points.size() << std::endl;

	//Regular reg;
	CGALTriangulation<TKernel>::Regular reg;
	reg.insert(points.begin(), points.end());
	std::cout << "Created reg tri with weights from file. Valid?=" << reg.is_valid() << std::endl;

	//handle infos: (finite vertex infos have been set when inserting them)
    reg.infinite_vertex()->info() = -1;
	int cnt = 0;
	for(auto it = reg.cells_begin(); it != reg.cells_end(); ++it)
	{
		if(reg.is_infinite(it)) it->info() = -1;
		else it->info() = cnt++;
	}

	return reg;
}

template<class TKernel>
void
CGALTriangulation<TKernel>::replaceMeshByRegular(double variance, std::vector<int> &innerShell, std::vector<int> &middleShell, std::vector<int> &outerShell, double minVolume, bool boundary_only, bool removeInner){

	CGALTriangulation<TKernel>::Regular reg = generateRandomRegular(variance);
	replaceMeshByRegular(reg,innerShell, middleShell, outerShell, minVolume, boundary_only, removeInner);
}

template<class TKernel>
void
CGALTriangulation<TKernel>::replaceMeshByRegular(Regular &reg, std::vector<int> &innerShell, std::vector<int> &middleShell, std::vector<int> &outerShell, double minVolume, bool boundary_only, bool removeInner){

	// Translate to IndexedTetmesh
    IndexedTetMesh ret;
	int nv = reg.number_of_vertices();

	ret.vertices.resize(nv);
	std::unordered_map<int, int> idconversion;
	std::unordered_map<int, int> idconversion_inverse;

	int inscounter = 0;
    for(auto it = reg.vertices_begin(); it != reg.vertices_end(); ++it) {
        if(it->info() != -1){
			ret.vertices[inscounter][0] = it->point().x();
			ret.vertices[inscounter][1] = it->point().y();
			ret.vertices[inscounter][2] = it->point().z();
			idconversion[it->info()] = inscounter;
			//idconversion_inverse[inscounter] = it->info();
			inscounter++;
		}
	}

    for(auto it: reg.finite_cell_handles()){

		bool addCell = true;

		// check each cell if (!boundary_only):
		if (removeInner) {
			bool innerCell=true;
			for (int i = 0; i < 4; ++i) {
				auto vh = it->vertex(i);
				if (std::find(innerShell.begin(), innerShell.end(), vh->info()) == innerShell.end()){
					innerCell=false;
					break;
				}
			}
			if (innerCell) addCell=false;
		}

		bool checkvolcell = false;
		if (minVolume > 0) {
			if (!boundary_only){
				std::cout << "Check all Cells" << std::endl;		
				checkvolcell = true;
			}

			// if boundary only: check if cell is on boundary
			if (!checkvolcell) {
				for (int i = 0; i < 4 ; ++i) {
					if (reg.mirror_vertex(it, i)->info() == -1) {
						//std::cout << "Boundary cell " << std::endl;
						checkvolcell = true;
					}
				}
			}

			if (checkvolcell) {
				auto tet = reg.tetrahedron(it);
				double vol = tet.volume();
				//std::cout << "Vol:    " << vol               << std::endl;
				//std::cout << "Minvol: " << minVolume << std::endl;
				if (vol < minVolume){
					//std::cout << "Raus damit!" << std::endl;	
					addCell = false;
				}
			}
		}

        if(it->info() == -1) addCell = false;

		if(addCell)
		{
            ret.tets.push_back(std::array<unsigned int, 4>{(unsigned int) idconversion[it->vertex(0)->info()],
														   (unsigned int) idconversion[it->vertex(1)->info()],
														   (unsigned int) idconversion[it->vertex(2)->info()],
														   (unsigned int) idconversion[it->vertex(3)->info()] });
		}
	}
    
	// replace triangulation with random reg triangulation
	ret.convert(*this);

	// reset vertex infos to old numbering	
    for(auto it = reg.vertices_begin(); it != reg.vertices_end(); ++it) {
        if(it->info() != -1){
			it->info() = idconversion[it->info()];
		}
	}

	std::vector<int> new_inner;
	for (int i: innerShell) {
		if (idconversion.find(i) != idconversion.end()) {
			new_inner.push_back(idconversion[i]);	
		}
	}
	innerShell = new_inner;

	std::vector<int> new_middle;
	for (int i: middleShell) {
		if (idconversion.find(i) != idconversion.end()) {
			new_middle.push_back(idconversion[i]);	
		}
	}
	middleShell = new_middle;

	std::vector<int> new_outer;
	for (int i: outerShell) {
		if (idconversion.find(i) != idconversion.end()) {
			new_outer.push_back(idconversion[i]);	
		}
	}
	outerShell = new_outer;
	
	/*
	for (int i; i < changeind_vecs.size(); ++i) {
		std::vector<int> new_indices;
		for (int i: changeind_vecs[i]) {
			if (idconversion.find(i) != idconversion.end()) {
				new_indices.push_back(idconversion[i]);	
			}
		}
		changeind_vecs[i] = new_indices;
	}
	*/
	/*
	for (auto indvec: changeind_vecs) {
		std::vector<int> new_indices;
		for (int i: indvec) {
			if (idconversion.find(i) != idconversion.end()) {
				new_indices.push_back(idconversion[i]);	
			}
		}
		indvec = new_indices;
	}
	*/

	/*
	 * DEPRECATED, delete soon
	// convert origin and orbit indices
	if (originind >= 0) {
		int new_originind = -1;
		if (idconversion.find(originind) != idconversion.end()) {
			originind = idconversion[originind];
		} else {
			std::cout << "origin lost during noising" << std::endl;	
		}
	}
	std::vector<int> new_orbitinds;
	for (int i: orbitinds) {
		if (idconversion.find(i) != idconversion.end()) {
			new_orbitinds.push_back(idconversion[i]);	
		}
	}
	orbitinds = new_orbitinds;
	*/
}
