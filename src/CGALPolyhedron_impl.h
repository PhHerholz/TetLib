#pragma once

#include <vector>
#include <fstream>

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>


template<class TKernel>
class CGALPolyhedron<TKernel>::Impl
{
    typedef typename TKernel::Point_3 Point3;
    typedef typename Polyhedron::HalfedgeDS             HalfedgeDS;
    typedef typename Polyhedron::Halfedge_handle        Halfedge_handle;
    
    
    // for distance and intersection queries we define an AABB tree.
    typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
    typedef CGAL::AABB_traits<TKernel, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;
    
    std::unique_ptr<Tree> tree;
    
    Polyhedron& poly;
    
    template <class HDS>
    class BuildMesh : public CGAL::Modifier_base<HDS>
    {
    public:
        
        std::string fname;
        std::ifstream instream;
        int NV, NF;

        bool valid = 1;
        
        BuildMesh(const std::string _fname) : fname(_fname)
        {
            // TODO: Better move this code to operator() so that the stream is not open for longer than necessary
            instream = std::ifstream(fname, std::ios::in | std::ios::binary);
            
            if(!instream.is_open())
            {
                valid = 0;
                return;
            }
            
            std::string str;
            instream >> str;
            
            if(str.compare("OFF") != 0)
            {
                valid = 0;
                return;
            }
            
            instream >> NV;
            instream >> NF;
            instream.ignore(INT_MAX, '\n');
        }
    
        ~BuildMesh()
        {}

        void operator()( HDS& hds)
        {
            CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
   
            B.begin_surface(NV, NF);
            typedef typename HDS::Vertex   Vertex;
            
            int cnt = 0;
            
            for(int i = 0; i < NV; ++i)
            {
                double x, y, z;
                instream >> x;
                instream >> y;
                instream >> z;
                B.add_vertex(Point3(x,y,z))->id = cnt++;
            }
            
            for(int l = 0; l < NF; ++l)
            {
                instream.ignore(128,'\n');
                instream.ignore(128,' ');
                
                int i,j,k;
                instream >> i;
                instream >> j;
                instream >> k;
                
                B.begin_facet();
                B.add_vertex_to_facet(i);
                B.add_vertex_to_facet(j);
                B.add_vertex_to_facet(k);
                B.end_facet();
            }
            
            B.end_surface();
            instream.close();
        }
    };

    using Matrix = Eigen::SparseMatrix<double>;
    using EigenPoints = Eigen::Matrix<double, -1, 3>;
    using Triplet = Eigen::Triplet<double>;
    using Triplets = std::vector<Triplet>;
    using Chol = Eigen::SimplicialLDLT<Matrix>;
    
    std::array<int, 3>
    triIds(typename Polyhedron::Facet_handle h)
    {
        auto he = h->halfedge();
        
        return {
            he->vertex()->id,
            he->next()->vertex()->id,
            he->opposite()->vertex()->id
        };
    }
    
    double faceArea(typename Polyhedron::Face_handle fh)
    {
        auto he = fh->halfedge();
        
        return sqrt(typename TKernel::Triangle_3(he->vertex()->point(),
                              he->next()->vertex()->point(),
                              he->next()->next()->vertex()->point()).squared_area());
    }
    
    
    double faceArea(typename Polyhedron::Face_handle fh, const EigenPoints& X)
    {
        const auto t = triIds(fh);
        return .5 * (X.row(t[1]) - X.row(t[0])).cross(X.row(t[2]) - X.row(t[0])).norm();
    }
    
    
    Triplets
    massTriplets(const Eigen::MatrixXd& X)
    {
        Triplets massTriplets;

        int fid = 0;
        for(auto it = poly.facets_begin(); it != poly.facets_end(); ++it)
        {
            double a = faceArea(it, X) / 12.;
            auto f = triIds(it);
            
            for(int i = 0; i < 3; ++i)
            {
                const int v1 = f[(i+1)%3];
                const int v2 = f[(i+2)%3];
                
                massTriplets.emplace_back(v1, v2, a);
                massTriplets.emplace_back(v2, v1, a);
              
                massTriplets.emplace_back(v1, v1, a);
                massTriplets.emplace_back(v2, v2, a);
            }
            
            ++fid;
        }
        
        return massTriplets;
    }
    
    Triplets
    cotanTriplets(const EigenPoints& X)
    {
        std::vector<Triplet> cotanTriplets ;
        
        for(auto it = poly.facets_begin(); it != poly.facets_end(); ++it)
        {
            auto t = triIds(it);
            
            for(int i = 0; i < 3; ++i)
            {
                const int v0 = t[i];
                const int v1 = t[(i+1)%3];
                const int v2 = t[(i+2)%3];
                
                const Eigen::Vector3d e0 = X.row(v1) - X.row(v0);
                const Eigen::Vector3d e1 = X.row(v2) - X.row(v0);
                
                double val = e0.dot(e1);
                
                const double fac = sqrt(e0.dot(e0) * e1.dot(e1) - val * val);
                val = abs(fac) > 1.e-10 ? val / fac : 1.;
                
                cotanTriplets.emplace_back(v1, v2, val);
                cotanTriplets.emplace_back(v2, v1, val);
                cotanTriplets.emplace_back(v1, v1, -val);
                cotanTriplets.emplace_back(v2, v2, -val);
            }
        }
        
        return cotanTriplets;
    }
    
    void rescale(EigenPoints& X)
    {
        Eigen::Vector3d mean = X.colwise().mean();
        
        double maxLen = .0;
        
        for(int i = 0; i < X.rows(); ++i)
        {
            X.row(i) -= mean.transpose();
            if(X.row(i).norm() > maxLen) maxLen = X.row(i).norm();
        }
        
        for(int i = 0; i < X.rows(); ++i)
            X.row(i) /= maxLen;
    }
    
    
    void buildTree()
    {
        if(poly.size_of_facets() == 0) return;
        
        tree.reset(new Tree(faces(poly).first, faces(poly).second, poly));
        tree->accelerate_distance_queries();
    }
    

public:
    
    int load(const std::string fname)
    {
        BuildMesh<HalfedgeDS> builder(fname);
        if(!builder.valid) return 0;
      
        poly.delegate(builder);
        
        buildTree();
        
        return 1;
    }
    
    int write(const std::string fname)
    {
        if(poly.size_of_vertices() == 0) return 0;
        
        std::ofstream file(fname);
        file << "OFF " << poly.size_of_vertices() << " " << poly.size_of_facets() << " 0\n";
        
        std::vector<Point3> pts(poly.size_of_vertices());
        
        for(auto it = poly.vertices_begin(); it != poly.vertices_end(); ++it)
            pts[it->id] = it->point();
        
        for(auto p : pts)
            file << p.x() << " " << p.y() << " " << p.z() << "\n";
        
        for(auto it = poly.facets_begin(); it != poly.facets_end(); ++it)
        {
            auto t = triIds(it);
            file << "3 " << t[0] << " " << t[1] << " " << t[2] << "\n";
        }
        
        file.close();
        
        return 1;
    }
    
    
    class PointWithId : public Point3
    {
    public:
        int id;
        
        PointWithId(const Point3 p, const int _id = -1)
        : Point3(p), id(_id)
        {}
        
        bool operator==(const PointWithId p2)
        {
            return this->x() == p2.x() && this->z() == p2.z() && this->z() == p2.z();
        }
        
        // lexicographic compare
        bool operator<(const PointWithId p2)
        {
            if(this->x() < p2.x()) return true;
            else if(this->x() < p2.x()) return false;
            
            if(this->y() < p2.y()) return true;
            else if(this->y() < p2.y()) return false;
            
            if(this->z() < p2.z()) return true;
            
            return false;
        }
    };
    
    std::vector<PointWithId>
    getPoints(Polyhedron& p)
    {
        std::vector<PointWithId> points;
        int cnt = 0;
       
        for(auto it = p.vertices_begin(); it != p.vertices_end(); ++it)
            points.emplace_back(it->point(), cnt++);
   
        return points;
    }
    
    template<class TPoint>
    void
    setPoints(Polyhedron& p, const std::vector<TPoint>& pts)
    {
        assert(p.size_of_vertices() == pts.size());
        auto it2 = pts.begin();
        
        for(auto it = p.vertices_begin(); it != p.vertices_end(); ++it)
        {
            it->point() = *it2;
            ++it2;
        }
        
        buildTree();
    }
    

    Impl(Impl&& inst)
    : poly(inst.poly)
    {
        buildTree();
    }
    
    Impl(Impl& inst)
    : poly(inst.poly)
    {
        buildTree();
    }
    
    
    Impl(Polyhedron& p)
    : poly(p)
    {
        buildTree();
    }
    
    ~Impl()
    {
        
    }
};


template<class TKernel>
CGALPolyhedron<TKernel>::CGALPolyhedron()
{
    poly.clear();
    impl.reset(new Impl(poly));
}

template<class TKernel>
CGALPolyhedron<TKernel>::CGALPolyhedron(Polyhedron& p)
: impl(new Impl(p)), poly(p)
{
}

template<class TKernel>
CGALPolyhedron<TKernel>::CGALPolyhedron(const CGALPolyhedron& p)
: impl(new Impl(*p.impl)), poly(p.poly)
{
    
}
   
template<class TKernel>
CGALPolyhedron<TKernel>::CGALPolyhedron(CGALPolyhedron&& p)
: impl(new Impl(std::move(*p.impl))), poly(p.poly)
{
    
}

template<class TKernel>
CGALPolyhedron<TKernel>::~CGALPolyhedron()
{}

template<class TKernel>
int CGALPolyhedron<TKernel>::load(const std::string fname)
{
   return impl->load(fname);
}

template<class TKernel>
int CGALPolyhedron<TKernel>::write(const std::string fname)
{
    return impl->write(fname);
}

