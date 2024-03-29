#pragma once

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <vector>
#include "IndexedTetMesh.hpp"

template<class TKernel>
class CGALTriangulation
{
public:
    typedef TKernel Kernel;
    typedef typename CGAL::Triangulation_vertex_base_with_info_3<int, Kernel> VB;
    typedef typename CGAL::Triangulation_cell_base_with_info_3<int, Kernel> CB;
    typedef CGAL::Triangulation_data_structure_3<VB, CB>  TriangulationDS;
    typedef CGAL::Triangulation_3<Kernel, TriangulationDS> Triangulation;
    typedef CGAL::Delaunay_triangulation_3<Kernel, TriangulationDS> Delaunay;
    typedef typename TKernel::Point_3 Point;
    typedef typename TKernel::Vector_3 Vector;
    typedef typename Triangulation::Vertex_handle Vertex_handle;
    
    // actual data is represented here
    Triangulation mesh;
   
    // get indizes of boundary vertices
    std::vector<int>
    surfaceVertices() const;
   
    // extract all faces that form the boundary of all cells that are completly on the positive side of the plane.
    // return indizes of original vertices. Output mesh consists of individual triangles.
    std::vector<int>
    cutMesh(const std::array<double, 4>& plane, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
    
    std::vector<int>
    surfaceMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F);
    
    // extract a level set of the piecewise linear function induced by x
    void
    marchingTets(const Eigen::VectorXd& x, Eigen::MatrixXd& V, Eigen::MatrixXi& F, const double val = 0.0);
    
    // set vertex and cell indizes. Infinite vertex and cells get index -1.
    void
    setIndizes();
        
    // construct the DEC Laplacian, optionally together with the circumcentric mass matrix
    void
    DECLaplacian(Eigen::SparseMatrix<double>& L, Eigen::SparseMatrix<double>* M = nullptr);
    
    // construct FEM Laplacian
    void
    FEMLaplacian(Eigen::SparseMatrix<double>& L);
    
    // construct the combinatoric Laplacian
    void
    umbrellaLaplacian(Eigen::SparseMatrix<double>& L);
    
    // build the gradient operator (3 |cells| x |vertices|)
    void
    gradientOperator(Eigen::SparseMatrix<double>& G);
    
    // lumped barycentric mass matrix, each entry is repeated three times (3 |cells| x 3 |cells|)
    void
    massMatrix3(Eigen::SparseMatrix<double>& M);
       
    // regular lumped barycentric mass matrix (|vertices| x |vertices|)
    void
    massMatrix(Eigen::SparseMatrix<double>& M);

    // convert to indexed tet set representation
    IndexedTetMesh
    toIndexed() const;
      
    void
    getPoints(Eigen::MatrixXd& V);
    
    void
    setPoints(const Eigen::MatrixXd& V);
    
    std::vector<double>
    tetrahedraVolumes();
    
    // find a vertex that is close to the mean of all others
    int
    centerVertex();
    
    // compute the squared mean of edge lengths (needed for geodesics in heat)
    double
    meanEdgeLengthSquared();
    
    CGALTriangulation();
    
    // intialize with the Delaunay triangulation of some points. The indizes will be given according to the order of 'pts'
    CGALTriangulation(const std::vector<typename TKernel::Point_3>& pts);
    
    CGALTriangulation(Triangulation&& mesh);
    
    ~CGALTriangulation();
    
    // write the polygons of the circumcentric dual (Voronoi diagram in case the triangulation is Delauany)
    void
    writeDual(const std::string fname);
    
    // IO functions
    // data format:
    // "number of vertices" "number of tetraheder"
    // x y z
    // x y z
    // ....
    // i j k l
    // i j k l
    // ...
    
    void
    write(std::string fname);
    
    void
    read(std::string fname);
    
};

#include "CGALTriangulation_impl.hpp"
