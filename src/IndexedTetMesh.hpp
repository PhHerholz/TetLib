#pragma once

#include <vector>
#include <array>
#include <Eigen/Sparse>
#include <utility>
#include <unordered_map>

template<typename T>
struct std::hash<std::array<T, 2>> {
public:
    size_t operator()(const std::array<T, 2>& a) const
    {
        return std::hash<T>()(a[2]) ^ (std::hash<T>()(a[1]) << 1);
    }
};

template<class Kernel>
class CGALTriangulation;



// Represents a tetrahedral mesh using a list of points and a list of tetrahedra vertex indizes.
// Neighbourhood information is kept explicitly.
// This class is usefull to load and save meshes as well as convert between different formats (e.g. Tetgen, CGAL Triangulation, CGAL embedded Triangulation)

class IndexedTetMesh
{
private:
    std::vector<int> flag;
    
    std::vector<int>
    sortRing(const std::vector<std::array<int, 2>>& ring);
    
    template<class TTriangulationDS>
    typename TTriangulationDS::Vertex_handle
    buildCgalMesh(TTriangulationDS& t);
    
public:
    // could use CGAL or Eigen datatype. The idea was to keep the class self contained. Not sure if this is a good idea.
    typedef typename std::array<double, 3> Point;

    std::vector<Point> vertices;
    std::vector<std::array<unsigned int, 4>> tets;
    std::vector<std::array<int, 4>> tetNeighbours;
    std::vector<std::vector<int>> vertexNeighbours;
    std::unordered_map<std::array<unsigned int, 2>, std::vector<int>> edgeCellNeighbours;
  
    IndexedTetMesh();
    
    ~IndexedTetMesh();
    
    void
    buildTetNeighbours();
    
    void
    buildVertexNeighbours();
   
    void
    buildEdgeNeighours();

    std::vector<char>
    boundaryVertexFlag();
    
    void
    write(std::string fname);
   
    void
    read(std::string fname);
    
    template<class Kernel>
    void
    convert(CGALTriangulation<Kernel>& tri);
};

#include "IndexedTetMesh_impl.h"

