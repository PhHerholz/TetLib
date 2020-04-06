#include "tetgen.h"
#include <string>
#include <chrono>
#include <iostream>


#include "IndexedTetMesh.hpp"

template<class TPolyhedron>
void
tetgenMeshPolyhedron(const TPolyhedron& p, CGALTriangulation<typename TPolyhedron::Kernel>& tri, const double area)
{
    using namespace std;
    using namespace std::chrono;
    
    tetgenio in, out;
    
    in.firstnumber = 0;
    in.numberofpoints = p.poly.size_of_vertices();
    in.pointlist = new REAL[in.numberofpoints * 3];
    
    
    auto it = p.poly.vertices_begin();
    
    for(int i = 0; i < in.numberofpoints; ++i, ++it)
        for(int j = 0; j < 3; ++j)
            in.pointlist[3 * i + j] = it->point()[j];
    
    in.numberoffacets = p.poly.size_of_facets();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    fill_n(in.facetmarkerlist, in.numberoffacets, 1);
    
    
    auto itf = p.poly.facets_begin();
    
    for(int i = 0; i < in.numberoffacets; ++i, ++itf)
    {
        auto f = &in.facetlist[i];
        f->numberofpolygons = 1;
        f->numberofholes = 0;
        f->polygonlist = new tetgenio::polygon;
        f->holelist = NULL;
        auto p = f->polygonlist;
        p->numberofvertices = 3;
        p->vertexlist = new int[3];
        
        p->vertexlist[0] = itf->facet_begin()->vertex()->id;
        p->vertexlist[1] = itf->facet_begin()->next()->vertex()->id;
        p->vertexlist[2] = itf->facet_begin()->opposite()->vertex()->id;
    }
    
    tetgenbehavior settings;
    string opts = string("pq1.414a") + to_string(area);
    
    settings.parse_commandline((char*)opts.c_str());
    settings.quiet = 1;
    
    tetrahedralize(&settings, &in, &out);
    
    
    IndexedTetMesh mesh;
    
    mesh.vertices.resize(out.numberofpoints);
    copy_n(out.pointlist, 3 * out.numberofpoints, (double*)mesh.vertices.data() );
    
    mesh.tets.resize(out.numberoftetrahedra);
    copy_n(out.tetrahedronlist, 4 * out.numberoftetrahedra, (int*)mesh.tets.data() );
    
    mesh.convert(tri);
}


