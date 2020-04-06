#pragma once

#include <memory>
#include <CGAL/Polyhedron_3.h>

template<class TKernel>
class CGALPolyhedron
{
    template <class Refs>
    class VertexWithId : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, typename TKernel::Point_3>
    {
        typedef typename CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, typename TKernel::Point_3> Base;
    public:
        
        int id;
        
        VertexWithId()
        : Base()
        {}
        
        VertexWithId(typename TKernel::Point_3 p)
        : Base(p)
        {}
        
        VertexWithId(typename TKernel::Point_3 p, const int _id)
        : Base(p), id(_id)
        {}
        
    };
    
    struct ItemsId : public CGAL::Polyhedron_items_3 {
        template <class Refs, class Traits>
        struct Vertex_wrapper {
            typedef VertexWithId<Refs> Vertex;
        };
    };
    

    class Impl;
    std::unique_ptr<Impl> impl;
    
public:
    typedef TKernel Kernel;
   
    typedef CGAL::Polyhedron_3<TKernel, ItemsId> Polyhedron;
    
    Polyhedron poly;
    
    CGALPolyhedron();
    
    ~CGALPolyhedron();
   
    CGALPolyhedron(Polyhedron& p);
   
    CGALPolyhedron(const CGALPolyhedron& p);
   
    CGALPolyhedron(CGALPolyhedron&& p);
    
    int
    load(const std::string fname);
    
    int
    write(const std::string fname);
    
    // TODO: add libigl compatible interface
};


#include "CGALPolyhedron_impl.h"

