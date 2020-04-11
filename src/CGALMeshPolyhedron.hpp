#pragma once

#include "CGALTriangulation.hpp"

template<class TPolyhedron>
void
meshPolyhedron(const TPolyhedron& p, IndexedTetMesh& indexed, const double cellSize = 0.05);

template<class TPolyhedron, class TKernel>
void
meshPolyhedron(const TPolyhedron& p, CGALTriangulation<TKernel>& tri, const double cellSize = 0.05);

template<class TKernel>
typename TKernel::FT
sphere_function (const typename TKernel::Point& p);

template<class TKernel>
void
meshSphere(IndexedTetMesh& indexed, const double cellSize = 0.05);

template<class TKernel2, class TKernel>
void
meshSphere(CGALTriangulation<TKernel>& tri, const double cellSize = 0.05);


#include "CGALMeshPolyhedron_impl.h"
