#pragma once

#include "CGALTriangulation.hpp"

template<class TPolyhedron>
void
meshPolyhedron(const TPolyhedron& p, IndexedTetMesh& indexed, const double cellSize = 0.05);

template<class TPolyhedron, class TKernel>
void
meshPolyhedron(const TPolyhedron& p, CGALTriangulation<TKernel>& tri, const double cellSize = 0.05);

#include "CGALMeshPolyhedron_impl.h"


