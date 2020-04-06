#pragma once

template<class TPolyhedron>
void
tetgenMeshPolyhedron(const TPolyhedron& p, CGALTriangulation<typename TPolyhedron::Kernel>& tri, const double area);

#include "TetgenMeshPolyhedron_impl.h"

