#pragma once


template<class TTriangulationDS>
typename TTriangulationDS::Vertex_handle
IndexedTetMesh::buildCgalMesh(TTriangulationDS& tds)
{
    using namespace std;
    typedef typename TTriangulationDS::Vertex_handle Vertex_handle;
    
    constexpr int faceIds[4][3]
    {
        {1,2,3},
        {0,3,2},
        {0,1,3},
        {0,2,1}
    };
    
    typedef TTriangulationDS TDS;
    typedef typename TDS::Vertex::Point Point3;
       
    tds.clear();
    tds.set_dimension(3);
    
    vector<Vertex_handle> handles;
    
    int cnt = 0;
    for(auto p : vertices)
    {
        auto vh = tds.create_vertex();
        vh->point() = Point3(p[0], p[1], p[2]);
        vh->info() = cnt++;
        handles.push_back(vh);
    }

    typename TTriangulationDS::Vertex_handle infiniteVertex = tds.create_vertex();
    
    handles.push_back(infiniteVertex);
    infiniteVertex->info() = -1;
    
    /* set cell neighbors */
    auto backup = tetNeighbours;
    
    if(tetNeighbours.size() != tets.size())
        buildTetNeighbours();
    
    /* insert dummy (infinite) tets and vertex into data structure */
    const int oldSize = tets.size();
    const unsigned int infVert = vertices.size();
    
    int j = 0;
    for(auto& nbh: tetNeighbours)
    {
        for(int i = 0; i < 4; ++i)
        {
            if(nbh[i] == -1)
                tets.push_back({tets[j][faceIds[i][0]], tets[j][faceIds[i][1]], tets[j][faceIds[i][2]], infVert});
        }
        
        ++j;
    }
    
    vertices.push_back({0,0,0});
    
    /* now set extended cell neighbors */
    buildTetNeighbours();
    
    vector<typename TDS::Cell_handle> cellHandles;
    cnt = 0;
    
    for(auto t : tets)
    {
        auto ch = tds.create_cell(handles[t[0]],
                                  handles[t[1]],
                                  handles[t[2]],
                                  handles[t[3]]);
        
        for(int j = 0; j < 4; ++j)
            handles[t[j]]->set_cell(ch);
        
        ch->info() = cnt >= oldSize ? -1 : cnt++;
        cellHandles.push_back(ch);
    }
    
    
    const int NT = tetNeighbours.size();
    
    for(int i = 0; i < NT; ++i)
    {
        for(int j = 0; j < 4; ++j)
        {
            int id = tetNeighbours[i][j];
            if(id != -1) cellHandles[i]->set_neighbor(j, cellHandles[id]);
            else assert(0);
        }
    }
      
    /* clean up: roll back data structures */
    tets.erase(tets.begin() + oldSize, tets.end());
    vertices.pop_back();
    
    tetNeighbours.swap(backup);
    
    return infiniteVertex;
}

template<class Kernel>
void
IndexedTetMesh::convert(CGALTriangulation<Kernel>& tri)
{
    typename CGALTriangulation<Kernel>::TriangulationDS triDS;
    tri.mesh.set_infinite_vertex(buildCgalMesh(tri.mesh.tds()));
}

