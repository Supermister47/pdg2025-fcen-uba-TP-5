//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-10-24 18:56:55 taubin>
//------------------------------------------------------------------------
//
// SceneGraphProcessor.cpp
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Brown University nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL GABRIEL TAUBIN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <cmath>
#include <iostream>
#include "SceneGraphProcessor.hpp"
#include "SceneGraphTraversal.hpp"
#include "Shape.hpp"
#include "IndexedLineSet.hpp"
#include "IndexedFaceSet.hpp"
#include "Appearance.hpp"
#include "Material.hpp"
#include "core/Graph.hpp"
#include "core/IsoSurf.hpp"
#include "core/HexGridPartition.hpp"

const int SceneGraphProcessor::_hexGridEdge[12][2] = {
  {0,4}, {1,5}, {2,6}, {3,7},
  {0,2}, {1,3}, {4,6}, {5,7},
  {0,1}, {2,3}, {4,5}, {6,7}
};

SceneGraphProcessor::SceneGraphProcessor(SceneGraph& wrl) :
  _wrl(wrl),
  _nGrid(0),
  _nPoints(0),
  _next(nullptr),
  _first(nullptr) {
}

SceneGraphProcessor::~SceneGraphProcessor() {
  _deletePartition();
}

int SceneGraphProcessor::_nCells() {
  return _nGrid * _nGrid * _nGrid;
}

int SceneGraphProcessor::_createPartition
(Vec3f& min, Vec3f& max, int depth, vector<float>& coord) {
  int nPointsInPartition = 0;
  _deletePartition();
  float dx = max.x - min.x;
  float dy = max.y - min.y;
  float dz = max.z - min.z;
  if (dx > 0.0f && dy > 0.0f && dz > 0.0f) {
    _nPoints = static_cast<int>((coord.size() / 3));
    if (depth < 0) depth = 0;
    _nGrid = 1 << depth;
    int nCells = _nCells();
    _next = new int[_nPoints];
    _first = new int[nCells];
    int iCell, iPoint, ix, iy, iz;
    float x, y, z, nG = static_cast<float>(_nGrid);
    for (iCell = 0;iCell < nCells;iCell++)
      _first[iCell] = -1;
    for (iPoint = 0;iPoint < _nPoints;iPoint++) {
      if ((x = coord[static_cast<size_t>(3 * iPoint)]) < min.x || x > max.x) continue;
      if ((y = coord[static_cast<size_t>(3 * iPoint + 1)]) < min.y || y > max.y) continue;
      if ((z = coord[static_cast<size_t>(3 * iPoint + 2)]) < min.z || z > max.z) continue;
      ix = static_cast<int>((nG * (x - min.x)) / dx);
      iy = static_cast<int>((nG * (y - min.y)) / dy);
      iz = static_cast<int>((nG * (z - min.z)) / dz);
      iCell = ix + _nGrid * (iy + _nGrid * iz);
      _next[iPoint] = _first[iCell];
      _first[iCell] = iPoint;
      nPointsInPartition++;
    }
  }
  return nPointsInPartition;
}

void SceneGraphProcessor::_deletePartition() {
  if (_first != nullptr) delete[] _first;
  if (_next != nullptr) delete[] _next;
  _nGrid = 0;
  _nPoints = 0;
  _next = nullptr;
  _first = nullptr;
}

int SceneGraphProcessor::numberOfShapeNodes() {
  int nShapes = 0;
  SceneGraphTraversal sgt(_wrl);
  Node* node;
  while ((node = sgt.next()) != (Node*)0) {
    if (node->isShape()) nShapes++;
  }
  return nShapes;
}

void SceneGraphProcessor::normalClear() {
  _applyToIndexedFaceSet(_normalClear);
}

void SceneGraphProcessor::normalInvert() {
  _applyToIndexedFaceSet(_normalInvert);
}

void SceneGraphProcessor::computeNormalPerFace() {
  _applyToIndexedFaceSet(_computeNormalPerFace);
}

void SceneGraphProcessor::computeNormalPerVertex() {
  _applyToIndexedFaceSet(_computeNormalPerVertex);
}

void SceneGraphProcessor::computeNormalPerCorner() {
  _applyToIndexedFaceSet(_computeNormalPerCorner);
}

void SceneGraphProcessor::_applyToIndexedFaceSet(IndexedFaceSet::Operator o) {
  SceneGraphTraversal sgt(_wrl);
  Node* node;
  while ((node = sgt.next()) != (Node*)0) {
    if (node->isShape()) {
      Shape* shape = (Shape*)node;
      node = shape->getGeometry();
      if (node != (Node*)0 && node->isIndexedFaceSet()) {
        IndexedFaceSet& ifs = *((IndexedFaceSet*)node);
        o(ifs);
      }
    }
  }
}

void SceneGraphProcessor::_normalClear(IndexedFaceSet& ifs) {
  vector<float>& normal = ifs.getNormal();
  vector<int>& normalIndex = ifs.getNormalIndex();
  ifs.setNormalPerVertex(true);
  normal.clear();
  normalIndex.clear();
}

void SceneGraphProcessor::_normalInvert(IndexedFaceSet& ifs) {
  vector<float>& normal = ifs.getNormal();
  for (int i = 0;i < (int)normal.size();i++)
    normal[i] = -normal[i];
}

void SceneGraphProcessor::_computeFaceNormal
(vector<float>& coord, vector<int>& coordIndex,
  int i0, int i1, Vec3f& n, bool normalize) {
  int niF, iV, i;
  Vec3f p, pi, ni, v1, v2;
  niF = i1 - i0; // number of face corners
  // n << 0,0,0;
  n[0] = n[1] = n[2] = 0.0f;
  if (niF == 3) { // triangle
    // triangle
    iV = coordIndex[i0];
    // p << coord[3*iV  ],coord[3*iV+1],coord[3*iV+2];
    p[0] = coord[3 * iV]; p[1] = coord[3 * iV + 1]; p[2] = coord[3 * iV + 2];
    iV = coordIndex[i0 + 1];
    // v1 << coord[3*iV  ],coord[3*iV+1],coord[3*iV+2];
    v1[0] = coord[3 * iV]; v1[1] = coord[3 * iV + 1]; v1[2] = coord[3 * iV + 2];
    iV = coordIndex[i0 + 2];
    // v2 << coord[3*iV  ],coord[3*iV+1],coord[3*iV+2];
    v2[0] = coord[3 * iV]; v2[1] = coord[3 * iV + 1]; v2[2] = coord[3 * iV + 2];

    // v1 -= p;
    v1[0] -= p[0]; v1[1] -= p[1]; v1[2] -= p[2];
    // v2 -= p;
    v2[0] -= p[0]; v2[1] -= p[1]; v2[2] -= p[2];
    // n = v1.cross(v2);
    n[0] = v1[1] * v2[2] - v1[2] * v2[1];
    n[1] = v1[2] * v2[0] - v1[0] * v2[2];
    n[2] = v1[0] * v2[1] - v1[1] * v2[0];

    if (normalize) {
      // n.normalize();
      float nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
      if (nn > 0.0f) {
        nn = (float)sqrt(nn);
        n[0] /= nn; n[1] /= nn; n[2] /= nn;
      }
    }
  }
  else if (niF > 3) { // polygon
    // compute face centroid
    // p << 0,0,0;
    p[0] = p[1] = p[2] = 0.0f;
    for (i = i0;i < i1;i++) {
      iV = coordIndex[i];
      // pi << coord[3*iV  ],coord[3*iV+1],coord[3*iV+2];
      pi[0] = coord[3 * iV]; pi[1] = coord[3 * iV + 1]; pi[2] = coord[3 * iV + 2];
      // p += pi;
      p[0] += pi[0]; p[1] += pi[1]; p[2] += pi[2];
    }
    // p /= ((float)niF);
    p[0] /= ((float)niF); p[1] /= ((float)niF); p[2] /= ((float)niF);
    // accumulate face normal
    iV = coordIndex[i1 - 1];
    // v1 << coord[3*iV  ],coord[3*iV+1],coord[3*iV+2];
    v1[0] = coord[3 * iV]; v1[1] = coord[3 * iV + 1]; v1[2] = coord[3 * iV + 2];
    // v1 -= p;
    v1[0] -= p[0]; v1[1] -= p[1]; v1[2] -= p[2];
    for (i = i0;i < i1;i++) {
      iV = coordIndex[i];
      // v2 << coord[3*iV  ],coord[3*iV+1],coord[3*iV+2];
      v2[0] = coord[3 * iV]; v2[1] = coord[3 * iV + 1]; v2[2] = coord[3 * iV + 2];
      // v2 -= p;
      v2[0] -= p[0]; v2[1] -= p[1]; v2[2] -= p[2];

      // ni = v1.cross(v2);
      ni[0] = v1[1] * v2[2] - v1[2] * v2[1];
      ni[1] = v1[2] * v2[0] - v1[0] * v2[2];
      ni[2] = v1[0] * v2[1] - v1[1] * v2[0];

      // n += ni;
      n[0] += ni[0]; n[1] += ni[1]; n[2] += ni[2];

      // v1.swap(v2);
      v1[0] = v2[0]; v1[1] = v2[1]; v1[2] = v2[2];
    }
    if (normalize) {
      // n.normalize();
      float nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
      if (nn > 0.0f) {
        nn = (float)sqrt(nn);
        n[0] /= nn; n[1] /= nn; n[2] /= nn;
      }
    }
  }
  else /* if(n<3) */ { // face with less than 3 vertices
    // throw exception ?
    // n == (0,0,0)
  }
}

void SceneGraphProcessor::_computeNormalPerFace(IndexedFaceSet& ifs) {
  if (ifs.getNormalBinding() == IndexedFaceSet::PB_PER_FACE) return;
  vector<float>& coord = ifs.getCoord();
  vector<int>& coordIndex = ifs.getCoordIndex();
  vector<float>& normal = ifs.getNormal();
  vector<int>& normalIndex = ifs.getNormalIndex();
  ifs.setNormalPerVertex(false);
  normal.clear();
  normalIndex.clear();
  Vec3f n;
  int /*iF,*/ i0, i1;
  for (i0 = i1 = 0;i1 < (int)coordIndex.size();i1++) {
    if (coordIndex[i1] < 0) {
      _computeFaceNormal(coord, coordIndex, i0, i1, n, true);
      normal.push_back((float)(n[0]));
      normal.push_back((float)(n[1]));
      normal.push_back((float)(n[2]));
      i0 = i1 + 1; /* iF++; */
    }
  }
}

void SceneGraphProcessor::_computeNormalPerVertex(IndexedFaceSet& ifs) {
  if (ifs.getNormalBinding() == IndexedFaceSet::PB_PER_VERTEX) return;
  vector<float>& coord = ifs.getCoord();
  vector<int>& coordIndex = ifs.getCoordIndex();
  vector<float>& normal = ifs.getNormal();
  vector<int>& normalIndex = ifs.getNormalIndex();
  ifs.setNormalPerVertex(true);
  normal.clear();
  normalIndex.clear();
  Vec3f n;
  int /*iF,*/ nV, i, i0, i1, iV;
  float x0, x1, x2;
  nV = (int)(coord.size() / 3);
  // initialize accumulators
  normal.insert(normal.end(), coord.size(), 0.0f);
  // accumulate face normals
  for (i0 = i1 = 0;i1 < (int)coordIndex.size();i1++) {
    if (coordIndex[i1] < 0) {
      _computeFaceNormal(coord, coordIndex, i0, i1, n, false);
      // accumulate
      for (i = i0;i < i1;i++) {
        iV = coordIndex[i];
        x0 = normal[3 * iV];
        x1 = normal[3 * iV + 1];
        x2 = normal[3 * iV + 2];
        normal[3 * iV] = x0 + ((float)(n[0]));
        normal[3 * iV + 1] = x1 + ((float)(n[1]));
        normal[3 * iV + 2] = x2 + ((float)(n[2]));
      }
      i0 = i1 + 1; /* iF++; */
    }
  }
  for (iV = 0;iV < nV;iV++) {
    n[0] = normal[3 * iV];
    n[1] = normal[3 * iV + 1];
    n[2] = normal[3 * iV + 2];
    float nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
    if (nn > 0.0f) {
      nn = (float)sqrt(nn);
      n[0] /= nn; n[1] /= nn; n[2] /= nn;
    }
    normal[3 * iV] = n[0];
    normal[3 * iV + 1] = n[1];
    normal[3 * iV + 2] = n[2];
  }
}

void SceneGraphProcessor::_computeNormalPerCorner(IndexedFaceSet& ifs) {
  if (ifs.getNormalBinding() == IndexedFaceSet::PB_PER_CORNER) return;

  vector<float>& coord = ifs.getCoord();
  vector<int>& coordIndex = ifs.getCoordIndex();

  vector<float>& normal = ifs.getNormal();
  vector<int>& normalIndex = ifs.getNormalIndex();
  ifs.setNormalPerVertex(true);
  normal.clear();
  normalIndex.clear();

  int nC = (int)coordIndex.size();

  Vec3f pP, p0, pN, vP, vN, n;
  int /*iF,*/iC, iC0, iC1, iCp, iCn, niF, iN, iVp, iV0, iVn;
  for (/*iF=*/iC0 = iC1 = 0;iC1 < nC;iC1++) {
    if (coordIndex[iC1] >= 0) continue;
    niF = iC1 - iC0; // number of face corners
    if (niF >= 3) { // polygon
      // n << 0,0,0;
      n[0] = n[1] = n[2] = 0.0f;
      for (iC = iC0;iC < iC1;iC++) {

        iCp = (iC0 < iC) ? iC - 1 : iC1 - 1;
        iCn = (iC + 1 < iC1) ? iC + 1 : iC0;

        iVp = coordIndex[iCp];
        pP[0] = coord[3 * iVp]; pP[1] = coord[3 * iVp + 1]; pP[2] = coord[3 * iVp + 2];
        iV0 = coordIndex[iC];
        p0[0] = coord[3 * iV0]; p0[1] = coord[3 * iV0 + 1]; p0[2] = coord[3 * iV0 + 2];
        iVn = coordIndex[iCn];
        pN[0] = coord[3 * iVn]; pN[1] = coord[3 * iVn + 1]; pN[2] = coord[3 * iVn + 2];

        vP[0] = pP[0] - p0[0]; vP[1] = pP[1] - p0[1]; vP[2] = pP[2] - p0[2];
        vN[0] = pN[0] - p0[0]; vN[1] = pN[1] - p0[1]; vN[2] = pN[2] - p0[2];

        // n = vN.cross(vP);
        n[0] = vN[1] * vP[2] - vN[2] * vP[1];
        n[1] = vN[2] * vP[0] - vN[0] * vP[2];
        n[2] = vN[0] * vP[1] - vN[1] * vP[0];

        // n.normalize();
        float nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
        if (nn > 0.0f) {
          nn = (float)sqrt(nn);
          n[0] /= nn; n[1] /= nn; n[2] /= nn;
        }

        iN = (int)(normal.size() / 3);
        normal.push_back(n[0]);
        normal.push_back(n[1]);
        normal.push_back(n[2]);
        normalIndex.push_back(iN);
      }
      normalIndex.push_back(-1);

    }
    else /* if(niF<3) */ { // face with less than 3 vertices
      // throw exception ?
      n[0] = n[1] = n[2] = 0.0f;
      for (iC = iC0;iC < iC1;iC++) {
        iN = (int)(normal.size() / 3);
        normal.push_back(n[0]);
        normal.push_back(n[1]);
        normal.push_back(n[2]);
        normalIndex.push_back(iN);
      }
      normalIndex.push_back(-1);
    }

    iC0 = iC1 + 1; /*iF++;*/
  }
}

void SceneGraphProcessor::gridAdd
(int depth, float scale, bool isCube) {
  const string name = "GRID";
  Shape* shape = (Shape*)0;
  const Node* node = _wrl.getChild(name);
  if (node == (Node*)0) {
    shape = new Shape();
    shape->setName(name);
    Appearance* appearance = new Appearance();
    shape->setAppearance(appearance);
    Material* material = new Material();
    // colors should be stored in WrlViewerData
    Color gridColor(0.5f, 0.3f, 0.0f);
    material->setDiffuseColor(gridColor);
    appearance->setMaterial(material);
    _wrl.addChild(shape);
  }
  else if (node->isShape()) {
    shape = (Shape*)node;
  }
  if (shape == (Shape*)0) { /* throw exception ??? */ return; }

  IndexedLineSet* ils = (IndexedLineSet*)0;
  node = shape->getGeometry();
  if (node == (Node*)0) {
    ils = new IndexedLineSet();
    shape->setGeometry(ils);
  }
  else if (node->isIndexedLineSet()) {
    ils = (IndexedLineSet*)node;
  }
  if (ils == (IndexedLineSet*)0) { /* throw exception ??? */ return; }

  vector<float>& coord = ils->getCoord();
  vector<int>& coordIndex = ils->getCoordIndex();
  vector<float>& color = ils->getColor();
  vector<int>& colorIndex = ils->getColorIndex();
  coord.clear();
  coordIndex.clear();
  color.clear();
  colorIndex.clear();
  ils->setColorPerVertex(true);

  _wrl.updateBBox();
  Vec3f& center = _wrl.getBBoxCenter();
  Vec3f& size = _wrl.getBBoxSize();

  float dx = size.x / 2.0f;
  float dy = size.y / 2.0f;
  float dz = size.z / 2.0f;
  if (isCube) {
    float dMax = dx; if (dy > dMax) dMax = dy; if (dz > dMax) dMax = dz;
    dx = dMax; dy = dMax; dz = dMax;
  }
  if (scale > 0.0f) {
    dx *= scale; dy *= scale; dz *= scale;
  }

  float x0 = center.x - dx; float y0 = center.y - dy; float z0 = center.z - dz;
  float x1 = center.x + dx; float y1 = center.y + dy; float z1 = center.z + dz;

  if (depth == 0) {

    // vertices
    coord.push_back(x0); coord.push_back(y0); coord.push_back(z0);
    coord.push_back(x0); coord.push_back(y0); coord.push_back(z1);
    coord.push_back(x0); coord.push_back(y1); coord.push_back(z0);
    coord.push_back(x0); coord.push_back(y1); coord.push_back(z1);
    coord.push_back(x1); coord.push_back(y0); coord.push_back(z0);
    coord.push_back(x1); coord.push_back(y0); coord.push_back(z1);
    coord.push_back(x1); coord.push_back(y1); coord.push_back(z0);
    coord.push_back(x1); coord.push_back(y1); coord.push_back(z1);

    // edges
    coordIndex.push_back(0); coordIndex.push_back(1); coordIndex.push_back(-1);
    coordIndex.push_back(2); coordIndex.push_back(3); coordIndex.push_back(-1);
    coordIndex.push_back(4); coordIndex.push_back(5); coordIndex.push_back(-1);
    coordIndex.push_back(6); coordIndex.push_back(7); coordIndex.push_back(-1);
    //
    coordIndex.push_back(0); coordIndex.push_back(2); coordIndex.push_back(-1);
    coordIndex.push_back(1); coordIndex.push_back(3); coordIndex.push_back(-1);
    coordIndex.push_back(4); coordIndex.push_back(6); coordIndex.push_back(-1);
    coordIndex.push_back(5); coordIndex.push_back(7); coordIndex.push_back(-1);
    // iz=0
    coordIndex.push_back(0); coordIndex.push_back(4); coordIndex.push_back(-1);
    coordIndex.push_back(1); coordIndex.push_back(5); coordIndex.push_back(-1);
    coordIndex.push_back(2); coordIndex.push_back(6); coordIndex.push_back(-1);
    coordIndex.push_back(3); coordIndex.push_back(7); coordIndex.push_back(-1);

  }
  else {

    int N = 1 << depth;

    // vertices
    float x, y, z;
    int ix, iy, iz, jx, jy, jz, iV0, iV1;

    for (iz = 0, jz = N;iz <= N;iz++, jz--) {
      z = (((float)jz) * z0 + ((float)iz) * z1) / ((float)N);
      for (iy = 0, jy = N;iy <= N;iy++, jy--) {
        y = (((float)jy) * y0 + ((float)iy) * y1) / ((float)N);
        for (ix = 0, jx = N;ix <= N;ix++, jx--) {
          x = (((float)jx) * x0 + ((float)ix) * x1) / ((float)N);
          coord.push_back(x); coord.push_back(y); coord.push_back(z);
        }
      }
    }

    // edges
    for (iz = 0;iz < N;iz++) {
      for (iy = 0;iy <= N;iy++) {
        for (ix = 0;ix <= N;ix++) {
          iV0 = (ix)+(N + 1) * ((iy)+(N + 1) * (iz));
          iV1 = (ix)+(N + 1) * ((iy)+(N + 1) * (iz + 1));
          coordIndex.push_back(iV0);
          coordIndex.push_back(iV1);
          coordIndex.push_back(-1);
        }
      }
    }
    for (iz = 0;iz <= N;iz++) {
      for (iy = 0;iy < N;iy++) {
        for (ix = 0;ix <= N;ix++) {
          iV0 = (ix)+(N + 1) * ((iy)+(N + 1) * (iz));
          iV1 = (ix)+(N + 1) * ((iy + 1) + (N + 1) * (iz));
          coordIndex.push_back(iV0);
          coordIndex.push_back(iV1);
          coordIndex.push_back(-1);
        }
      }
    }
    for (iz = 0;iz <= N;iz++) {
      for (iy = 0;iy <= N;iy++) {
        for (ix = 0;ix < N;ix++) {
          iV0 = (ix)+(N + 1) * ((iy)+(N + 1) * (iz));
          iV1 = (ix + 1) + (N + 1) * ((iy)+(N + 1) * (iz));
          coordIndex.push_back(iV0);
          coordIndex.push_back(iV1);
          coordIndex.push_back(-1);
        }
      }
    }

  }

}

void SceneGraphProcessor::gridAdd(HexGridPartition& hgp) {
  const string name = "GRID";
  Shape* shape = (Shape*)0;
  const Node* node = _wrl.getChild(name);
  if (node == (Node*)0) {
    shape = new Shape();
    shape->setName(name);
    Appearance* appearance = new Appearance();
    shape->setAppearance(appearance);
    Material* material = new Material();
    // colors should be stored in WrlViewerData
    Color gridColor(0.5f, 0.3f, 0.0f);
    material->setDiffuseColor(gridColor);
    appearance->setMaterial(material);
    _wrl.addChild(shape);
  }
  else if (node->isShape()) {
    shape = (Shape*)node;
  }
  if (shape == (Shape*)0) { /* throw exception ??? */ return; }

  IndexedLineSet* ils = (IndexedLineSet*)0;
  node = shape->getGeometry();
  if (node == (Node*)0) {
    ils = new IndexedLineSet();
    shape->setGeometry(ils);
  }
  else if (node->isIndexedLineSet()) {
    ils = (IndexedLineSet*)node;
  }
  if (ils == (IndexedLineSet*)0) { /* throw exception ??? */ return; }

  vector<float>& coord = ils->getCoord();
  vector<int>& coordIndex = ils->getCoordIndex();
  vector<float>& color = ils->getColor();
  vector<int>& colorIndex = ils->getColorIndex();
  coord.clear();
  coordIndex.clear();
  color.clear();
  colorIndex.clear();
  ils->setColorPerVertex(true);

  Vec3f gridMin = hgp.getMin();
  Vec3f gridMax = hgp.getMax();
  int   N = hgp.getResolution();

  float x, x0, x1, y, y0, y1, z, z0, z1;
  x0 = gridMin.x; y0 = gridMin.y; z0 = gridMin.z;
  x1 = gridMax.x; y1 = gridMax.y; z1 = gridMax.z;

  map<int, int> vMap;
  Graph graph;

  // int iPointFirst;
  int h, iCell, iC, iCx, iCy, iCz, iVx, iVy, iVz, iVertex, iV[8], iV0, iV1, iE;
  map<int, int> first = hgp.getFirstMap();
  map<int, int>::iterator iMap;

  for (iMap = first.begin();iMap != first.end();iMap++) {
    iCell = iMap->first;
    iC = iCell; iCx = iC % N; iC /= N; iCy = iC % N; iC /= N; iCz = iC;

    // iPointFirst = iMap->second;

    // vertices
    //
    // 0 ---- 1
    // |\     |\
    // | 2 ---- 3
    // 4 +--- 5 |
    //  \|     \|
    //   6 ---- 7
    //
    // 0 : (iCx  ,iCy  ,iCz  )
    // 1 : (iCx+1,iCy  ,iCz  )
    // 2 : (iCx  ,iCy+1,iCz  )
    // 3 : (iCx+1,iCy+1,iCz  )
    // 4 : (iCx  ,iCy  ,iCz+1)
    // 5 : (iCx+1,iCy  ,iCz+1)
    // 6 : (iCx  ,iCy+1,iCz+1)
    // 7 : (iCx+1,iCy+1,iCz+1)
    for (h = 0;h < 8;h++) {
      iVx = iCx + ((h) % 2);
      iVy = iCy + ((h / 2) % 2);
      iVz = iCz + ((h / 4) % 2);
      iVertex = iVx + (N + 1) * (iVy + (N + 1) * iVz);
      if (vMap.count(iVertex) == 0) {
        iV[h] = static_cast<int>(coord.size() / 3);
        x = (((float)(N - iVx)) * x0 + ((float)iVx) * x1) / ((float)N);
        coord.push_back(x);
        y = (((float)(N - iVy)) * y0 + ((float)iVy) * y1) / ((float)N);
        coord.push_back(y);
        z = (((float)(N - iVz)) * z0 + ((float)iVz) * z1) / ((float)N);
        coord.push_back(z);
        vMap[iVertex] = iV[h];
      }
      else /* if(vMap.count(iVertex])==0) */ {
        iV[h] = vMap[iVertex];
      }
    }

    for (h = 0;h < 12;h++) {
      iV0 = iV[_hexGridEdge[h][0]];
      iV1 = iV[_hexGridEdge[h][1]];
      iE = graph.getEdge(iV0, iV1);
      if (iE < 0) {
        graph.insertEdge(iV0, iV1, true);
        coordIndex.push_back(iV0);
        coordIndex.push_back(iV1);
        coordIndex.push_back(-1);
      }
    }

  } // for(iMap=first.begin();iMap!=first.end();iMap++)
}

void SceneGraphProcessor::gridRemove() {
  vector<pNode>& children = _wrl.getChildren();
  vector<pNode>::iterator i;
  for (i = children.begin();i != children.end();i++)
    if ((*i)->nameEquals("GRID"))
      break;
  if (i != children.end())
    children.erase(i);
}

void SceneGraphProcessor::edgesAdd() {
  SceneGraphTraversal sgt(_wrl);
  const Node* node;
  while ((node = sgt.next()) != (Node*)0) {
    if (node->isShape()) {
      Shape* shape = (Shape*)node;
      const Node* parent = shape->getParent();
      Group* group = (Group*)parent;

      node = shape->getGeometry();
      if (node != (Node*)0 && node->isIndexedFaceSet()) {
        IndexedFaceSet* ifs = (IndexedFaceSet*)node;

        shape->setShow(false);

        // compose the node name ???
        string name = "EDGES";
        node = group->getChild(name);
        if (node == (Node*)0) {
          shape = new Shape();
          shape->setName(name);
          Appearance* appearance = new Appearance();
          shape->setAppearance(appearance);
          Material* material = new Material();
          // colors should be stored in WrlViewerData
          Color edgeColor(1.0f, 0.5f, 0.0f);
          material->setDiffuseColor(edgeColor);
          appearance->setMaterial(material);
          group->addChild(shape);
        }
        else if (node->isShape()) {
          shape = (Shape*)node;
        }
        else /* if(node!=(Node*)0 && node->isShape()==false */ {
          // throw exception ???
        }
        if (shape == (Shape*)0) { /* throw exception ??? */ return; }

        IndexedLineSet* ils = (IndexedLineSet*)0;
        node = shape->getGeometry();
        if (node == (Node*)0) {
          ils = new IndexedLineSet();
          shape->setGeometry(ils);
        }
        else if (node->isIndexedLineSet()) {
          ils = (IndexedLineSet*)node;
        }
        else /* if(node!=(Node*)0 && node->isIndexedLineSet()==false) */ {
          // throw exception ???
        }

        if (ils == (IndexedLineSet*)0) { /* throw exception ??? */ return; }

        ils->clear();

        vector<float>& coordIfs = ifs->getCoord();
        vector<int>& coordIndexIfs = ifs->getCoordIndex();

        vector<float>& coordIls = ils->getCoord();
        vector<int>& coordIndexIls = ils->getCoordIndex();

        coordIls.insert(coordIls.end(),
          coordIfs.begin(), coordIfs.end());

        int i, i0, i1, nV, iV, iV0, iV1/*,iF*/;
        nV = static_cast<int>(coordIfs.size() / 3);
        for (iV = 0;iV < nV;iV++) {
          coordIls.push_back(coordIfs[3 * iV]);
          coordIls.push_back(coordIfs[3 * iV + 1]);
          coordIls.push_back(coordIfs[3 * iV + 2]);
        }

        for (/*iF=*/i0 = i1 = 0;i1 < (int)coordIndexIfs.size();i1++) {
          if (coordIndexIfs[i1] < 0) {
            iV0 = coordIndexIfs[i1 - 1];
            for (i = i0;i < i1;i++) {
              iV1 = coordIndexIfs[i];
              coordIndexIls.push_back(iV0);
              coordIndexIls.push_back(iV1);
              coordIndexIls.push_back(-1);
              iV0 = iV1;
            }
            i0 = i1 + 1; /* iF++; */
          }
        }


      }
    }
  }
}

void SceneGraphProcessor::edgesRemove() {
  SceneGraphTraversal sgt(_wrl);
  Node* node;
  while ((node = sgt.next()) != (Node*)0) {
    if (node->isShape()) {
      Shape* shape = (Shape*)node;
      const Node* parent = shape->getParent();
      Group* group = (Group*)parent;
      vector<pNode>& children = group->getChildren();
      vector<pNode>::iterator i;
      do {
        for (i = children.begin();i != children.end();i++)
          if ((*i)->nameEquals("EDGES"))
            break;
        if (i != children.end()) {
          children.erase(i);
          i = children.begin();
        }
      } while (i != children.end());
    }
  }
}

void SceneGraphProcessor::shapeIndexedFaceSetShow() {
  SceneGraphTraversal sgt(_wrl);
  Node* node;
  while ((node = sgt.next()) != (Node*)0) {
    if (node->isShape()) {
      Shape* shape = (Shape*)node;
      node = shape->getGeometry();
      if (node != (Node*)0 && node->isIndexedFaceSet()) {
        shape->setShow(true);
      }
    }
  }
}

void SceneGraphProcessor::shapeIndexedFaceSetHide() {
  SceneGraphTraversal sgt(_wrl);
  Node* node;
  while ((node = sgt.next()) != (Node*)0) {
    if (node->isShape()) {
      Shape* shape = (Shape*)node;
      node = shape->getGeometry();
      if (node != (Node*)0 && node->isIndexedFaceSet()) {
        shape->setShow(false);
      }
    }
  }
}

void SceneGraphProcessor::shapeIndexedLineSetShow() {
  SceneGraphTraversal sgt(_wrl);
  Node* node;
  while ((node = sgt.next()) != (Node*)0) {
    if (node->isShape()) {
      Shape* shape = (Shape*)node;
      node = shape->getGeometry();
      if (node != (Node*)0 && node->isIndexedLineSet()) {
        shape->setShow(true);
      }
    }
  }
}

void SceneGraphProcessor::shapeIndexedLineSetHide() {
  SceneGraphTraversal sgt(_wrl);
  Node* node;
  while ((node = sgt.next()) != (Node*)0) {
    if (node->isShape()) {
      Shape* shape = (Shape*)node;
      node = shape->getGeometry();
      if (node != (Node*)0 && node->isIndexedLineSet()) {
        shape->setShow(false);
      }
    }
  }
}

bool SceneGraphProcessor::hasGrid() {
  return _wrl.getChild("GRID") != (Node*)0;
}

bool SceneGraphProcessor::_hasShapeProperty(Shape::Property p) {
  bool value = false;
  SceneGraphTraversal sgt(_wrl);
  Node* node;
  while (value == false && (node = sgt.next()) != (Node*)0) {
    if (node->isShape()) {
      Shape& shape = *(Shape*)node;
      value = p(shape);
    }
  }
  return value;
}

bool SceneGraphProcessor::_hasIndexedFaceSetProperty(IndexedFaceSet::Property p) {
  bool value = false;
  SceneGraphTraversal sgt(_wrl);
  Node* node;
  while (value == false && (node = sgt.next()) != (Node*)0) {
    if (node->isShape()) {
      Shape* shape = (Shape*)node;
      if (shape->hasGeometryIndexedFaceSet()) {
        IndexedFaceSet& ifs = *(IndexedFaceSet*)(shape->getGeometry());
        value = p(ifs);
      }
    }
  }
  return value;
}

bool SceneGraphProcessor::_hasIndexedLineSetProperty(IndexedLineSet::Property p) {
  bool value = false;
  SceneGraphTraversal sgt(_wrl);
  Node* node;
  while (value == false && (node = sgt.next()) != (Node*)0) {
    if (node->isShape()) {
      Shape* shape = (Shape*)node;
      if (shape->hasGeometryIndexedLineSet()) {
        IndexedLineSet& ils = *(IndexedLineSet*)(shape->getGeometry());
        value = p(ils);
      }
    }
  }
  return value;
}

bool SceneGraphProcessor::_hasFaces(IndexedFaceSet& ifs) {
  return (ifs.getNumberOfCoord() > 0 && ifs.getNumberOfFaces() > 0);
}

bool SceneGraphProcessor::_hasNormalNone(IndexedFaceSet& ifs) {
  return
    (ifs.getNumberOfCoord() == 0 && ifs.getNumberOfFaces() == 0) ||
    (ifs.getNormalBinding() == IndexedFaceSet::PB_NONE);
}

bool SceneGraphProcessor::_hasNormalPerFace(IndexedFaceSet& ifs) {
  return
    (ifs.getNumberOfFaces() > 0) &&
    (ifs.getNormalBinding() == IndexedFaceSet::PB_PER_FACE);
}

bool SceneGraphProcessor::_hasNormalPerVertex(IndexedFaceSet& ifs) {
  return
    (ifs.getNumberOfCoord() > 0) &&
    (ifs.getNormalBinding() == IndexedFaceSet::PB_PER_VERTEX);
}

bool SceneGraphProcessor::_hasNormalPerCorner(IndexedFaceSet& ifs) {
  return
    (ifs.getNumberOfFaces() > 0) &&
    (ifs.getNormalBinding() == IndexedFaceSet::PB_PER_CORNER);
}

bool SceneGraphProcessor::hasIndexedFaceSetFaces() {
  return _hasIndexedFaceSetProperty(_hasFaces);
}

bool SceneGraphProcessor::hasIndexedFaceSetNormalNone() {
  return _hasIndexedFaceSetProperty(_hasNormalNone);
}

bool SceneGraphProcessor::hasIndexedFaceSetNormalPerFace() {
  return _hasIndexedFaceSetProperty(_hasNormalPerFace);
}

bool SceneGraphProcessor::hasIndexedFaceSetNormalPerVertex() {
  return _hasIndexedFaceSetProperty(_hasNormalPerVertex);
}

bool SceneGraphProcessor::hasIndexedFaceSetNormalPerCorner() {
  return _hasIndexedFaceSetProperty(_hasNormalPerCorner);
}

// VRML'97
//
// If the color field is not NULL, it shall contain a Color node, and
// the colours are applied to the line(s) as follows:
//
// If colorPerVertex is FALSE:
//
// If the colorIndex field is not empty, then one colour is used for
// each polyline of the IndexedLineSet. There must be at least as many
// indices in the colorIndex field as there are polylines in the
// IndexedLineSet. If the greatest index in the colorIndex field is N,
// then there must be N+1 colours in the Color node. The colorIndex
// field must not contain any negative entries.
//
// If the colorIndex field is empty, then the colours from the Color
// node are applied to each polyline of the IndexedLineSet in
// order. There must be at least as many colours in the Color node as
// there are polylines.
//
// If colorPerVertex is TRUE:
//
// If the colorIndex field is not empty, then colours are applied to
// each vertex of the IndexedLineSet in exactly the same manner that
// the coordIndex field is used to supply coordinates for each vertex
// from the Coordinate node. The colorIndex field must contain at
// least as many indices as the coordIndex field and must contain
// end-of-polyline markers (-1) in exactly the same places as the
// coordIndex field. If the greatest index in the colorIndex field is
// N, then there must be N+1 colours in the Color node.
//
// If the colorIndex field is empty, then the coordIndex field is used
// to choose colours from the Color node. If the greatest index in the
// coordIndex field is N, then there must be N+1 colours in the Color
// node.
//
// If the color field is NULL and there is a Material defined for the
// Appearance affecting this IndexedLineSet, the emissiveColor of the
// Material shall be used to draw the lines.

bool SceneGraphProcessor::_hasColorNone(IndexedLineSet& ils) {
  vector<float>& color = ils.getColor();
  return (color.size() == 0);
}

bool SceneGraphProcessor::_hasColorPerVertex(IndexedLineSet& ils) {
  vector<float>& color = ils.getColor();
  // vector<int>&   colorIndex    = ils.getColorIndex();
  bool           colorPerVerex = ils.getColorPerVertex();
  // not testing for errors, but
  // if(colorIndex.size()==0)
  //   we should have color.size()/3 == ils.getNumberOfCoord()
  // else
  //   we should have colorIndex.size() == ils.getNumberOfCoord()
  return color.size() > 0 && colorPerVerex == false;
}

bool SceneGraphProcessor::_hasColorPerPolyline(IndexedLineSet& ils) {
  vector<float>& color = ils.getColor();
  // vector<int>&   colorIndex    = ils.getColorIndex();
  bool           colorPerVerex = ils.getColorPerVertex();
  // not testing for errors, but
  // if(colorIndex.size()==0)
  //   we should have color.size()/3 == ils.getNumberOfPolylines()
  // else
  //   we should have colorIndex.size() == ils.getNumberOfPolylines()
  return color.size() > 0 && colorPerVerex == false;
}

bool SceneGraphProcessor::hasIndexedLineSetColorNone() {
  return _hasIndexedLineSetProperty(_hasColorNone);
}

bool SceneGraphProcessor::hasIndexedLineSetColorPerVertex() {
  return _hasIndexedLineSetProperty(_hasColorPerVertex);
}

bool SceneGraphProcessor::hasIndexedLineSetColorPerPolyline() {
  return _hasIndexedLineSetProperty(_hasColorPerPolyline);
}

bool SceneGraphProcessor::_hasEdges(Shape& shape) {
  return shape.nameEquals("EDGES");
}

bool SceneGraphProcessor::_hasIndexedFaceSetShown(Shape& shape) {
  return shape.hasGeometryIndexedFaceSet() && shape.getShow() == true;
}

bool SceneGraphProcessor::_hasIndexedFaceSetHidden(Shape& shape) {
  return shape.hasGeometryIndexedFaceSet() && shape.getShow() == false;
}

bool SceneGraphProcessor::_hasIndexedLineSetShown(Shape& shape) {
  return shape.hasGeometryIndexedLineSet() && shape.getShow() == true;
}

bool SceneGraphProcessor::_hasIndexedLineSetHidden(Shape& shape) {
  return shape.hasGeometryIndexedLineSet() && shape.getShow() == false;
}

bool SceneGraphProcessor::hasEdges() {
  return _hasShapeProperty(_hasEdges);
}

bool SceneGraphProcessor::hasIndexedFaceSetShown() {
  return _hasShapeProperty(_hasIndexedFaceSetShown);
}

bool SceneGraphProcessor::hasIndexedFaceSetHidden() {
  return _hasShapeProperty(_hasIndexedFaceSetHidden);
}

bool SceneGraphProcessor::hasIndexedLineSetShown() {
  return _hasShapeProperty(_hasIndexedLineSetShown);
}

bool SceneGraphProcessor::hasIndexedLineSetHidden() {
  return _hasShapeProperty(_hasIndexedLineSetHidden);
}

void SceneGraphProcessor::removeSceneGraphChild(const string& name) {
  vector<pNode>& children = _wrl.getChildren();
  vector<pNode>::iterator i;
  for (i = children.begin();i != children.end();i++)
    if ((*i)->nameEquals(name))
      break;
  if (i != children.end())
    children.erase(i);
}

void SceneGraphProcessor::pointsRemove() {
  removeSceneGraphChild("POINTS");
}

void SceneGraphProcessor::samplesRemove() {
  removeSceneGraphChild("SAMPLES");
}

void SceneGraphProcessor::planeRemove() {
  removeSceneGraphChild("PLANE");
}

void SceneGraphProcessor::surfaceRemove() {
  removeSceneGraphChild("SURFACE");
}

IndexedFaceSet* SceneGraphProcessor::_getNamedShapeIFS
(const string& name, bool create) {
  IndexedFaceSet* ifs = (IndexedFaceSet*)0;
  Node* node = _wrl.getChild(name);
  if (node != (Node*)0 && node->isShape()) {
    Shape* shape = (Shape*)node;
    node = shape->getGeometry();
    if (node != (Node*)0 && node->isIndexedFaceSet()) {
      ifs = (IndexedFaceSet*)node;
    }
  }
  if (ifs == (IndexedFaceSet*)0 && create) {
    Shape* shape = new Shape();
    shape->setName(name);
    _wrl.addChild(shape);
    Appearance* appearance = new Appearance();
    shape->setAppearance(appearance);
    Material* material = new Material();
    appearance->setMaterial(material);
    ifs = new IndexedFaceSet();
    shape->setGeometry(ifs);
  }
  return ifs;
}

//////////////////////////////////////////////////////////////////////
void meanFit(const vector<float>& coordPoints,
  const vector<float>& normalPoints,
  const Vec3f& min, const Vec3f& max, Vec4f& f) {

  int nPoints = (int)((coordPoints.size()) / 3);
  float x0 = min.x, x1 = max.x, dx = x1 - x0;
  float y0 = min.y, y1 = max.y, dy = y1 - y0;
  float z0 = min.z, z1 = max.z, dz = z1 - z0;
  float dMax = dx; if (dy > dMax) dMax = dy; if (dz > dMax) dMax = dz;

  // compute the mean of the points contained in the 
  double x, y, z, nn;
  double xMean = 0.0f;
  double yMean = 0.0f;
  double zMean = 0.0f;
  double nxMean = 0.0f;
  double nyMean = 0.0f;
  double nzMean = 0.0f;
  int    nMean = 0;
  for (int iP = 0;iP < nPoints;iP++) {
    x = (double)(coordPoints[3 * iP]);
    y = (double)(coordPoints[3 * iP + 1]);
    z = (double)(coordPoints[3 * iP + 2]);
    if (x0 <= x && x <= x1 && y0 <= y && y <= y1 && z0 <= z && z <= z1) {
      xMean += x;
      yMean += y;
      zMean += z;
      x = (double)(normalPoints[3 * iP]);
      y = (double)(normalPoints[3 * iP + 1]);
      z = (double)(normalPoints[3 * iP + 2]);
      nxMean += x;
      nyMean += y;
      nzMean += z;
      nMean++;
    }
  }

  if (nMean == 0) {
    // no points found inside the bounding box
    // throw exception ??
    f.x = f.y = f.z = f.w = 0.0f;
    return;
  }

  // normalize the point mean
  xMean /= ((double)nMean);
  yMean /= ((double)nMean);
  zMean /= ((double)nMean);
  // normalize the normal mean to unit length
  nn = nxMean * nxMean + nyMean * nyMean + nzMean * nzMean;
  if (nn > 0.0) {
    nn = sqrt(nn); nxMean /= nn; nyMean /= nn; nzMean /= nn;
  }
  // set the linear function coefficients
  f.x = (float)(nxMean);
  f.y = (float)(nyMean);
  f.z = (float)(nzMean);
  f.w = -(float)(nxMean * xMean + nyMean * yMean + nzMean * zMean);
}

//////////////////////////////////////////////////////////////////////
void SceneGraphProcessor::fitSinglePlane
(const Vec3f& center, const Vec3f& size,
  const float scale, const bool cube, Vec4f& f) {

  // make sure that the bounding box is not empty
  if (size.x <= 0.0f || size.y <= 0.0f || size.z <= 0.0f || scale <= 0.0f) {
    // throw exception ?
    return;
  }

  // find the input point set in the scene graph
  IndexedFaceSet* points = _getNamedShapeIFS("POINTS", false);
  if (points == (IndexedFaceSet*)0) return;
  if (points->getNormalBinding() != IndexedFaceSet::PB_PER_VERTEX) return;
  const vector<float>& coordPoints = points->getCoord();
  const vector<float>& normalPoints = points->getNormal();
  // int nPoints = (int)(coordPoints.size()/3);
  // cout << "  nPoints = " << nPoints << endl;

  // compute the coordinates of the bounding box corners
  float dx = size.x / 2.0f, dy = size.y / 2.0f, dz = size.z / 2.0f;
  float dMax = dx; if (dy > dMax) dMax = dy; if (dz > dMax) dMax = dz;
  if (cube) { dx = dy = dz = dMax; }
  if (scale > 0.0f) { dx *= scale; dy *= scale; dz *= scale; }
  float x0 = center.x - dx; float y0 = center.y - dy; float z0 = center.z - dz;
  float x1 = center.x + dx; float y1 = center.y + dy; float z1 = center.z + dz;
  Vec3f min(x0, y0, z0);
  Vec3f max(x1, y1, z1);
  float x, y, z;
  Vec3f v[8]; // bbox corner coordinates
  for (int i = 0;i < 8;i++) {
    v[i].z = (((i >> 0) & 0x1) == 0) ? z0 : z1;
    v[i].y = (((i >> 1) & 0x1) == 0) ? y0 : y1;
    v[i].x = (((i >> 2) & 0x1) == 0) ? x0 : x1;
  }

  // fit a plane to the points contained in the bounding box
  // eigenFit(coordPoints,normalPoints,min,max,f);
  meanFit(coordPoints, normalPoints, min, max, f);

  // find or create the output surface in the scene graph 
  IndexedFaceSet* plane = _getNamedShapeIFS("PLANE", true);
  plane->clear();
  vector<float>& coordIfs = plane->getCoord();

  // evaluate the linear function at bounding box corners
  float F[8]; // function values at bbox corners
  bool  b[8]; // function is positive or negative ?
  for (int i = 0;i < 8;i++) {
    x = v[i].x; y = v[i].y; z = v[i].z;
    b[i] = (F[i] = x * f.x + y * f.y + z * f.z + f.w) < 0.0f;
  }

  // cout << "//                6 ----- 7 = (x1,y1,z1)" << endl;
  // cout << "//               /|      /|             " << endl;
  // cout << "//              4 ----- 5 |             " << endl;
  // cout << "//              | |     | |             " << endl;
  // cout << "//              | 2 ----| 3             " << endl;
  // cout << "//              |/      |/              " << endl;
  // cout << "// (x0,y0,z0) = 0 ----- 1               " << endl;

  //////////////////////////////////////////////////////////////////////
  //
  //    vertices      //    edges                 //    faces
  //      6-----7     //        [6]---11---[7]    //        1
  //     /|    /|     //        /|         /|     //        | 3
  //    4-----5 |     //       6 2        7 3     //        |/
  //    | 2---|-3     //      /  |       /  |     //    4---+---5
  //    |/    |/      //    [4]---10---[5]  |     //       /|
  //    0-----1       //     |   |      |   |     //      2 |
  //                  //     |  [2]--9--|--[3]    //        0
  //                  //     0  /       1  /      //
  //                  //     | 4        | 5       //
  //                  //     |/         |/        //
  //                  //    [0]---8----[1]        //
  //

  const int (*edge)[2] = IsoSurf::getEdgeTable();

  // compute the isovertex coordinates
  float tj, tk;
  int   iE[12], iV, i, j, k;
  for (i = 0;i < 12;i++) {
    iV = -1;
    j = edge[i][0];
    k = edge[i][1];
    if (b[j] != b[k]) {
      // create new vertex index
      iV = (int)((coordIfs.size() / 3));
      // isovertex coordinates
      tk = F[j] / (F[j] - F[k]);
      tj = F[k] / (F[k] - F[j]);
      x = tj * v[j].x + tk * v[k].x;
      y = tj * v[j].y + tk * v[k].y;
      z = tj * v[j].z + tk * v[k].z;
      coordIfs.push_back(x);
      coordIfs.push_back(y);
      coordIfs.push_back(z);
    }
    iE[i] = iV;
  }

  // create isosurface faces
  vector<int>& coordIndex = plane->getCoordIndex(); // coordIndex.size()==0
  /* int nFaces = */ IsoSurf::makeCellFaces(b, iE, coordIndex);
  // cout << "  nFaces = " << nFaces << endl;

  // save plane normal vector as face normal
  plane->setNormalPerVertex(false);
  vector<float>& normal = plane->getNormal();
  // assert(normal.size()==0);
  normal.push_back(f.x);
  normal.push_back(f.y);
  normal.push_back(f.z);
}

void SceneGraphProcessor::_fitPlane
(const Vec3f& center, const Vec3f& size,
  const float scale, const bool cube, Vec4f& f,
  const vector<float>& coordPoints, const vector<float>& normalPoints, float* F, bool* b,
  bool showPlane) {

  // make sure that the bounding box is not empty
  if (size.x <= 0.0f || size.y <= 0.0f || size.z <= 0.0f || scale <= 0.0f) {
    // throw exception ?
    return;
  }

  // find the input point set in the scene graph
  //IndexedFaceSet* points  = _getNamedShapeIFS("POINTS",false);

  // int nPoints = (int)(coordPoints.size()/3);
  // cout << "  nPoints = " << nPoints << endl;

  // compute the coordinates of the bounding box corners
  float dx = size.x / 2.0f, dy = size.y / 2.0f, dz = size.z / 2.0f;
  float dMax = dx; if (dy > dMax) dMax = dy; if (dz > dMax) dMax = dz;
  if (cube) { dx = dy = dz = dMax; }
  if (scale > 0.0f) { dx *= scale; dy *= scale; dz *= scale; }
  float x0 = center.x - dx; float y0 = center.y - dy; float z0 = center.z - dz;
  float x1 = center.x + dx; float y1 = center.y + dy; float z1 = center.z + dz;
  Vec3f min(x0, y0, z0);
  Vec3f max(x1, y1, z1);
  float x, y, z;
  Vec3f v[8]; // bbox corner coordinates
  for (int i = 0;i < 8;i++) {
    v[i].z = (((i >> 0) & 0x1) == 0) ? z0 : z1;
    v[i].y = (((i >> 1) & 0x1) == 0) ? y0 : y1;
    v[i].x = (((i >> 2) & 0x1) == 0) ? x0 : x1;
  }

  // fit a plane to the points contained in the bounding box
  // eigenFit(coordPoints,normalPoints,min,max,f);
  meanFit(coordPoints, normalPoints, min, max, f);

  // find or create the output surface in the scene graph


  // evaluate the linear function at bounding box corners
  if (F == nullptr) {
    F = new float[8]; // function values at bbox corners
  }
  if (b == nullptr) {
    b = new bool[8]; // function is positive or negative ?
  }
  // F.clear();
  // b.clear();
  // F.resize(8);
  // b.resize(8);
  // float F[8]; // function values at bbox corners
  // bool  b[8]; // function is positive or negative ?
  for (int i = 0;i < 8;i++) {
    x = v[i].x; y = v[i].y; z = v[i].z;
    b[i] = (F[i] = x * f.x + y * f.y + z * f.z + f.w) < 0.0f;
  }

  // cout << "//                6 ----- 7 = (x1,y1,z1)" << endl;
  // cout << "//               /|      /|             " << endl;
  // cout << "//              4 ----- 5 |             " << endl;
  // cout << "//              | |     | |             " << endl;
  // cout << "//              | 2 ----| 3             " << endl;
  // cout << "//              |/      |/              " << endl;
  // cout << "// (x0,y0,z0) = 0 ----- 1               " << endl;

  //////////////////////////////////////////////////////////////////////
  //
  //    vertices      //    edges                 //    faces
  //      6-----7     //        [6]---11---[7]    //        1
  //     /|    /|     //        /|         /|     //        | 3
  //    4-----5 |     //       6 2        7 3     //        |/
  //    | 2---|-3     //      /  |       /  |     //    4---+---5
  //    |/    |/      //    [4]---10---[5]  |     //       /|
  //    0-----1       //     |   |      |   |     //      2 |
  //                  //     |  [2]--9--|--[3]    //        0
  //                  //     0  /       1  /      //
  //                  //     | 4        | 5       //
  //                  //     |/         |/        //
  //                  //    [0]---8----[1]        //
  //

  if (!showPlane) {
    return;
  }

  const int (*edge)[2] = IsoSurf::getEdgeTable();
  IndexedFaceSet* plane = _getNamedShapeIFS("PLANE", false);
  vector<float>& coordIfs = plane->getCoord();

  // compute the isovertex coordinates
  float tj, tk;
  int   iE[12], iV, i, j, k;
  for (i = 0;i < 12;i++) {
    iV = -1;
    j = edge[i][0];
    k = edge[i][1];
    if (b[j] != b[k]) {
      // create new vertex index
      iV = (int)((coordIfs.size() / 3));
      // isovertex coordinates
      tk = F[j] / (F[j] - F[k]);
      tj = F[k] / (F[k] - F[j]);
      x = tj * v[j].x + tk * v[k].x;
      y = tj * v[j].y + tk * v[k].y;
      z = tj * v[j].z + tk * v[k].z;
      coordIfs.push_back(x);
      coordIfs.push_back(y);
      coordIfs.push_back(z);
    }
    iE[i] = iV;
  }

  // create isosurface faces
  vector<int>& coordIndex = plane->getCoordIndex();


  int nIndexBefore = (int)coordIndex.size();
  IsoSurf::makeCellFaces(b, iE, coordIndex);
  int nFaces = ((int)coordIndex.size() - nIndexBefore) / 4;

  // save plane normal vector as face normal
  if (nFaces > 0) {
    plane->setNormalPerVertex(false);
    vector<float>& normal = plane->getNormal();
    for (int k = 0; k < nFaces; k++) {
      normal.push_back(f.x);
      normal.push_back(f.y);
      normal.push_back(f.z);
    }
  }
}

//////////////////////////////////////////////////////////////////////
void SceneGraphProcessor::fitMultiplePlanes
(const Vec3f& center, const Vec3f& size,
  const int depth, const float scale, const bool cube,
  vector<float>& f) {
  // make sure that the bounding box is not empty
  if (size.x <= 0.0f || size.y <= 0.0f || size.z <= 0.0f || scale <= 0.0f) {
    // throw exception ?
    return;
  }

  IndexedFaceSet* planes = _getNamedShapeIFS("SURFACE", true);
  planes->clear();

  // find the input point set in the scene graph
  IndexedFaceSet* points = _getNamedShapeIFS("POINTS", false);
  if (points == (IndexedFaceSet*)0) return;
  if (points->getNormalBinding() != IndexedFaceSet::PB_PER_VERTEX) return;
  vector<float>& coordPoints = points->getCoord();
  const vector<float>& normalPoints = points->getNormal();
  // int nPoints = (int)(coordPoints.size()/3);
  // cout << "  nPoints = " << nPoints << endl;

  // compute the coordinates of the bounding box corners to create the partition
  float dx = size.x / 2.0f, dy = size.y / 2.0f, dz = size.z / 2.0f;
  float dMax = dx; if (dy > dMax) dMax = dy; if (dz > dMax) dMax = dz;
  if (cube) { dx = dy = dz = dMax; }
  if (scale > 0.0f) { dx *= scale; dy *= scale; dz *= scale; }
  float x0 = center.x - dx; float y0 = center.y - dy; float z0 = center.z - dz;
  float x1 = center.x + dx; float y1 = center.y + dy; float z1 = center.z + dz;
  Vec3f min(x0, y0, z0);
  Vec3f max(x1, y1, z1);

  _createPartition(min, max, depth, coordPoints);

  int N = _nGrid; // N = 2^depth


  vector<float>& coordIfs = planes->getCoord();
  vector<float> coordCell;
  vector<float> normalCell;

  // Compute all corners values
  for (int ix = 0; ix < N; ix++) {
    for (int iy = 0; iy < N; iy++) {
      for (int iz = 0; iz < N; iz++) {
        Vec3f minCell, maxCell;
        // compute the bounding box of the current cell
        minCell.x = ((N - ix) * min.x + (ix)*max.x) / (N);
        maxCell.x = ((N - ix - 1) * min.x + (ix + 1) * max.x) / (N);
        minCell.y = ((N - iy) * min.y + (iy)*max.y) / (N);
        maxCell.y = ((N - iy - 1) * min.y + (iy + 1) * max.y) / (N);
        minCell.z = ((N - iz) * min.z + (iz)*max.z) / (N);
        maxCell.z = ((N - iz - 1) * min.z + (iz + 1) * max.z) / (N);

        int iCell = ix + N * (iy + N * iz);
        int iPoint = _first[iCell];

        // skip empty cells
        if (iPoint == -1)
          continue;

        coordCell.clear();
        normalCell.clear();
        // Obtain the points contained in the current cell
        for (int iP = iPoint; iP != -1; iP = _next[iP]) {
          float x = coordPoints[3 * iP];
          float y = coordPoints[3 * iP + 1];
          float z = coordPoints[3 * iP + 2];
          coordCell.push_back(x);
          coordCell.push_back(y);
          coordCell.push_back(z);

          if (normalPoints.size() > 0) {
            x = normalPoints[3 * iP];
            y = normalPoints[3 * iP + 1];
            z = normalPoints[3 * iP + 2];

            normalCell.push_back(x);
            normalCell.push_back(y);
            normalCell.push_back(z);
          }
        }


        // fit a plane to the points contained in the bounding box
        Vec4f fCell;
        meanFit(coordCell, normalCell, minCell, maxCell, fCell);


        // Now, for the current linear function, create the isosurface in the current cell
        float x, y, z;
        Vec3f v[8]; // bbox corner coordinates
        for (int i = 0;i < 8;i++) {
          v[i].z = (((i >> 0) & 0x1) == 0) ? minCell.z : maxCell.z;
          v[i].y = (((i >> 1) & 0x1) == 0) ? minCell.y : maxCell.y;
          v[i].x = (((i >> 2) & 0x1) == 0) ? minCell.x : maxCell.x;
        }



        // evaluate the linear function at bounding box corners
        float F[8]; // function values at bbox corners
        bool  b[8]; // function is positive or negative ?
        for (int i = 0;i < 8;i++) {
          x = v[i].x; y = v[i].y; z = v[i].z;
          b[i] = (F[i] = x * fCell.x + y * fCell.y + z * fCell.z + fCell.w) < 0.0f;
        }

        // cout << "//                6 ----- 7 = (x1,y1,z1)" << endl;
        // cout << "//               /|      /|             " << endl;
        // cout << "//              4 ----- 5 |             " << endl;
        // cout << "//              | |     | |             " << endl;
        // cout << "//              | 2 ----| 3             " << endl;
        // cout << "//              |/      |/              " << endl;
        // cout << "// (x0,y0,z0) = 0 ----- 1               " << endl;

        //////////////////////////////////////////////////////////////////////
        //
        //    vertices      //    edges                 //    faces
        //      6-----7     //        [6]---11---[7]    //        1
        //     /|    /|     //        /|         /|     //        | 3
        //    4-----5 |     //       6 2        7 3     //        |/
        //    | 2---|-3     //      /  |       /  |     //    4---+---5
        //    |/    |/      //    [4]---10---[5]  |     //       /|
        //    0-----1       //     |   |      |   |     //      2 |
        //                  //     |  [2]--9--|--[3]    //        0
        //                  //     0  /       1  /      //
        //                  //     | 4        | 5       //
        //                  //     |/         |/        //
        //                  //    [0]---8----[1]        //
        //

        const int (*edge)[2] = IsoSurf::getEdgeTable();

        // compute the isovertex coordinates
        float tj, tk;
        int   iE[12], iV, i, j, k;
        for (i = 0;i < 12;i++) {
          iV = -1;
          j = edge[i][0];
          k = edge[i][1];
          if (b[j] != b[k]) {
            // create new vertex index
            iV = (int)((coordIfs.size() / 3));
            // isovertex coordinates
            tk = F[j] / (F[j] - F[k]);
            tj = F[k] / (F[k] - F[j]);
            x = tj * v[j].x + tk * v[k].x;
            y = tj * v[j].y + tk * v[k].y;
            z = tj * v[j].z + tk * v[k].z;
            coordIfs.push_back(x);
            coordIfs.push_back(y);
            coordIfs.push_back(z);
          }
          iE[i] = iV;
        }

        // create isosurface faces
        vector<int>& coordIndex = planes->getCoordIndex(); // coordIndex.size()==0
        /* int nFaces = */ IsoSurf::makeCellFaces(b, iE, coordIndex);
        // cout << "  nFaces = " << nFaces << endl;

      }
    }
  }

  // Recompute normals
  planes->setNormalPerVertex(false);
  _computeNormalPerFace(*planes);


  _deletePartition();
}


//////////////////////////////////////////////////////////////////////
void SceneGraphProcessor::fitContinuous
(const Vec3f& center, const Vec3f& size,
  const int depth, const float scale, const bool cube,
  vector<float>& fGrid) {

  IndexedFaceSet* planes = _getNamedShapeIFS("SURFACE", true);
  planes->clear();

  // find the input point set in the scene graph
  IndexedFaceSet* points = _getNamedShapeIFS("POINTS", false);
  if (points == (IndexedFaceSet*)0) return;
  if (points->getNormalBinding() != IndexedFaceSet::PB_PER_VERTEX) return;
  vector<float>& coordPoints = points->getCoord();
  const vector<float>& normalPoints = points->getNormal();
  // int nPoints = (int)(coordPoints.size()/3);
  // cout << "  nPoints = " << nPoints << endl;

  // compute the coordinates of the bounding box corners to create the partition
  float dx = size.x / 2.0f, dy = size.y / 2.0f, dz = size.z / 2.0f;
  float dMax = dx; if (dy > dMax) dMax = dy; if (dz > dMax) dMax = dz;
  if (cube) { dx = dy = dz = dMax; }
  if (scale > 0.0f) { dx *= scale; dy *= scale; dz *= scale; }
  float x0 = center.x - dx; float y0 = center.y - dy; float z0 = center.z - dz;
  float x1 = center.x + dx; float y1 = center.y + dy; float z1 = center.z + dz;
  Vec3f min(x0, y0, z0);
  Vec3f max(x1, y1, z1);

  _createPartition(min, max, depth, coordPoints);

  int N = _nGrid; // N = 2^depth
  int nGridPoints = (N + 1) * (N + 1) * (N + 1);
  int nCells = N * N * N;



  vector<float> coordCell;
  vector<float> normalCell;

  // This vector maps each cell corner to the global grid vertex
  vector<int> cornerToVertex(8 * nCells, -1);
  // FGrid holds the corresponding function values at grid corners
  vector<float> FGrid(8 * nCells);

  Vec3f minCell, maxCell;

  vector<float> fCells(4 * nCells, 0.0f);

  for (int ix = 0; ix < N; ix++) {
    for (int iy = 0; iy < N; iy++) {
      for (int iz = 0; iz < N; iz++) {

        // compute the bounding box of the current cell
        minCell.x = ((N - ix) * min.x + (ix)*max.x) / (N);
        maxCell.x = ((N - ix - 1) * min.x + (ix + 1) * max.x) / (N);
        minCell.y = ((N - iy) * min.y + (iy)*max.y) / (N);
        maxCell.y = ((N - iy - 1) * min.y + (iy + 1) * max.y) / (N);
        minCell.z = ((N - iz) * min.z + (iz)*max.z) / (N);
        maxCell.z = ((N - iz - 1) * min.z + (iz + 1) * max.z) / (N);

        int iCell = ix + N * (iy + N * iz);
        int iPoint = _first[iCell];
        // skip empty cells
        if (iPoint == -1)
          continue;

        coordCell.clear();
        normalCell.clear();
        // Obtain the points contained in the current cell
        for (int iP = iPoint; iP != -1; iP = _next[iP]) {
          float x = coordPoints[3 * iP];
          float y = coordPoints[3 * iP + 1];
          float z = coordPoints[3 * iP + 2];
          coordCell.push_back(x);
          coordCell.push_back(y);
          coordCell.push_back(z);

          if (normalPoints.size() > 0) {
            x = normalPoints[3 * iP];
            y = normalPoints[3 * iP + 1];
            z = normalPoints[3 * iP + 2];

            normalCell.push_back(x);
            normalCell.push_back(y);
            normalCell.push_back(z);
          }
        }


        // fit a plane to the points contained in the bounding box
        Vec4f fCell;
        meanFit(coordCell, normalCell, minCell, maxCell, fCell);

        float x, y, z;
        Vec3f v[8]; // bbox corner coordinates
        for (int i = 0;i < 8;i++) {
          v[i].z = (((i >> 0) & 0x1) == 0) ? minCell.z : maxCell.z;
          v[i].y = (((i >> 1) & 0x1) == 0) ? minCell.y : maxCell.y;
          v[i].x = (((i >> 2) & 0x1) == 0) ? minCell.x : maxCell.x;
          // map cell corner to global grid vertex
          int ixV = ix + ((i >> 2) & 0x1);;
          int iyV = iy + ((i >> 1) & 0x1);;
          int izV = iz + ((i >> 0) & 0x1);
          int iVGlobal = ixV + (N + 1) * (iyV + (N + 1) * izV);
          cornerToVertex[iCell * 8 + i] = iVGlobal;
        }



        // evaluate the linear function at bounding box corners
        float F[8]; // function values at bbox corners
        bool  b[8]; // function is positive or negative ?
        for (int i = 0;i < 8;i++) {
          x = v[i].x; y = v[i].y; z = v[i].z;
          b[i] = (F[i] = x * fCell.x + y * fCell.y + z * fCell.z + fCell.w) < 0.0f;
          // store the function value at the grid corner
          FGrid[iCell * 8 + i] = F[i];
        }

        // Finally, save the linear function coefficients
        fGrid.push_back(fCell.x);
        fGrid.push_back(fCell.y);
        fGrid.push_back(fCell.z);
        fGrid.push_back(fCell.w);
      }
    }
  }





  const int (*edge)[2] = IsoSurf::getEdgeTable();
  vector<float>& coordIfs = planes->getCoord();
  vector<int>& coordIndex = planes->getCoordIndex();

  int nGridVerts = (N + 1) * (N + 1) * (N + 1);
  vector<float> globalSum(nGridVerts, 0.0f);
  vector<int>   globalCount(nGridVerts, 0);

  // Accumulate F values at grid vertices
  for (int ix = 0; ix < N; ix++) {
    for (int iy = 0; iy < N; iy++) {
      for (int iz = 0; iz < N; iz++) {

        int iCell = ix + N * (iy + N * iz);

        // skip empty cells
        if (_first[iCell] == -1)
          continue;

        for (int k = 0; k < 8; k++) {
          int iVCell = cornerToVertex[iCell * 8 + k];
          globalSum[iVCell] += FGrid[iCell * 8 + k];
          globalCount[iVCell] += 1;
        }
      }
    }
  }

  // compute the mean F values at grid vertices
  for (int i = 0; i < nGridVerts; i++) {
    if (globalCount[i] > 0) {
      globalSum[i] /= (float)globalCount[i];
    }
    else {
      globalSum[i] = 0;
    }
  }




  // Finally, create the isosurface in each cell using the mean F values
  for (int ix = 0; ix < N; ix++) {
    for (int iy = 0; iy < N; iy++) {
      for (int iz = 0; iz < N; iz++) {

        int iCell = ix + N * (iy + N * iz);
        // skip empty cells
        if (_first[iCell] == -1)
          continue;

        // Recompute the bounding box of the current cell
        minCell.x = ((N - ix) * min.x + (ix)*max.x) / (N);
        maxCell.x = ((N - ix - 1) * min.x + (ix + 1) * max.x) / (N);
        minCell.y = ((N - iy) * min.y + (iy)*max.y) / (N);
        maxCell.y = ((N - iy - 1) * min.y + (iy + 1) * max.y) / (N);
        minCell.z = ((N - iz) * min.z + (iz)*max.z) / (N);
        maxCell.z = ((N - iz - 1) * min.z + (iz + 1) * max.z) / (N);

        float x, y, z;
        Vec3f v[8]; // bbox corner coordinates
        for (int i = 0;i < 8;i++) {
          v[i].z = (((i >> 0) & 0x1) == 0) ? minCell.z : maxCell.z;
          v[i].y = (((i >> 1) & 0x1) == 0) ? minCell.y : maxCell.y;
          v[i].x = (((i >> 2) & 0x1) == 0) ? minCell.x : maxCell.x;
        }


        float F[8];
        bool  b[8];
        // Retrieve the mean F values at cell corners
        for (int k = 0; k < 8; k++) {
          int iVCell = cornerToVertex[iCell * 8 + k];
          F[k] = globalSum[iVCell];
          b[k] = (globalSum[iVCell] < 0.0f);
        }

        // compute the isovertex coordinates
        float tj, tk;
        int   iE[12], iV, i, j, k;
        for (i = 0;i < 12;i++) {
          iV = -1;
          j = edge[i][0];
          k = edge[i][1];
          if (b[j] != b[k]) {
            // create new vertex index
            iV = (int)((coordIfs.size() / 3));
            // isovertex coordinates
            tk = F[j] / (F[j] - F[k]);
            tj = F[k] / (F[k] - F[j]);
            float x = tj * v[j].x + tk * v[k].x;
            float y = tj * v[j].y + tk * v[k].y;
            float z = tj * v[j].z + tk * v[k].z;
            coordIfs.push_back(x);
            coordIfs.push_back(y);
            coordIfs.push_back(z);
          }
          iE[i] = iV;
        }


        int nIndexBefore = (int)coordIndex.size();
        IsoSurf::makeCellFaces(b, iE, coordIndex);
        int nFaces = ((int)coordIndex.size() - nIndexBefore) / 4;

        
      }
    }
  }


  
  // Recompute normals
  planes->setNormalPerVertex(false);
  _computeNormalPerFace(*planes);

  _deletePartition();



}

//////////////////////////////////////////////////////////////////////
// implementation of the following methods
// is OPTIONAL for extra credit

//////////////////////////////////////////////////////////////////////
void SceneGraphProcessor::fitWatertight
(const Vec3f& center, const Vec3f& size,
  const int depth, const float scale, const bool isCube,
  vector<float>& fGrid) {

  std::cerr << "SceneGraphProcessor::fitWatertight() {" << endl;

  if (size.x <= 0.0f || size.y <= 0.0f || size.z <= 0.0f || scale <= 0.0f) return;

  // TODO
  // More details in the file DG2025-FCEN-TP5-Notes.pdf

  // 1) get the scene graph node named "POINTS"
  IndexedFaceSet* points = _getNamedShapeIFS("POINTS", false);

  // 2) if there is no such node return without doing anything
  if (points == (IndexedFaceSet*)0) return;

  // 3) if the node is found but it is empty, or it does not have
  // normals per vertex, also return without doing anything
  if (points->getNormalBinding() != IndexedFaceSet::PB_PER_VERTEX) return;

  // 4) get the coord and normal vectors from the points node 
  vector<float>& coordPoints = points->getCoord();
  const vector<float>& normalPoints = points->getNormal();

  // 5) from the center, size, and scale arguments, compute the min &
  // max corners of the bounding box
  // Vec3f min;
  // Vec3f max;
  float dx = size.x / 2.0f, dy = size.y / 2.0f, dz = size.z / 2.0f;
  float dMax = dx; if (dy > dMax) dMax = dy; if (dz > dMax) dMax = dz;
  if (isCube) { dx = dy = dz = dMax; }
  if (scale > 0.0f) { dx *= scale; dy *= scale; dz *= scale; }
  float x0 = center.x - dx; float y0 = center.y - dy; float z0 = center.z - dz;
  float x1 = center.x + dx; float y1 = center.y + dy; float z1 = center.z + dz;
  Vec3f min(x0, y0, z0);
  Vec3f max(x1, y1, z1);

  // 6) create a partition of the points as an array of linked lists
  _createPartition(min, max, depth, coordPoints);

  // 7) determine the total number of grid vertices as a function of
  // _nGrid 
  int N = _nGrid; // N = 2^depth
  int nGridVertices = (N + 1) * (N + 1) * (N + 1);

  // 7) initialize fGrid with zeros
  // fGrid.clear();
  // fGrid.insert(fGrid.end(),nGridVertices,0.0f);
  fGrid.clear();
  fGrid.resize(nGridVertices, 0.0f);

  // 8) create a temporary array of the same size to accumulate weights
  // vector<float> wGrid;
  // wGrid.insert(wGrid.end(),nGridVertices,0.0f);
  vector<float> wGrid(nGridVertices, 0.0f);

  // 9) allocate arrays to do the local linear fit
  // Vec4f fCell;
  // Vec3f minCell,maxCell;
  // vector<float> coordCell;
  // vector<float> normalCell;
  Vec4f fCell;
  Vec3f minCell, maxCell;
  vector<float> coordCell;
  vector<float> normalCell;

  // 10) accumulate grid vertex funtion values
  N = _nGrid;
  Vec3f nCell;
  int nonEmptyCells = 0;

  for (int ix = 0; ix < N; ix++) {
    for (int iy = 0; iy < N; iy++) {
      for (int iz = 0; iz < N; iz++) {
        // compute the bounding box of the current cell
        minCell.x = ((N - ix) * min.x + (ix)*max.x) / (N);
        maxCell.x = ((N - ix - 1) * min.x + (ix + 1) * max.x) / (N);
        minCell.y = ((N - iy) * min.y + (iy)*max.y) / (N);
        maxCell.y = ((N - iy - 1) * min.y + (iy + 1) * max.y) / (N);
        minCell.z = ((N - iz) * min.z + (iz)*max.z) / (N);
        maxCell.z = ((N - iz - 1) * min.z + (iz + 1) * max.z) / (N);

        int iCell = ix + N * (iy + N * iz);
        int iPoint = _first[iCell];
        // skip empty cells
        if (iPoint == -1)
          continue;

        // 11) copy coord and normal values of points in cell onto cell arrays
        coordCell.clear();
        normalCell.clear();
        nCell.x = nCell.y = nCell.z = N;

        for (int iP = iPoint; iP != -1; iP = _next[iP]) {
          float x = coordPoints[3 * iP];
          float y = coordPoints[3 * iP + 1];
          float z = coordPoints[3 * iP + 2];
          coordCell.push_back(x);
          coordCell.push_back(y);
          coordCell.push_back(z);

          if (normalPoints.size() > 0) {
            x = normalPoints[3 * iP];
            y = normalPoints[3 * iP + 1];
            z = normalPoints[3 * iP + 2];

            normalCell.push_back(x);
            normalCell.push_back(y);
            normalCell.push_back(z);
          }
          //nPoints++;
        }

        // 12) fit a linear function to the points contained in the
        // cell using the meanFit() method
        meanFit(coordCell, normalCell, minCell, maxCell, fCell);


        // 13) evaluate the local linear function at each of the
        // cell's eight vertices, accumulate the values in the
        // corresponding entries of the fGrid vector, and increment
        // the corresponding entry of the weights vector wGrid by 1
        for (int k = 0;k < 8;k++) {
          int ixV = ix + ((k >> 2) & 0x1);;
          int iyV = iy + ((k >> 1) & 0x1);;
          int izV = iz + ((k >> 0) & 0x1);
          int iV = ixV + (N + 1) * (iyV + (N + 1) * izV);
          float x = ((k & 0x4) == 0) ? minCell.x : maxCell.x;
          float y = ((k & 0x2) == 0) ? minCell.y : maxCell.y;
          float z = ((k & 0x1) == 0) ? minCell.z : maxCell.z;
          float fValue = x * fCell.x + y * fCell.y + z * fCell.z + fCell.w;
          fGrid[iV] += fValue;
          wGrid[iV] += 1.0f;


        }
      }
    }
  }


  // 14) arrays needded to implement the wavefront propagation algorithm
  vector<int> src;
  vector<int> dst;

  // 15) normalize the fGrid values
  for (int iV = 0;iV < nGridVertices;iV++)
    if (wGrid[iV] > 0.0f) { // only for vertices of occupied cells !
      fGrid[iV] /= wGrid[iV];
      // since the wGrid[iV] value is no longer needed, we will use it
      // to indicate which vertices have defined function values
      // (WGrid[iV]==-2), and which ones have undefined values
      // (wGrid[iV]==0)
      wGrid[iV] = -2.0f;
      // save the grid vertex indices of occupied cells in a list 
      src.push_back(iV);
    }

  // 16) extend vertex function values to all vertices by wavefront propagation
  int iV0, iV1, iV, ix, iy, iz;
  float fV0;

  // for each vertex in the wavefront
  while (src.size() > 0) {
    iV0 = iV = src.back(); src.pop_back(); // wGrid[iV0]<0
    // 17) compute the x,y,z integer coordinates of the vertex
    // inverting the formula
    // iV = ix+(N+1)*(iy+(N+1)*iz);
    iz = iV0 / ((N + 1) * (N + 1));
    iy = (iV0 - iz * (N + 1) * (N + 1)) / (N + 1);
    ix = iV0 - iz * (N + 1) * (N + 1) - iy * (N + 1);

    // 18) get the vertex function falue from fGrid
    fV0 = fGrid[iV0];

    // 19) for each of the 6 neighbors jV of vertex iV in the grid,
    // if vertex jV is inside the grid, and it does not have a
    // defined function value yet (i.e. wGrid[jV]>=0), then


    for (int k = 0;k < 6;k++) {
      // neighboors in x
      int ixV = ix + ((k == 0) ? -1 : (k == 1) ? 1 : 0);
      // neighboors in y
      int iyV = iy + ((k == 2) ? -1 : (k == 3) ? 1 : 0);
      // neighboors in z
      int izV = iz + ((k == 4) ? -1 : (k == 5) ? 1 : 0);

      if (ixV >= 0 && ixV <= N && iyV >= 0 && iyV <= N && izV >= 0 && izV <= N) {
        int jV = ixV + (N + 1) * (iyV + (N + 1) * izV);


        // The vertex doesn't have a defined function value yet
        if (abs(wGrid[jV]) < 1e-6) {
          // 19.1) add the fGrid iV function value to the fGrid jV function
          // value
          fGrid[jV] += fV0;
          // 19.2) increment the wGrid jV weight value by one
          wGrid[jV] += 1.0f;
          // 19.3) if this is the first visit to this grid vertex (i.e. if
          // the wGrid jV value is equal to zero), save the jV index in
          // the output wavefront list dst
          if (abs(wGrid[jV] - 1.0f) < 1e-6)
            dst.push_back(jV);
        }
      }

    }
    // 20) note that at this point the src list is empty, and the dst
    // list in general is not empty

    // 20) normalize new function values for the grid vertices in the
    // new wavefront, and create new wavefront
    while (dst.size() > 0) {
      iV1 = dst.back(); dst.pop_back();
      // if(wGrid[iV1]>0.0f) {
      fGrid[iV1] /= wGrid[iV1];
      // we use -2 to indicate the vertices of the occupied cells, and
      // -1 the vertices of cells where the function has been defined
      // by the wavefront propagation algorithm
      wGrid[iV1] = -1.0f;
      src.push_back(iV1);
      // }
    } // while(dst.size()>0)

  } // while(source.size()>0)

  // create output surface
  IndexedFaceSet* surface = _getNamedShapeIFS("SURFACE", true);
  surface->clear();
  IsoSurf::computeIsosurface(center, size, depth, scale, isCube, fGrid, *surface);

  // Recompute normals
  surface->setNormalPerVertex(false);
  _computeNormalPerFace(*surface);

  // we no longer need the point partition
  _deletePartition();

  std::cerr << "}" << endl;
}


//////////////////////////////////////////////////////////////////////
void SceneGraphProcessor::fitOptimalJacobi
(const Vec3f& center, const Vec3f& size,
  const int depth, const float scale, const bool cube,
  vector<float>& fGrid /* input & output */) {

  cerr << "SceneGraphProcessor::fitOptimalJacobi() {" << endl;

  if (size.x <= 0.0f || size.y <= 0.0f || size.z <= 0.0f || scale <= 0.0f) return;

  // TODO
  // More details in the file DG2025-FCEN-TP5-Notes.pdf

  // steps 1) to 7) same as in fitWatertight()

  // int n = nGridVertices;

  // 8) these parameters should be exposed in the user interface, but
  // we will keep them here for now
  //
  // int    nIter  = 20;
  // float  lambda = 0.10f;
  // float  mu     = 0.5f;

  // 9) allocate arrays for the function values, the displacements,
  // and the weights; we need a separate array for the function
  // values, because fGrid contains the input function values which
  // are needed to compute the displacements; when we finish we
  // replace the fGrid values by the f values and return
  //
  // vector<float> f;
  // vector<float> df;
  // vector<float> wf;
  // float fErr,wErr;

  // // 10) initialize the function values using the input values
  // f.insert(f.end(),fGrid.begin(),fGrid.end());

  // // 11) iterate nIter times
  // for(int iIter=0;iIter<nIter;iIter++) {
  //
  //   // 12) zero accumulators df and wf
  //
  //   // Data term 
  //
  //   // 13) for each grid vertex iV of an occupied cell accumulate
  //   // (fGrid[iV]-f[iV]) in the corresponding displacement, and
  //   // increment the corresponding weight
  //
  //   // Regularization term
  //
  //   // 14) accumulate Laplacian
  //   // for each grid edge (jV,kV) {
  //   //   - accumulate lambda*(f[kV]-f[jV]) in df[jV] and increment the
  //   //     corresponding weight
  //   //   - accumulate lambda*(f[jV]-f[kV]) in df[kV] and increment the
  //   //     corresponding weight
  //   // }
  //
  //   // 15) normalize the displacements
  //   for(iV=0;iV<n;iV++)
  //     df[iV] /= wf[iV];
  //
  //   // 16) update the function values
  //   for(iV=0;iV<n;iV++)
  //     f[iV] += mu*df[iV];
  //
  //   // 17) measure and report the error
  //   // ...
  //   cerr << "    err = " << fErr << endl;
  // }
  //
  // // 18) save result; fGrid becomes the output
  // fGrid.clear();
  // fGrid.insert(fGrid.end(),f.begin(),f.end());

  // // 19) update the output surface
  // computeIsosurface(center,size,depth,scale,cube,fGrid);

  // // 20) we no longer need the point partition
  // _deletePartition();

  cerr << "}" << endl;
}
