//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-10-22 15:23:34 taubin>
//------------------------------------------------------------------------
//
// SceneGraphProcessor.hpp
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

#ifndef SCENE_GRAPHPROCESSOR_HPP
#define SCENE_GRAPHPROCESSOR_HPP

#include <iostream>
#include <QObject>
#include "SceneGraph.hpp"
#include "Shape.hpp"
#include "IndexedFaceSet.hpp"
#include "IndexedLineSet.hpp"
#include "core/HexGridPartition.hpp"

class SceneGraphProcessor : public QObject {
Q_OBJECT

public:

  SceneGraphProcessor(SceneGraph& wrl);
  ~SceneGraphProcessor();

  int numberOfShapeNodes();

  void normalClear();
  void normalInvert();
  void computeNormalPerFace();
  void computeNormalPerVertex();
  void computeNormalPerCorner();

  void gridAdd(int depth=0, float scale=1.0f, bool isCube=true);
  void gridAdd(HexGridPartition& hgp);
  void gridRemove();
  bool hasGrid();

  void edgesAdd();
  void edgesRemove();
  bool hasEdges();

  bool hasIndexedFaceSetFaces();
  bool hasIndexedFaceSetNormalNone();
  bool hasIndexedFaceSetNormalPerFace();
  bool hasIndexedFaceSetNormalPerVertex();
  bool hasIndexedFaceSetNormalPerCorner();

  bool hasIndexedLineSetColorNone();
  bool hasIndexedLineSetColorPerVertex();
  bool hasIndexedLineSetColorPerPolyline();

  bool hasIndexedFaceSetShown();
  bool hasIndexedFaceSetHidden();
  bool hasIndexedLineSetShown();
  bool hasIndexedLineSetHidden();

  void shapeIndexedFaceSetShow();
  void shapeIndexedFaceSetHide();
  void shapeIndexedLineSetShow();
  void shapeIndexedLineSetHide();

  void removeSceneGraphChild(const string& name);
  void pointsRemove();
  void samplesRemove();
  void planeRemove();
  void surfaceRemove();

  void fitSinglePlane
  (const Vec3f& center, const Vec3f& size, const float scale, const bool cube,
   Vec4f& f /* temporary output */ );


  void fitMultiplePlanes
  (const Vec3f& center, const Vec3f& size,
   const int depth, const float scale, const bool cube,
   vector<float>& fVec);

  void fitContinuous
  (const Vec3f& center, const Vec3f& size,
   const int depth, const float scale, const bool cube,
   vector<float>& fVec );

  void fitWatertight
  (const Vec3f& center, const Vec3f& size,
   const int depth, const float scale, const bool cube,
   vector<float>& fGrid);

  void fitOptimalJacobi
  (const Vec3f& center, const Vec3f& size,
   const int depth, const float scale, const bool cube,
   vector<float>& fGrid /* input & output */);

signals:

  void progressReset();
  void progressSetRange(int min, int max);
  void progressSetValue(int value);

private:

  // edges of hexahedron
  //
  // 0 ---- 1
  // |\     |\
  // | 2 ---- 3
  // 4 +--- 5 |
  //  \|     \|
  //   6 ---- 7
  //
  static const int _hexGridEdge[12][2];

  SceneGraph&    _wrl;

  // computation grid and point set partition
  int         _nGrid;
  int         _nPoints;
  int*        _next;
  int*        _first;

  int         _nCells();
  int         _createPartition
              (Vec3f& min, Vec3f& max, int depth, vector<float>&coord);
  void        _deletePartition();
  void        _applyToIndexedFaceSet(IndexedFaceSet::Operator p);

  // IndexedFaceSet::Operator
  static void _normalClear(IndexedFaceSet& ifs);
  static void _normalInvert(IndexedFaceSet& ifs);
  static void _computeNormalPerFace(IndexedFaceSet& ifs);
  static void _computeNormalPerVertex(IndexedFaceSet& ifs);
  static void _computeNormalPerCorner(IndexedFaceSet& ifs);

  static void _computeFaceNormal
              (vector<float>& coord, vector<int>&   coordIndex,
               int i0, int i1, Vec3f& n, bool normalize);

  bool        _hasShapeProperty(Shape::Property p);
  bool        _hasIndexedFaceSetProperty(IndexedFaceSet::Property p);
  bool        _hasIndexedLineSetProperty(IndexedLineSet::Property p);

  // Shape::Property
  static bool _hasEdges(Shape& shape);
  static bool _hasIndexedFaceSetShown(Shape& shape);
  static bool _hasIndexedFaceSetHidden(Shape& shape);
  static bool _hasIndexedLineSetShown(Shape& shape);
  static bool _hasIndexedLineSetHidden(Shape& shape);

  // IndexedFaceSet::Property
  static bool _hasFaces(IndexedFaceSet& ifs);
  static bool _hasNormalNone(IndexedFaceSet& ifs);
  static bool _hasNormalPerFace(IndexedFaceSet& ifs);
  static bool _hasNormalPerVertex(IndexedFaceSet& ifs);
  static bool _hasNormalPerCorner(IndexedFaceSet& ifs);
    
  // IndexedLineSet::Property
  static bool _hasColorNone(IndexedLineSet& ils);
  static bool _hasColorPerVertex(IndexedLineSet& ils);
  static bool _hasColorPerPolyline(IndexedLineSet& ils);

  IndexedFaceSet* _getNamedShapeIFS(const string& name, bool create);

  void _fitPlane (const Vec3f& center, const Vec3f& size,
    const float scale, const bool cube, Vec4f& f,
    const vector<float>& coordPoints, const vector<float>& normalPoints,
    float* FGrid=nullptr, bool* bGrid=nullptr,
    bool showPlane=true);

};

#endif // SCENE_GRAPHPROCESSOR_HPP
