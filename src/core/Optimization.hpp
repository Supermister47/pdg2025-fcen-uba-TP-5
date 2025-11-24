//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-07 20:22:24 taubin>
//------------------------------------------------------------------------
//
// Optimization.hpp
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//     * Redistributions of source code must retain the above
//       copyright notice, this list of conditions and the following
//       disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials
//       provided with the distribution.
//     * Neither the name of the Brown University nor the names of its
//       contributors may be used to endorse or promote products
//       derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GABRIEL
// TAUBIN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
// OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#ifndef _OPTIMIZATION_HPP_
#define _OPTIMIZATION_HPP_

#include <QObject>
#include <wrl/IndexedFaceSet.hpp>
// #include "PolygonMesh.hpp"

class Optimization : public QObject {
    Q_OBJECT

public:

    Optimization();

    void clear();

    IndexedFaceSet* getInput();
    IndexedFaceSet* getOptimized();

    void  setInput(IndexedFaceSet* ifsInput);
    void  setOptimized(IndexedFaceSet* ifsOptimized, const bool reset=false);
    void  saveOptimized();

    int   getSteps();
    void  setSteps(const int value);

    float getLambda();
    void  setLambda(const float value);
    float getMu();
    void  setMu(const float value);
    float getKappa();
    void  setKappa(const float value);

    void  setSelectedVertexIndex(const int value);
    void  setSelectedEdgeIndex(const int value);
    void  setSelectedFaceIndex(const int value);

    // Jacobi iteration for energy function
    //
    // E(x) = \sum_{0<=iE<nE} {1} \| x_{iV0}-x_{iV1}^0\|^2
    //
    void  laplacianSmoothingVertexCoordinatesRun();

    // Jacobi iteration for energy function
    //
    // E(n) = \sum_{0<=iE<nED} {1} \| n_{iF0}-n_{iF1}^0\|^2
    //
    void  laplacianSmoothingFaceNormalsRun
        (const bool normalize=true);

    float getJacobiWeightData();
    void  setJacobiWeightData(const float value);
    float getJacobiWeightSmoothing();
    void  setJacobiWeightSmoothing(const float value);

    // Jacobi iteration for energy function
    //
    // E(x) = (1-sigma)*\sum_{0<=iV<nV} \| x_{iV}-x_{iV}^0\|^2 +
    //        (  sigma)*\sum_{0<=iE<nE} \| x_{iV0}-x_{iV1}^0\|^2
    //
    void  jacobiRun();

    vector<float>& getEdgeLengths();

    float getMinEdgeLength();
    float getMaxEdgeLength();
    float getTargetEdgeLength();

    void  setMinEdgeLength(const float value);
    void  setMaxEdgeLength(const float value);
    void  setTargetEdgeLength(const float value);

    void clusterVerticesApply
        (const Vec3f& min, const Vec3f& max, const int resolution);
    void clusterVerticesApply();
    int  getQuantizationResolution();
    void setQuantizationResolution
        (const int resolution);
    void setQuantizationBox
        (const Vec3f& min, const Vec3f& max);
    int getNumberOfClusters();

    enum EdgeCollapseIndependentSet {
        VERTICES_2, VERTICES_4, VERTICES_8
    };

    void  collapseEdgesShow
        (const EdgeCollapseIndependentSet indepSet);
    void  collapseEdgesApply
        (const EdgeCollapseIndependentSet indepSet);

    enum SplitEdgesMode {
        LONG, ALL, SELECTED
    };

    void  adaptiveSubdivisionShow
        (const SplitEdgesMode mode);
    void  adaptiveSubdivisionApply
        (const SplitEdgesMode mode, const bool colorFaces=false);

    void  equalizeValencesShow();
    void  equalizeValencesApply();

signals:

    // void refreshSignal(bool waiting);

    void progressReset();
    void progressSetRange(int min, int max);
    void progressSetValue(int value);
    void progressSetError(QString errorStr);

    // public slots:
    // void setStop(bool value);

private:

    void  _collapseEdgesSelect
        (const EdgeCollapseIndependentSet indepSet,
         vector<int>&                     edgeSelection,
         bool                             colorIncidentFaces);

    void  _adaptiveSubdivisionSelect
        (vector<int>& vertexSelection,
         vector<int>& edgeSelection,
         const SplitEdgesMode mode,
         bool colorIncidentFaces);

    void  _equalizeValencesSelect
        (vector<int>& edgeSelection, bool colorIncidentFaces);

private:

    // bool         _stop;

    IndexedFaceSet* _ifsInput;
    IndexedFaceSet* _ifsOptimized;
    int             _steps;
    float           _lambda;
    float           _mu;
    float           _jacobiWeightSmoothing;

    vector<float>   _edgeLengths;

    float           _minEdgeLength;
    float           _maxEdgeLength;
    float           _targetEdgeLength;

    Vec3f           _clusterMin;
    Vec3f           _clusterMax;
    int             _clusterResolution;
    int             _nClusters;

    int             _vSelIndex;
    int             _eSelIndex;
    int             _fSelIndex;
};

#endif /* _OPTIMIZATION_HPP_ */
