//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-07 20:22:17 taubin>
//------------------------------------------------------------------------
//
// HexGridPartition.hpp
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
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

// This class is used to partition a set of 3D points into buckets
// corresponding to the cells of a regular hexahedral grid.
// The coordinates of the points are quantized to the desired
// resolution with respect to a bounding box, and inserted into
// the corresponding buckets. The buckets are represented as lists
// using an internal vector<int> next array. The first element of
// the bucket lists are accessed throup a map<int,int>. where the
// first index of the bucket (ix,iy,iz) is the scan order index
// ix+N*(iy+N*iz), where N is the number of cells each axis of
// the bounding box is subdivided into. The second index is the
// point index of the first point in the bucket list. Subsequent
// points are found using the next array. The last point in each
// bucket list has a -1 as its next value.
//
// We will use it to implement a simple version of the vertex clustering
// simplification method, which was pending since we worked on loading
// STL files. The Optimization panel in the GUI has a new section named
// "CLUSTER VERTICES", and the class core/Optimization has some new
// methods to support this function, which you also have to
// implement. Details below.
//
// We will also use the core/HexGridPartition to subsample a point cloud
// by selecting one sample point from each occupied cell. Since choosing
// the original point closest to the mean of the points contained in each
// cell has been shown to produce good results in practice, you will have
// to implement this sample() method of the class. There is a new
// Implicit panel in the GUI where this sampling method is used to
// generate a new node named SAMPLE in the scene graph. The Implicit
// panel allows you to chose to run NCH Surface reconstruction method
// either on the full point cloud, or on the subsampled point cloud to
// compare speeds. Again, more details below.
//
// Finally, we will also use the occupied cells of the
// core/HexGridPartition class to implicitly generate sparse regular
// grids in an attempt to save time in the computation of isosurfaces for
// surface reconstruction. More details below.

#include <cmath>
#include <iostream>
#include "HexGridPartition.hpp"

HexGridPartition::HexGridPartition
(const Vec3f& min, const Vec3f& max, const int resolution):
  _min(min),
  _max(max),
  _resolution(resolution),
  _nPointsInsideBox(0),
  _nPointsOutsideBox(0),
  _next(),
  _first(),
  _coord(nullptr),
  _normal(nullptr) {

}

HexGridPartition::HexGridPartition
(const Vec3f& center, const Vec3f& size,
 const int resolution, const float scale, const bool cube):
  _min(),
  _max(),
  _resolution(resolution),
  _nPointsInsideBox(0),
  _nPointsOutsideBox(0),
  _next(),
  _first(),
  _coord(nullptr),
  _normal(nullptr) {
  float dx=size.x/2.0f, dy=size.y/2.0f, dz=size.z/2.0f;
  float dMax = dx; if(dy>dMax) dMax=dy; if(dz>dMax) dMax=dz;
  if(cube      ) { dx = dy = dz = dMax; }
  if(scale>0.0f) { dx *= scale; dy *= scale; dz *= scale; }
  _min.x = center.x-dx; _min.y = center.y-dy; _min.z = center.z-dz;
  _max.x = center.x+dx; _max.y = center.y+dy; _max.z = center.z+dz;
}

bool HexGridPartition::insertPoints(const vector<float>& coord) {
  _nPointsInsideBox = 0;
  _nPointsOutsideBox = 0;
  // reset the lists
  _next.clear();
  _first.clear();

  // first some error checking ...
  if(_resolution<1) return false; // failure
  int nPoints = static_cast<int>(coord.size()/3);
  // _min & _max define the bounding box
  float dx = _max.x-_min.x;
  float dy = _max.y-_min.y;
  float dz = _max.z-_min.z;
  if(nPoints<=0 || dx<=0.0f || dy<=0.0f || dz<=0.0f) {
    _coord = nullptr;
    _normal = nullptr;
    return false; // failure
  }

  _coord  = &coord;
  // initialize the lists
  _next.resize(nPoints,-1);

  int N  = _resolution;
  int n0 = N;
  int n1 = N;
  int n2 = N;
  int iCell,iPoint,ix,iy,iz;
  float x,y,z;
  for(iPoint=0;iPoint<nPoints;iPoint++) {
    // 1) get the x,y,z coordinates of the iPoint
    x = coord[3*iPoint+0];
    y = coord[3*iPoint+1];
    z = coord[3*iPoint+2];
    // 2) if the points is not contained in the bounding box
    //    - increase the _nPointsOutsideBox variable
    //    - do not insert the point in any list
    if (x<_min.x || x>=_max.x ||
        y<_min.y || y>=_max.y ||
        z<_min.z || z>=_max.z) {
      _nPointsOutsideBox++;
      continue;
    }
    // 3) quantize the x,y,z coordinates to ix,iy,iz
    //    - so that [_min.x:_max.x) --> 0<=ix<N
    //    - same for iy, and iz
    ix = static_cast<int>(floorf( (x - _min.x) / dx * N ));
    iy = static_cast<int>(floorf( (y - _min.y) / dy * N ));
    iz = static_cast<int>(floorf( (z - _min.z) / dz * N ));
    if (ix<0) ix=0; else if (ix>=N) ix=N-1;
    if (iy<0) iy=0; else if (iy>=N) iy=N-1;
    if (iz<0) iz=0; else if (iz>=N) iz=N-1;

    // 4) increment _nPointsInsideBox
    _nPointsInsideBox++;
    // 5) compute the scan order index iCell cooresponding to
    //    the cell index (ix,iy,iz)
    int i = ix + N*(iy + N*iz);
    
    // 6) insert iPoint in the iCell list
    //    - Note: use _first.count(iCell) to find out whether or not
    //    - the iCell list is empty or not; otherwise you will end up
    //    - with _first being a dense map 
    if(_first.count(i) == 0) {
      _first[i] = iPoint;
      _next[iPoint] = -1;
    }
    else {
      int iNext = _first[i];
      // traverse the list to find the last point
      while (_next[iNext] != -1) {
        iNext = _next[iNext];
      }

      _next[iNext] = iPoint;
      _next[iPoint] = -1;
    }
  }

  return true; // success
}

bool HexGridPartition::sample
(vector<float>& coordSample, vector<int>* vMap) {

  vector<float> normalSample;
  return sample(coordSample,normalSample,vMap);
}

bool HexGridPartition::insertPoints
(const vector<float>& coord, const vector<float>& normal) {
  bool success = false;
  int nPoints = static_cast<int>(coord.size()/3);
  int nNormal = static_cast<int>(normal.size()/3);
  if(nPoints>0 && nNormal==nPoints) {
    success = insertPoints(coord);
  }
  if(success)  {
    _coord  = &coord;
    _normal = &normal;
  } else {
    _coord = nullptr;
    _normal = nullptr;
  }
  return success;
}

bool HexGridPartition::sample
(vector<float>& coordSample, vector<float>& normalSample, vector<int>* vMap) {

  // same as
  // sample(vector<float>& coordSample, vector<int>* vMap)
  //
  // but also clearing normalSample at the begining
  //
  // and pushing back the normal coordinates of iPointMin onto the
  // normalSample array

    if(_coord==nullptr) return false;

  coordSample.clear();
  // if vMap!=nullpointer, on output *vMap maps input point indices
  // onto output sample indices
  if(vMap!=nullptr) {
    // resize and initialize the output vMap
    vMap->clear();
    int nPoints = static_cast<int>(_coord->size()/3);
    vMap->resize(nPoints,-1);
  }

  normalSample.clear();
  
  // int N=_resolution,n0=N,n1=N,n2=N,ix,iy,iz,iCell;
  int k,iPointFirst,iPoint,nPointsCell,iPointMin,iSamplePoint;
  float xMean,yMean,zMean,fPointsInCell,dx,dy,dz,dPoint2,dMin2;

  // traverse the cell of the _first map
  map<int,int>::iterator iMap = _first.begin();
  for(;iMap!=_first.end();iMap++) {
    int iCell = iMap->first;

    // if the points were properly inserted this value should be >=0 
    iPointFirst = iMap->second;
    //assert(iPointFirst>=0);
    // get the next available coordSample index
    iSamplePoint = static_cast<int>(coordSample.size()/3);

    // count number of points in cell
    // and mean coords of points contained in the cell
    nPointsCell = 0;
    xMean=yMean=zMean=0.0f;
    for(iPoint=iPointFirst;iPoint>=0;iPoint=_next[iPoint]) {
      // 1) get x,y,z coords of iPoint
      float x = _coord->at(3*iPoint+0);
      float y = _coord->at(3*iPoint+1);
      float z = _coord->at(3*iPoint+2);
      // 2) accumulate xMean,yMean,zMean
      xMean += x;
      yMean += y;
      zMean += z;
      // 3) increment nPointsCell
      nPointsCell++;
      // 4) if(vMap!=nullptr)
      //    - map iPoint onto iSamplePoint
      if (vMap != nullptr) {
        vMap->at(iPoint) = iSamplePoint;
      }
    }
    // normalize xMean,yMean,zMean
    xMean /= static_cast<float>(nPointsCell);
    yMean /= static_cast<float>(nPointsCell);
    zMean /= static_cast<float>(nPointsCell);

    // find amongst the points contained in the cell, the closest to the mean
    dMin2 = std::numeric_limits<float>::max();
    iPointMin = -1;
    for(iPoint=iPointFirst;iPoint>=0;iPoint=_next[iPoint]) {
      float x = _coord->at(3*iPoint+0);
      float y = _coord->at(3*iPoint+1);
      float z = _coord->at(3*iPoint+2);

      // compute squared distance to mean
      float dx = xMean - x;
      float dy = yMean - y;
      float dz = zMean - z;
      float d1 = dx*dx + dy*dy + dz*dz;
      float d2 = d1*d1;
      if (d2 < dMin2) {
        dMin2 = d2;
        iPointMin = iPoint;
      }
    }
    coordSample.push_back(_coord->at(3*iPointMin+0));
    coordSample.push_back(_coord->at(3*iPointMin+1));
    coordSample.push_back(_coord->at(3*iPointMin+2));

    // if normals where provided, push back the normal of iPointMin
    if (_normal != nullptr) {
        normalSample.push_back(_normal->at(3*iPointMin+0));
        normalSample.push_back(_normal->at(3*iPointMin+1));
        normalSample.push_back(_normal->at(3*iPointMin+2));
    }

  }
  return true;
}

int HexGridPartition::getNumberOfPoints() {
  if(_coord==nullptr || _normal==nullptr || _coord->size()!=_normal->size())
    return 0;
  else
    return static_cast<int>(_coord->size()/3);
}

int HexGridPartition::getFirst(const int ix, const int iy, const int iz) {
  int N  = _resolution, n0 = N, n1 = N /*, n2 = N*/;
  if(ix<0 || ix>=n0 || iy<0 || iy>=n0 || iz<0 || iz>=n0) return -1;
  int iCell = ix+n0*(iy+n1*iz); // 0<=iCell<(N+1)*(N+1)*(N+1)

  // if _first does not have an entry for iCell, and we call
  // _first[iCell] directly here, an entry will be inserted using the
  // default value; this can be prevented using _first.count(iCell) instead
  
  if(_first.count(iCell)==0) return -1;
  return _first[iCell];
}

int HexGridPartition::getNext(const int iPoint) {
  int nPoints = getNumberOfPoints();
  return (iPoint<0 || iPoint>=nPoints)?-1:_next[iPoint];
}


int HexGridPartition::getNumberOfVertices() {
  int nV = 0;
  int h,h0,h1,h2,iCell,iC,iCx,iCy,iCz,iVx,iVy,iVz,iVertex;
  int N = _resolution;

  // - this map is used to count each cell vertex only once
  //
  // - it also enumerates the regular vertex indices which are
  // - vertices of at least one cell, but that information is not used
  //   here, though
  // - check how it is used in
  //   void IsoSurf::computeIsosurface
  //     (HexGridPartition& hgp, vector<float>& fGrid, IndexedFaceSet& isoSurface) 
  //
  map<int,int> vMap;

  map<int,int>::iterator iMap;
  for(iMap=_first.begin();iMap!=_first.end();iMap++) {
    iCell       = iMap->first;
    iC = iCell; iCx = iC%N; iC/=N; iCy = iC%N; iC/=N; iCz = iC;
    //
    // 0 ---- 1
    // |\     |\
    // | 2 ---- 3
    // 4 +--- 5 |
    //  \|     \|
    //   6 ---- 7
    //
    for(h=0;h<8;h++) {
      h0 = (h  )%2; iVx = iCx+h0;
      h1 = (h/2)%2; iVy = iCy+h1;
      h2 = (h/4)%2; iVz = iCz+h2;
      iVertex = iVx+(N+1)*(iVy+(N+1)*iVz);
      // count this vertex only on the first visit
      if(vMap.count(iVertex)==0) {
        vMap[iVertex] = nV++;
      }
    }
  }
  
  return nV;
}
