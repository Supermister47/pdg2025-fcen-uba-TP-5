//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-07 20:24:37 taubin>
//------------------------------------------------------------------------
//
// NchThread.cpp
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

// ASSIGNMENT 5
//
//  search for TODO comments in this file
//
#include <iostream>
#include "NchThread.hpp"
#include "NchProcessor.hpp"
#include <QMutex>

int NchThread::_progressValue = 0;

NchThread::NchThread
(const int idThread,
 const int nThreads,
 const std::vector<float>& coord,
 const std::vector<float>& normal,
 std::vector<float>& rho):
  _id(idThread),
  _nThreads(nThreads),
  _coord(coord),
  _normal(normal),
  _rho(rho),
  _x0(0),_x1(0),_y0(0),
  _y1(0),_z0(0),_z1(0),
  _nGrid(0),_fGrid(nullptr),
  _mode(ESTIMATE) {
}

NchThread::NchThread
(const int idThread,
 const int nThreads,
 const std::vector<float>& coord,
 const std::vector<float>& normal,
 std::vector<float>& rho,
 const Vec3f& min, const Vec3f& max,
 const int nGrid,
 std::vector<float>& fGrid):
  _id(idThread),
  _nThreads(nThreads),
  _coord(coord),
  _normal(normal),
  _rho(rho),
  _x0(min.x),_x1(max.x),
  _y0(min.y),_y1(max.y),
  _z0(min.z),_z1(max.z),
  _nGrid(nGrid),
  _coordGridPtr(nullptr),
  _fGrid(&fGrid),
  _mode(EVALUATE) {
}

NchThread::NchThread
  (const int idThread,
   const int nThreads,
   const std::vector<float>& coord,
   const std::vector<float>& normal,
   std::vector<float>& rho,
   vector<float>& coordGrid,
   std::vector<float>& fGrid):
  _id(idThread),
  _nThreads(nThreads),
  _coord(coord),
  _normal(normal),
  _rho(rho),
  _x0(0),_x1(0),
  _y0(0),_y1(0),
  _z0(0),_z1(0),
  _nGrid(0),
  _coordGridPtr(&coordGrid),
  _fGrid(&fGrid),
  _mode(EVALUATE) {
}

void NchThread::_runEstimate() {
  int nPoints = static_cast<int>(_coord.size()/3);

  // Hint: use the pseudo-code from the lecture slides

  int   i,j;

  // the partition of work into threads is based on point index
float pix,piy,piz,pjx,pjy,pjz,nix,niy,niz,d0,d1,d2,a,b,rho_i;

    for(i=0;i<nPoints;i++) {

      // get the (x,y,z) coordinates to the i-th point from the _coord array
      pix = _coord[3*i+0];
      piy = _coord[3*i+1];
      piz = _coord[3*i+2];

      // Get the normal components
      nix = _normal[3*i+0];
      niy = _normal[3*i+1];
      niz = _normal[3*i+2];

      rho_i  = 0.0f;
      for(j=0;j<nPoints;j++) {
        if(j==i) continue;
        // get the (x,y,z) coordinates to the j-th point from the _coord array
        pjx = _coord[3*j+0];
        pjy = _coord[3*j+1];
        pjz = _coord[3*j+2];

        // Compute a = n_i^t = (pj - pi)
        a =  nix*(pjx - pix) +
             niy*(pjy - piy) +
             niz*(pjz - piz);

        // Compute b = ||pj - pi||^2
        d0 = pjx - pix;
        d1 = pjy - piy;
        d2 = pjz - piz;
        b = d0*d0 + d1*d1 + d2*d2;

        
        if ((a - rho_i * b) > 0.0f) {
            if (b > 0.0f) {
                rho_i = a / b;
            }
        }
        
      }
    _rho[i] = rho_i;

    QMutex mutex;
    mutex.lock();
    if(++_progressValue%1000==0) {
      emit progressSetValue(_progressValue);    
      // getApp()->processEvents();
    }
    mutex.unlock();
  }
}

void NchThread::_runEvaluateRegular() {

  // TODO
  //
  // you can call the static function
  //
  //   float NchProcessor::nchEvaluate
  //   (const vector<float>& coord,
  //    const vector<float>& normal,
  //    const vector<float>& nchRhos,
  //    const float x0, const float x1, const float x2);
  //
  // or only copy the inner loop

  int nPoints = static_cast<int>(_coord.size()/3);
  int N = _nGrid;

  float x,y,z,f_iV,p0i,p1i,p2i,n0i,n1i,n2i,d0,d1,d2,a,b,c;
  int i,iV,ix,iy,iz;

  // the partition of work into threads is based on the z layer index
  for(iz=_id;iz<=N;iz+=_nThreads) {
    // z0<=z<=z1
    z = (((float)(N-iz  ))*_z0+((float)(iz  ))*_z1)/((float)N);
  
    for(iy=0;iy<=N;iy++) {
      // y0<=y<=y1
      y = (((float)(N-iy  ))*_y0+((float)(iy  ))*_y1)/((float)N);

      for(ix=0;ix<=N;ix++/*,iV++*/) {
        // x0<=x<=x1
        x = (((float)(N-ix  ))*_x0+((float)(ix  ))*_x1)/((float)N);

        // compute the scan order grid vertex index
        iV   = ix+(N+1)*(iy+(N+1)*iz);

        // f_iV = max_{i} f_i(x,y,z)
        // where
        // f_i(x) =n_i^t((x,y,z)^t-p_i)-_rho[i]*\|(x,y,z)^t-p_i\|^2
        
        f_iV = NchProcessor::nchEvaluate(_coord, _normal, _rho, x, y, z);

        (*_fGrid)[iV] = f_iV;

        QMutex mutex;
        mutex.lock();
        if(++_progressValue%1000==0) {
          emit progressSetValue(_progressValue);    
          // getApp()->processEvents();
        }
        mutex.unlock();
      }
    }
  }
}

void NchThread::_runEvaluateAdaptive() {

  //
  // you can call the static function
  //
  //   float NchProcessor::nchEvaluate
  //   (const vector<float>& coord,
  //    const vector<float>& normal,
  //    const vector<float>& nchRhos,
  //    const float x0, const float x1, const float x2);
  //
  // or only copy the inner loop

  int nPoints = static_cast<int>(_coord.size()/3);
  int nGridVertices = static_cast<int>(_coordGridPtr->size()/3);

  float x,y,z,f_iV,p0i,p1i,p2i,n0i,n1i,n2i,d0,d1,d2,a,b,c;
  int i=0,iV=0,iVgrid=0;
  
  // the partition of work into threads is based on the iVgrid index
  for(iVgrid=_id;iVgrid<nGridVertices;iVgrid+=_nThreads) {
    
    // get the (x,y,z) coordinates of the iVgrid vertex from the
    // (*_coordGridPtr) array
    x = (*_coordGridPtr)[3*iVgrid+0];
    y = (*_coordGridPtr)[3*iVgrid+1];
    z = (*_coordGridPtr)[3*iVgrid+2];

    // inner loop same as in _runEvaluateRegular()

    // f_iV = max_{i} f_i(x,y,z)
    // where
    // f_i(x) =n_i^t((x,y,z)^t-p_i)-_rho[i]*\|(x,y,z)^t-p_i\|^2
        
    f_iV = NchProcessor::nchEvaluate(_coord, _normal, _rho, x, y, z);

    (*_fGrid)[iVgrid] = f_iV;

    QMutex mutex;
    mutex.lock();
    if(++_progressValue%1000==0) {
      emit progressSetValue(_progressValue);
      // getApp()->processEvents();
    }
    mutex.unlock();
  }
}

void NchThread::run() {
  if(_mode==ESTIMATE)
    _runEstimate();
  else if(_mode==EVALUATE) {
    if(_coordGridPtr==nullptr)
      _runEvaluateRegular();
    else
      _runEvaluateAdaptive();
  }
}
