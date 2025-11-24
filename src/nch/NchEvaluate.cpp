//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-07 20:24:33 taubin>
//------------------------------------------------------------------------
//
// NchEvaluate.cpp
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
// recycle code from your implementation of
//
// float NchProcessor::nchEvaluate
// (const vector<float>& coord,
//  const vector<float>& normal,
//  const vector<float>& nchRhos,
//  const float x0, const float x1, const float x2);


#include <cmath>
#include <iostream>
#include "NchEvaluate.hpp"
#include "NchThread.hpp"
#include "core/IsoSurf.hpp"

// static
QProgressBar* NchEvaluate::_progressBar = nullptr;
// static
void NchEvaluate::setProgressBar(QProgressBar* progressBar) {
  _progressBar = progressBar;
}

//////////////////////////////////////////////////////////////////////
// constructor for regular grid
NchEvaluate::NchEvaluate
(const vector<float>& coord,
 const vector<float>& normal,
 vector<float>& rho,
 const Vec3f& center, const Vec3f& size,
 const int depth, const float scale, const bool cube,
 vector<float>& fGrid,
 IndexedFaceSet* surface,
 const bool multiThreaded,
 QObject *parent):
  QThread(parent),
  _coord(coord),
  _normal(normal),
  _rho(rho),
  // for regular grid
  _center(center),
  _size(size),
  _depth(depth),
  _scale(scale),
  _cube(cube),
  // for adaptive grid
  _hgpPtr(nullptr),
  //
  _fGrid(fGrid),
  _surface(surface),
  _multiThreaded(multiThreaded) {
}

//////////////////////////////////////////////////////////////////////
// constructor for adaptive rid
NchEvaluate::NchEvaluate
(const vector<float>& coord,
 const vector<float>& normal,
 vector<float>& rho,
 HexGridPartition& hgp,
 vector<float>& fGrid,
 IndexedFaceSet* surface,
 const bool multiThreaded,
 QObject *parent):
  QThread(parent),
  _coord(coord),
  _normal(normal),
  _rho(rho),
  // for regular grid
  _center(),
  _size(),
  _depth(0),
  _scale(1.0),
  _cube(true),
  // for adaptive grid
  _hgpPtr(&hgp),
  //
  _fGrid(fGrid),
  _surface(surface),
  _multiThreaded(multiThreaded) {
}

//////////////////////////////////////////////////////////////////////
void NchEvaluate::run() {

  _fGrid.clear();

  if(_hgpPtr==nullptr) { // regular grid

    _runRegular();

    IsoSurf::computeIsosurface
      (_center,_size,_depth,_scale,_cube,_fGrid,*_surface);

  } else { // adaptive grid

    _runAdaptive();
    
    IsoSurf::computeIsosurface
      (*_hgpPtr,_fGrid,*_surface);
  }
}

//////////////////////////////////////////////////////////////////////
void NchEvaluate::_runAdaptive() {

  int nPoints  = (int)(_coord.size()/3);
  int nNormals = (int)(_normal.size()/3);
  int nRhos    = (int)(_rho.size());
  if(nPoints<=0 || nPoints!=nNormals || nPoints!=nRhos) return;

  HexGridPartition& hgp = *_hgpPtr;
  Vec3f& min = hgp.getMin();
  Vec3f& max = hgp.getMax();
  int    N   = hgp.getResolution();
  int    nGridVertices = hgp.getNumberOfVertices();
  
  float neg_infty = -std::numeric_limits<float>::max();
  _fGrid.resize(nGridVertices,neg_infty);

  if(_progressBar!=(QProgressBar*)0) {
    connect(this,SIGNAL(progressReset()),
            _progressBar,SLOT(reset()));
    connect(this,SIGNAL(progressSetRange(int,int)),
            _progressBar,SLOT(setRange(int,int)));

    emit progressReset();
    emit progressSetRange(0,nPoints+1);
  }

  // create the vertex coord array
  vector<float> coordGrid;
  {
    float x,y,z;
    int h,h0,h1,h2,iCell,iC,iCx,iCy,iCz,iVx,iVy,iVz,iVertex,iVgrid;
    map<int,int> vMap;
    map<int,int> first =  hgp.getFirstMap();
    map<int,int>::iterator iMap;

    for(iMap=first.begin();iMap!=first.end();iMap++) {
      iCell = iMap->first;
      iC = iCell; iCx = iC%N; iC/=N; iCy = iC%N; iC/=N; iCz = iC;
      // vertices
      //
      // 0 ---- 1
      // |\     |\
      // | 2 ---- 3
      // 4 +--- 5 |
      //  \|     \|
      //   6 ---- 7
      //
      for(h=0;h<8;h++) {
        h0 = (h  )%2; iVx   = iCx+h0;
        h1 = (h/2)%2; iVy   = iCy+h1;
        h2 = (h/4)%2; iVz   = iCz+h2;      
        iVertex = iVx+(N+1)*(iVy+(N+1)*iVz);
        if(vMap.count(iVertex)==0) {
          iVgrid = static_cast<int>(coordGrid.size()/3);
          vMap[iVertex] = iVgrid;
          x = (((float)(N-iVx))*min.x+((float)iVx)*max.x)/((float)N);
          y = (((float)(N-iVy))*min.y+((float)iVy)*max.y)/((float)N);
          z = (((float)(N-iVz))*min.z+((float)iVz)*max.z)/((float)N);
          coordGrid.push_back(x);
          coordGrid.push_back(y);
          coordGrid.push_back(z);
        }
      }
    }
  }

  if(_multiThreaded) {

    // _progressValue = 0;  
    NchThread::resetProgressValue();

    // create threads
    int nThreads = QThread::idealThreadCount()-2;
    NchThread** thread = new NchThread*[nThreads];

    for(int th=0;th<nThreads;th++) {
      thread[th] =
        new NchThread
        (th,nThreads,_coord,_normal,_rho,coordGrid,_fGrid);
      
      if(_progressBar!=(QProgressBar*)0) {
        connect(thread[th],SIGNAL(progressSetValue(int)),
                _progressBar,SLOT(setValue(int)));
      }
    }

    // start threads
    for(int th = 0; th < nThreads; th++)
      thread[th]->start(QThread::LowPriority);

    // wait for all the threads to finish
    for(int th = 0; th < nThreads; th++) {
      thread[th]->wait();
    }

    // disconnect each thread from the _progressBar and delete it
    for(int th = 0; th < nThreads; th++) {
      if(_progressBar!=(QProgressBar*)0) {
        disconnect(thread[th],SIGNAL(progressSetValue(int)),
                   _progressBar,SLOT(setValue(int)));
      }
      delete thread[th];
    }

  } else /* if(singleThread) */ { /////////////////////////////////

    if(_progressBar!=(QProgressBar*)0) {
      connect(this,SIGNAL(progressSetValue(int)),
              _progressBar,SLOT(setValue(int)));
    }

    int iVgrid;
    float x,y,z,f_xyz;
    float p0i,p1i,p2i,n0i,n1i,n2i,rho_i,d0,d1,d2,f_i;
    float a, b, c;
    for(iVgrid=0;iVgrid<nGridVertices;iVgrid++) {
      x = coordGrid[3*iVgrid  ];
      y = coordGrid[3*iVgrid+1];
      z = coordGrid[3*iVgrid+2];
      
      f_xyz = -std::numeric_limits<float>::max();
      for(int i=0;i<nPoints;i++) {

          p0i = _coord[3*i  ];
          p1i = _coord[3*i+1];
          p2i = _coord[3*i+2];

          n0i = _normal[3*i  ];
          n1i = _normal[3*i+1];
          n2i = _normal[3*i+2];

          d0 = x - p0i;
          d1 = y - p1i;
          d2 = z - p2i;
          // Compute n_i^t (x - p_i)
          a = n0i*d0 + n1i*d1 + n2i*d2;

          // Compute ||x - p_i||^2
          b = d0*d0 + d1*d1 + d2*d2;

          // Compute a - rho_i b
          rho_i = _rho[i];
          c = a - rho_i * b;
          if (c > f_xyz) {
              f_xyz = c;
          }

      }
      _fGrid[iVgrid] = f_xyz;
    
      if(_progressBar!=(QProgressBar*)0) {
        emit progressSetValue(iVgrid);
      }
    }

    if(_progressBar!=(QProgressBar*)0) {
      disconnect(this,SIGNAL(progressReset()),
                 _progressBar,SLOT(reset()));
    }
  }
  
  if(_progressBar!=(QProgressBar*)0) {
    disconnect(this,SIGNAL(progressSetValue(int)),
            _progressBar,SLOT(setValue(int)));
    disconnect(this,SIGNAL(progressSetRange(int,int)),
            _progressBar,SLOT(setRange(int,int)));
  }
}

//////////////////////////////////////////////////////////////////////
void NchEvaluate::_runRegular() {

  if(_depth<0) return;
  if(_size.x<=0.0f || _size.y<=0.0f || _size.z<=0.0f || _scale<=0.0f) {
    return;
  }
  // int nGrid = 1<<_depth;
  
  int nPoints  = (int)(_coord.size()/3);
  int nNormals = (int)(_normal.size()/3);
  int nRhos    = (int)(_rho.size());
  if(nPoints<=0 || nPoints!=nNormals || nPoints!=nRhos) return;

  // compute the coordinates of the bounding box corners
  float dx=_size.x/2.0f, dy=_size.y/2.0f, dz=_size.z/2.0f;
  float dMax = dx; if(dy>dMax) dMax=dy; if(dz>dMax) dMax=dz;
  if(_cube      ) { dx = dy = dz = dMax; }
  if(_scale>0.0f) { dx *= _scale; dy *= _scale; dz *= _scale; }

  Vec3f min(_center.x-dx,_center.y-dy,_center.z-dz);
  Vec3f max(_center.x+dx,_center.y+dy,_center.z+dz);

  int N = 1<<_depth; // _nGrid;
  int nGridVertices = (N+1)*(N+1)*(N+1);
  float neg_infty = -std::numeric_limits<float>::max();
  _fGrid.resize(nGridVertices,neg_infty);

  if(_progressBar!=(QProgressBar*)0) {
    connect(this,SIGNAL(progressReset()),
            _progressBar,SLOT(reset()));
    connect(this,SIGNAL(progressSetRange(int,int)),
            _progressBar,SLOT(setRange(int,int)));

    emit progressReset();
    emit progressSetRange(0,nPoints+1);
  }

  if(_multiThreaded) {

    // _progressValue = 0;  
    NchThread::resetProgressValue();

    // create threads
    int nThreads = QThread::idealThreadCount()-2;
    NchThread** thread = new NchThread*[nThreads];

    for(int th=0;th<nThreads;th++) {
      thread[th] =
        new NchThread
        (th,nThreads,_coord,_normal,_rho,min,max,N,_fGrid);
      
      if(_progressBar!=(QProgressBar*)0) {
        connect(thread[th],SIGNAL(progressSetValue(int)),
                _progressBar,SLOT(setValue(int)));
      }
    }

    // start threads
    for(int th = 0; th < nThreads; th++)
      thread[th]->start(QThread::LowPriority);

    // wait for all the threads to finish
    for(int th = 0; th < nThreads; th++) {
      thread[th]->wait();
    }

    // disconnect each thread from the _progressBar and delete it
    for(int th = 0; th < nThreads; th++) {
      if(_progressBar!=(QProgressBar*)0) {
        disconnect(thread[th],SIGNAL(progressSetValue(int)),
                   _progressBar,SLOT(setValue(int)));
      }
      delete thread[th];
    }

  } else /* if(singleThread==false) */ {

    if(_progressBar!=(QProgressBar*)0) {
      connect(this,SIGNAL(progressSetValue(int)),
              _progressBar,SLOT(setValue(int)));
    }
  
    int iGridVertex,ix,iy,iz,i;
    float x,y,z;
    float p0i,p1i,p2i,n0i,n1i,n2i,rho_i,d0,d1,d2,f_i,f_xyz;
    float a, b, c;
    for(iGridVertex=iz=0;iz<=N;iz++) {
      z = (((float)(N-iz  ))*min.z+((float)(iz  ))*max.z)/((float)N);
      for(iy=0;iy<=N;iy++) {
        y = (((float)(N-iy  ))*min.y+((float)(iy  ))*max.y)/((float)N);      
        for(ix=0;ix<=N;ix++,iGridVertex++) {
          x = (((float)(N-ix  ))*min.x+((float)(ix  ))*max.x)/((float)N);

          f_xyz = -std::numeric_limits<float>::max();
          for(i=0;i<nPoints;i++) {

            p0i = _coord[3*i  ];
            p1i = _coord[3*i+1];
            p2i = _coord[3*i+2];

            n0i = _normal[3*i  ];
            n1i = _normal[3*i+1];
            n2i = _normal[3*i+2];
            
            d0 = x - p0i;
            d1 = y - p1i;
            d2 = z - p2i;
            // Compute n_i^t (x - p_i)
            a = n0i*d0 + n1i*d1 + n2i*d2;

            // Compute ||x - p_i||^2
            b = d0*d0 + d1*d1 + d2*d2;

            // Compute a - rho_i b
            rho_i = _rho[i];
            c = a - rho_i * b;
            if (c > f_xyz) {
              f_xyz = c;
            }

          }
          _fGrid[iGridVertex] = f_xyz;

          if(_progressBar!=(QProgressBar*)0) {
            emit progressSetValue(iGridVertex);
          }
        }
      }
    }

    if(_progressBar!=(QProgressBar*)0) {
      disconnect(this,SIGNAL(progressReset()),
                 _progressBar,SLOT(reset()));
    }
    
  }
  
  if(_progressBar!=(QProgressBar*)0) {
    disconnect(this,SIGNAL(progressSetValue(int)),
            _progressBar,SLOT(setValue(int)));
    disconnect(this,SIGNAL(progressSetRange(int,int)),
            _progressBar,SLOT(setRange(int,int)));
  }
}
