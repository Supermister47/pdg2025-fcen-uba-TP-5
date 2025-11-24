//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-07 20:24:31 taubin>
//------------------------------------------------------------------------
//
// NchEstimate.cpp
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
// - search for TODO comments

#include <cmath>
#include <iostream>
#include "NchEstimate.hpp"
#include "NchThread.hpp"

// static
QProgressBar* NchEstimate::_progressBar = nullptr;
// static
void NchEstimate::setProgressBar(QProgressBar* progressBar) {
  _progressBar = progressBar;
}

//////////////////////////////////////////////////////////////////////
NchEstimate::NchEstimate
(const vector<float>& coord,
 const vector<float>& normal,
 vector<float>& rho,
 const bool multiThreaded,
 QObject *parent):
  QThread(parent),
  _coord(coord),_normal(normal),_rho(rho),_multiThreaded(multiThreaded) {
}

//////////////////////////////////////////////////////////////////////
// NchEstimate::~NchEstimate() {
//
// }

//////////////////////////////////////////////////////////////////////
void NchEstimate::run() {
  
  int nPoints  = (int)(_coord.size()/3);
  int nNormals = (int)(_normal.size()/3);
  if(nPoints<=0 || nPoints!=nNormals) return;

  _rho.clear();
  _rho.resize(nPoints,0.0f);

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
        new NchThread(th,nThreads,_coord,_normal,_rho);
      
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

    // Hint: use the pseudo-code from the lecture slides
    
    int i,j;
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
            rho_i = a / b;
        }
        
      }
      _rho[i] = rho_i;

      if(_progressBar!=(QProgressBar*)0) {
        emit progressSetValue(i);
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
