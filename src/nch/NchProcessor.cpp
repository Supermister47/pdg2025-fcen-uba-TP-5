//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-07 20:24:35 taubin>
//------------------------------------------------------------------------
//
// NchProcessor.cpp
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

// ASSIGNMENT 5 --- TODO
// in this file, only complete implementation of
//
// float NchProcessor::nchEvaluate
// (const vector<float>& coord,
//  const vector<float>& normal,
//  const vector<float>& nchRhos,
//  const float x0, const float x1, const float x2);

#include <cmath>
#include <iostream>
#include "NchProcessor.hpp"
#include "Shape.hpp"
#include "IndexedFaceSet.hpp"
#include "Appearance.hpp"
#include "Material.hpp"
#include "io/StrException.hpp"

// static
QProgressBar* NchProcessor::_progressBar = nullptr;
// static
void NchProcessor::setProgressBar(QProgressBar* progressBar) {
  _progressBar = progressBar;
  NchEstimate::setProgressBar(_progressBar);
}

NchProcessor::NchProcessor
(SceneGraph& wrl):
  _wrl(wrl),
  _source("POINTS"),
  _nchEstimateThread(nullptr),
  _nchEvaluateThread(nullptr) {
}

NchProcessor::~NchProcessor() {
}

void NchProcessor::setSource(QString& source) {
  _source = source;
}

IndexedFaceSet* NchProcessor::_getNamedShapeIFS
(const string& name, bool create) {
  IndexedFaceSet* ifs = (IndexedFaceSet*)0;
  Node* node = _wrl.getChild(name);
  if(node!=(Node*)0 && node->isShape()) {
    Shape* shape = (Shape*)node;
    node = shape->getGeometry();
    if(node!=(Node*)0 && node->isIndexedFaceSet()) {
      ifs = (IndexedFaceSet*)node;
    }
  }
  if(ifs==(IndexedFaceSet*)0 && create) {
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
// estimate the NCH radii for the SOURCE node
void NchProcessor::nchEstimateRhos(const bool multiThreaded) {
  
  try {

    // find the input point set in the scene graph
    // _source should be equal to "POINTS" or "SAMPLES"
    IndexedFaceSet* points  = _getNamedShapeIFS(_source.toStdString(),false);
    if(points==(IndexedFaceSet*)0)
      throw new StrException("points==nullptr");

    if(points->getNormalBinding()!=IndexedFaceSet::PB_PER_VERTEX)
      throw new StrException("point cloud does not have normal vectors");

    const vector<float>& coord  = points->getCoord();
    const vector<float>& normal = points->getNormal();

    Variable* varNchRhos = points->getVariable("nchRhos");
    if(varNchRhos==nullptr) {
      varNchRhos = new VariableVecFloat("nchRhos"); // initial size == 0
      points->setVariable(varNchRhos);
    }
    void* value = varNchRhos->getValue();
    vector<float>* nchRhosPtr =
      static_cast< vector<float>* >(value);
    if(nchRhosPtr==nullptr)
      throw new StrException("nchRhosPtr==nullptr");

    if(_nchEstimateThread!=(NchEstimate*)0) {
      if(_nchEstimateThread->isRunning()) return;
      delete _nchEstimateThread;
    }

    _nchEstimateThread =
      new NchEstimate(coord,normal,*nchRhosPtr,multiThreaded);
    
    _nchEstimateThread->run();

  } catch(StrException* e) {
    cout << "  ERROR : " << e->what() << endl;
    delete e;
  }
}

//////////////////////////////////////////////////////////////////////
// evaluate the NCH function at a 3D point (x0,x1,x2)
// rho values are already computed and passed as argument
// static
float NchProcessor::nchEvaluate
(const vector<float>& coord,
 const vector<float>& normal,
 const vector<float>& nchRhos,
 const float x0, const float x1, const float x2) {

  int nPoints  = (int)(coord.size()/3);  
  int nNormals = (int)(normal.size()/3);  
  int nRhos   = (int) nchRhos.size();
  if(nPoints<=0 || nPoints!=nNormals || nPoints!=nRhos) {
    return std::nanf(""); // ERROR
  }

  // f(x) = max_{i} f_i(x)
  // where
  // f_i(x) =n_i^t(x-p_i)-ncRhos[i]*\|x-p_i\|^2
  // n_i = (normal[3*i  ],normal[3*i+1],normal[3*i+2])
  // p_i = ( coord[3*i  ], coord[3*i+1], coord[3*i+2])

  // Hint: use the pseudo-code from the lecture slides

  // initialize fx to the most negative number
  float p0i,p1i,p2i,n0i,n1i,n2i,rho_i,d0,d1,d2,f_x;
  float a,b,c;

  float fx = -std::numeric_limits<float>::max();
  for(int iPoint=0;iPoint<nPoints;iPoint++) {

    // Get point p_i coordinates
    float px = coord[3*iPoint  ];
    float py = coord[3*iPoint+1];
    float pz = coord[3*iPoint+2];

    // Get normal n_i components
    float nx = normal[3*iPoint  ];
    float ny = normal[3*iPoint+1];
    float nz = normal[3*iPoint+2];

    // Get rho_i
    float rho = nchRhos[iPoint];

    // Vector d = x - p_i
    float dx = x0 - px;
    float dy = x1 - py;
    float dz = x2 - pz;

    // Compute n_i^t * (x - p_i)
    float a = nx*dx + ny*dy + nz*dz;

    // Compute ||x - p_i||^2
    float b = dx*dx + dy*dy + dz*dz;

    // f_i(x) = dot - rho * dist2
    float fi = a - rho * b;

    // Update maximum
    if(fi > fx) {
      fx = fi;
    }

  }
  
  return fx;
}

//////////////////////////////////////////////////////////////////////
void NchProcessor::nchReconstruct
(const bool multiThreaded,
 const Vec3f& center, const Vec3f& size,
 const int depth, const float scale, const bool cube,
 vector<float>& fGrid) {

  try {
    
    IndexedFaceSet* points = _getNamedShapeIFS(_source.toStdString(),false);
    if(points==(IndexedFaceSet*)0) return;
    if(points->getNormalBinding()!=IndexedFaceSet::PB_PER_VERTEX) return;
    vector<float>& coord  = points->getCoord();
    vector<float>& normal = points->getNormal();
    
    int nPoints  = (int)(coord.size()/3);  
    int nNormals = (int)(normal.size()/3);  
    if(nPoints<=0 || nPoints!=nNormals) return;
    
    Variable* varNchRhos = points->getVariable("nchRhos");
    if(varNchRhos==nullptr) return;
    void* value = varNchRhos->getValue();
    vector<float>* nchRhosPtr =
      static_cast< vector<float>* >(value);
    if(nchRhosPtr==nullptr) return;

    IndexedFaceSet* surface = _getNamedShapeIFS("SURFACE",true);
 
    if(_nchEvaluateThread!=(NchEvaluate*)0) {
      if(_nchEvaluateThread->isRunning()) return;
      delete _nchEvaluateThread;
    }
    
    _nchEvaluateThread =
      new NchEvaluate
      (coord,normal,*nchRhosPtr,
       center,size,depth,scale,cube,fGrid,surface,multiThreaded);
    
    _nchEvaluateThread->run();

  } catch(StrException* e) {
    cout << "  ERROR : " << e->what() << endl;
    delete e;
  }

}

void NchProcessor::nchReconstruct
(const bool multiThreaded, HexGridPartition& hgp, vector<float>& fGrid) {

  try {
    
    IndexedFaceSet* points = _getNamedShapeIFS(_source.toStdString(),false);
    if(points==(IndexedFaceSet*)0) return;
    if(points->getNormalBinding()!=IndexedFaceSet::PB_PER_VERTEX) return;
    vector<float>& coord  = points->getCoord();
    vector<float>& normal = points->getNormal();
    
    int nPoints  = (int)(coord.size()/3);  
    int nNormals = (int)(normal.size()/3);  
    if(nPoints<=0 || nPoints!=nNormals) return;
    
    Variable* varNchRhos = points->getVariable("nchRhos");
    if(varNchRhos==nullptr) return;
    void* value = varNchRhos->getValue();
    vector<float>* nchRhosPtr =
      static_cast< vector<float>* >(value);
    if(nchRhosPtr==nullptr) return;

    IndexedFaceSet* surface = _getNamedShapeIFS("SURFACE",true);
 
    if(_nchEvaluateThread!=(NchEvaluate*)0) {
      if(_nchEvaluateThread->isRunning()) return;
      delete _nchEvaluateThread;
    }
    
    _nchEvaluateThread =
      new NchEvaluate
      (coord,normal,*nchRhosPtr,hgp,fGrid,surface,multiThreaded);
    
    _nchEvaluateThread->run();

  } catch(StrException* e) {
    cout << "  ERROR : " << e->what() << endl;
    delete e;
  }
}
