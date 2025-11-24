//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-07 20:30:47 taubin>
//------------------------------------------------------------------------
//
// SaverStl.cpp
//
// Written by: <Your Name>
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

#include <cstring>

#include "SaverStl.hpp"
#include "StrException.hpp"
#include <stdint.h>

#include "wrl/Shape.hpp"
// #include "wrl/Appearance.hpp"
// #include "wrl/Material.hpp"
// #include "core/Faces.hpp"

const char* SaverStl::_ext = "stl";
SaverStl::FileType SaverStl::_fileType = SaverStl::FileType::ASCII;

//////////////////////////////////////////////////////////////////////
// static
void SaverStl::setFileType(const SaverStl::FileType ft) {
  _fileType = ft;
}
  
//////////////////////////////////////////////////////////////////////
bool SaverStl::_saveAscii
(FILE* fp, const char* solidname, IndexedFaceSet& ifs) const {

  int nF = ifs.getNumberOfFaces();
  vector<float>& coord       = ifs.getCoord();
  vector<int>&   coordIndex  = ifs.getCoordIndex();
  vector<float>& normal      = ifs.getNormal();
  vector<int>&   normalIndex = ifs.getNormalIndex();
  // already checked that ifs.getNormalPerVertex()==false
  bool           npf_indexed = (static_cast<int>(normalIndex.size())==nF);

  fprintf(fp,"solid %s\n",solidname);
    
  int iF,iV0,iV1,iV2,iN;
  float x0,x1,x2,n0,n1,n2;
  for(iF=0;iF<nF;iF++) { // for each face ...
    // use fprintf() to print formatted text

    iN = (npf_indexed)?normalIndex[iF]:iF;
    n0 = normal[3*iN  ];
    n1 = normal[3*iN+1];
    n2 = normal[3*iN+2];

    fprintf(fp,"facet normal %e %e %e\n",n0,n1,n2);
    fprintf(fp,"  outer loop\n");

    iV0 = coordIndex[4*iF+0];
    x0  = coord[3*iV0  ];
    x1  = coord[3*iV0+1];
    x2  = coord[3*iV0+2];

    fprintf(fp,"    vertex %e %e %e\n",x0,x1,x2);

    iV1 = coordIndex[4*iF+1];
    x0  = coord[3*iV1  ];
    x1  = coord[3*iV1+1];
    x2  = coord[3*iV1+2];

    fprintf(fp,"    vertex %e %e %e\n",x0,x1,x2);

    iV2 = coordIndex[4*iF+2];
    x0  = coord[3*iV2  ];
    x1  = coord[3*iV2+1];
    x2  = coord[3*iV2+2];

    fprintf(fp,"    vertex %e %e %e\n",x0,x1,x2);

    fprintf(fp,"  endloop\n");
    fprintf(fp,"endfacet\n");

  }

  return true;
}

//////////////////////////////////////////////////////////////////////
bool SaverStl::_saveBinary
(FILE* fp, const char* solidname, IndexedFaceSet& ifs) const {

  int nF = ifs.getNumberOfFaces();
  vector<float>& coord       = ifs.getCoord();
  vector<int>&   coordIndex  = ifs.getCoordIndex();
  vector<float>& normal      = ifs.getNormal();
  vector<int>&   normalIndex = ifs.getNormalIndex();
  // already checked that ifs.getNormalPerVertex()==false
  bool           npf_indexed = (static_cast<int>(normalIndex.size())==nF);

  size_t written = 0;

  // allocate header and initialize to zero
  char header[80];
  memset(header,0x00,80);
  snprintf(header,80,"BINARY STL %s Exported by DGP2025",solidname);

  written = fwrite(header,1,80,fp);
  if(written!=80)
    throw new StrException("unable to write binary STL header");

  uint32_t nTriangles = static_cast<uint32_t>(nF);
  written = fwrite((void*)&nTriangles,1,4,fp);
  if(written<4)
    throw new StrException("unable to write number of triangles");

  uint16_t abc = 0x0000; // attribute byte count

  int iF,iV0,iV1,iV2,iN;
  float n[3],v[3];
  for(iF=0;iF<nF;iF++) {

    iN   = (npf_indexed)?normalIndex[iF]:iF;
    n[0] = normal[3*iN  ];
    n[1] = normal[3*iN+1];
    n[2] = normal[3*iN+2];
    written = fwrite(n,1,12,fp);
    if(written!=12)
      throw new StrException("unable to write normal vector");

    iV0  = coordIndex[4*iF+0];
    v[0] = coord[3*iV0  ];
    v[1] = coord[3*iV0+1];
    v[2] = coord[3*iV0+2];
    written = fwrite(v,1,12,fp);
    if(written!=12)
      throw new StrException("unable to write vertex 0");

    iV1  = coordIndex[4*iF+1];
    v[0] = coord[3*iV1  ];
    v[1] = coord[3*iV1+1];
    v[2] = coord[3*iV1+2];
    written = fwrite(v,1,12,fp);
    if(written!=12)
      throw new StrException("unable to write vertex 1");

    iV2  = coordIndex[4*iF+2];
    v[0] = coord[3*iV2  ];
    v[1] = coord[3*iV2+1];
    v[2] = coord[3*iV2+2];
    written = fwrite(v,1,12,fp);
    if(written!=12)
      throw new StrException("unable to write vertex 2");

    written = fwrite(&abc,1,2,fp);
    if(written<2)
      throw new StrException("unable to write atribute byte count");

  }

  return true;
}

//////////////////////////////////////////////////////////////////////
bool SaverStl::save(const char* filename, SceneGraph& wrl) const {
  bool success = false;
  FILE* fp = (FILE*)0;
  try {
    // Check these conditions
    if(filename==(char*)0)
      throw new StrException("empty filename");
    // 1) the SceneGraph should have a single child
    if(wrl.getNumberOfChildren()!=1)
      throw new StrException("number of SceneGraph children != 1");
    // 2) the child should be a Shape node
    Node* child_0 = wrl[0];
    Shape* shape = dynamic_cast<Shape*>(child_0);
    if(shape==(Shape*)0)
      throw new StrException("first SceneGraph child not a Shape node");
    // 3) the geometry of the Shape node should be an IndexedFaceSet node
    Node* geometry = shape->getGeometry();
    IndexedFaceSet* ifs = dynamic_cast<IndexedFaceSet*>(geometry);
    if(ifs==(IndexedFaceSet*)0)
      throw new StrException("Shape geometry not an IndexedFaceSet");
    // - construct an instance of the Faces class from the IndexedFaceSet
    // int nV = ifs->getNumberOfCoord();
    // vector<float>& coord      = ifs->getCoord();
    vector<int>&   coordIndex = ifs->getCoordIndex();

    // 4) the IndexedFaceSet should be a triangle mesh
    // - use the Faces class, or directly the coordIndex array to
    //   verify that all the faces are triangles
    // - if you find a face with more than thre vertices
    //   throw an exception

    // Faces faces(nV,coordIndex);
    // int nF = faces.getNumberOfFaces();
    // if(nF<1)
    //   throw new StrException("has no faces");
    // for(int iF=0;iF<nF;iF++) {
    //   if(faces.getFaceSize(iF)!=3)
    //     throw new StrException("is not a triangle mesh");
    // }

    int i0,i1,nFs;
    for(i0=i1=0;i0<static_cast<int>(coordIndex.size());i1++) {
      if(coordIndex[i1]>=0) continue;
      nFs = i1-i0; // size of face
      if(nFs!=3)
        throw new StrException("is not a triangle mesh");
      i0=i1+1;
    }

    // 5) verify that the IndexedFaceSet has normals per face
    IndexedFaceSet::Binding nb = ifs->getNormalBinding();
    bool npf_non_indexed = (nb==IndexedFaceSet::Binding::PB_PER_FACE);
    bool npf_indexed     = (nb==IndexedFaceSet::Binding::PB_PER_FACE_INDEXED);
    if(npf_non_indexed==false && npf_indexed==false)
        throw new StrException("does not have normals per face");

    // default solid name
    char solidname[256] = "solidname";
    string ifs_name = ifs->getName();
    if(ifs_name.empty()==false) {
      snprintf(solidname,256,"%s",ifs_name.c_str());
    } else {
      // otherwise use filename, but first remove directory and extension
      // TODO
    }

    if(_fileType==SaverStl::FileType::ASCII) { ///////////////////////

      // if (all the conditions are satisfied) try to open the file
      fp = fopen(filename,"w");
      if( fp==(FILE*)0)
        throw new StrException("unable to open ASCII STL outputfile");

      if(_saveAscii(fp,solidname,*ifs)==false)
        throw new StrException("unable to save ASCII STL outputfile");
    
      fclose(fp);

    } else /* if(_fileType==FileType::BINARY) */ { ///////////////////

      // if (all the conditions are satisfied) try to open the file
      fp = fopen(filename,"wb");
      if( fp==(FILE*)0)
        throw new StrException("unable to open BINARY STL outputfile");

      if(_saveBinary(fp,solidname,*ifs)==false)
        throw new StrException("unable to save BINARY STL outputfile");

      fclose(fp);

    } ////////////////////////////////////////////////////////////////
    
    success = true;
    
  } catch(StrException* e) { 
    
    if(fp!=(FILE*)0) fclose(fp);
    fprintf(stderr,"SaverStl | ERROR | %s\n",e->what());
    delete e;
  }
  return success;
}
