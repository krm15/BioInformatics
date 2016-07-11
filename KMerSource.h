/*=========================================================================
 Authors: Kishore Mosaliganti
 at Megason Lab, Systems biology, Harvard Medical school, 2009

 Copyright (c) 2009, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#ifndef __KMerSource_h
#define __KMerSource_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>
#include <cmath>



class KMerSource
{
public:
  KMerSource(int id)
  {
    iKMerSize = id;
    AC_GT = new bool[iKMerSize];
    AG_CT = new bool[iKMerSize];
  }

  KMerSource(std::string p)
  {
    iKMerSize = p.size();
    AC_GT = new bool[iKMerSize];
    AG_CT = new bool[iKMerSize];

    for( unsigned int i = 0; i < iKMerSize; i++ )
    {
      switch(p[i])
      {
        case 'A' :
          AC_GT[i] = 0;
          AG_CT[i] = 0;
          break;
        case 'C' :
          AC_GT[i] = 0;
          AG_CT[i] = 1;
          break;
        case 'G' :
          AC_GT[i] = 1;
          AG_CT[i] = 0;
          break;
        case 'T' :
          AC_GT[i] = 1;
          AG_CT[i] = 1;
      }
    }
  }

  ~KMerSource(){}

  const std::string GetNTide()
  {
    std::string p;
    p.reserve( iKMerSize );
    for( unsigned int i = 0; i < iKMerSize; i++ )
    {
      switch( 2 * int(AC_GT[i]) + AG_CT[i] )
      {
        case 0 :
          p.insert( i, 1, 'A' );
          break;
        case 1 :
          p.insert( i, 1, 'C' );
          break;
        case 2 :
          p.insert( i, 1, 'G' );
          break;
        case 3 :
          p.insert( i, 1, 'T' );
          break;
      }
    }

    return p;
  }

  bool operator <(const KMerSource& rhs) const
  {
    for( unsigned int i = 0; i < iKMerSize; i++ )
    {
      if ( 2 * int(AC_GT[i]) + AG_CT[i] < 2 * int(rhs.AC_GT[i]) + rhs.AG_CT[i] )
      {
        return true;
      }
      else if ( 2 * int(AC_GT[i]) + AG_CT[i] > 2 * int(rhs.AC_GT[i]) + rhs.AG_CT[i] )
      {
        return false;
      }
    }
    return false;
  }

private:
  unsigned char iKMerSize;
  bool *AC_GT;
  bool *AG_CT;
};

#endif
