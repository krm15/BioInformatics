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
#include "KMerSource.h"
#include <itkVector.h>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "itkTimeProbe.h"


// Multipass algorithm that exhaustively searches the k-mer space
// using available computer memory with optimal time.
// Each k-mer encountered is broken down into two prefixes and one suffix.
// The first prefix is explored in a single pass through the data.
// The second prefix is an index into std::map array.
// The suffix is binned in the corresponding map.
// Upon finishing each pass, the top candidates are collected


// Convert a value in base-4 to an index
std::string GetString( size_t val, unsigned int kMerSize )
{
  std::string s;
  s.resize( kMerSize );

  size_t rem;

  for( unsigned int i = 0; i < kMerSize; i++ )
  {
    rem = val%4;
    switch( rem )
    {
      case 0 :
        s.insert( kMerSize-1-i, 1, 'A' );
        break;
      case 1 :
        s.insert( kMerSize-1-i, 1, 'C' );
        break;
      case 2 :
        s.insert( kMerSize-1-i, 1, 'G' );
        break;
      case 3 :
        s.insert( kMerSize-1-i, 1, 'T' );
        break;
    }
    val = val/4;
  }

  return s;
}

// Convert a string to an value in base-4
size_t GetIndex( std::string kmer, unsigned int kMerSize )
{
  size_t s = 0;
  unsigned int val;
  for( unsigned int i = 0; i < kMerSize; i++ )
  {
    switch(kmer[i])
    {
      case 'A' :
        val = 0;
        break;
      case 'C' :
        val = 1;
        break;
      case 'G' :
        val = 2;
        break;
      case 'T' :
        val = 3;
        break;
      default :
        val = 0;
    }

    size_t power(1);
    for( unsigned int j = 0; j < kMerSize-1-i; j++ )
    {
      power *= 4;
    }
    s += power * val;
  }

  return s;
}

// Function to flip std::pair
template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
    return std::pair<B,A>(p.second, p.first);
}

// Function to flip a std::map to a std::multimap
template<typename A, typename B>
std::multimap< B, A, std::greater<B> > flip_map(const std::map<A,B> &src)
{
    std::multimap<B,A, std::greater<B> > dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()),
                   flip_pair<A,B>);
    return dst;
}


int main ( int argc, char* argv[] )
{
  // Start timer to measure performance
  itk::TimeProbe cputimer;
  cputimer.Start();

  if ( argc < 4 )
  {
    std::cout << "Usage: " << std::endl;
    std::cout << argv[0] << " iFilename iKMerSize iTopCount <PassLength> <KMerPrefixLength>" << std::endl;
    return EXIT_FAILURE;
  }

  typedef std::map<std::string, size_t> KMerMapType;
  typedef KMerMapType::iterator KMerMapIteratorType;
  typedef std::multimap<size_t, std::string, std::greater< size_t > > AccumulatorMapType;
  typedef AccumulatorMapType::iterator AccumulatorMapIteratorType;

  // User-input parameter for k-mer length
  unsigned int iKMerSize = atoi( argv[2] );

  // User-input parameter for number of k-mer candidates
  unsigned int iTopCount = atoi( argv[3] );


  if ( iKMerSize == 0 )
  {
    std::cout << "K-Mer size needs to be larger than 0" << std::endl;
  }

  if ( iTopCount == 0 )
  {
    std::cout << "iTopCount size needs to be larger than 0" << std::endl;
  }

  // Prefix length to be explored in each pass of the data
  // Determine this from computer memory available and file size
  // Pass is 0 for small file sizes and large memories
  // 2 is suitable for file sizes of 15 GB and computer memory of 10 GB
  unsigned int PassLength = 2;
  if ( argc > 4 )
  {
    PassLength = atoi( argv[4] );
  }

  // Detemine prefix length
  unsigned int kMerSizePrefix = 0.5*(iKMerSize - PassLength);
  if ( argc > 5 )
  {
    kMerSizePrefix = atoi( argv[5] );
  }

  unsigned int kMerSizeSuffix = iKMerSize - PassLength - kMerSizePrefix;

  size_t counterSize(1);
  for( unsigned int j = 0; j < kMerSizePrefix; j++ )
  {
    counterSize *= 4;
  }

  size_t index;
  std::string line, kMer, prefix, suffix, passprefix;
  AccumulatorMapType KMerCounterAccumulator;
  std::fstream inFile( argv[1], std::ios::in );
  unsigned int len, p;

  if( !inFile )
  {
    std::cout << "File not opened..." << std::endl;
    return EXIT_FAILURE;
  }

  unsigned int numOfPasses(1);
  for( unsigned int j = 0; j < PassLength; j++ )
  {
    numOfPasses *= 4;
  }

  for( unsigned int pass = 0; pass < numOfPasses; pass++ )
  {
    std::cout << "Pass ID: " << pass << " of " << numOfPasses << std::endl;

    // Initialize an array of k-mers
    KMerMapType *KMerCounter;
    KMerCounter = new KMerMapType[ counterSize ];

    size_t readCounter(0);
    while( !inFile.eof() )
    {
      if ( readCounter%100000 == 0 )
      {
        std::cout << "Finished reading " << readCounter/100000 << "K lines" << std::endl;
      }

      // Read in first two lines
      std::getline(inFile, line);

      // Read second line and determine its length
      std::getline(inFile, line);
      len = std::strlen( line.c_str() );

      // length of the sequence needs to be larger or equal to the kmer size
      if ( len >= iKMerSize )
      {
        // Extract k-mers and insert into the appropriate std::map
        for( unsigned int i = 0; i < len-iKMerSize+1; i++ )
        {
          kMer = line.substr( i, iKMerSize );
          passprefix = kMer.substr( 0, PassLength );
          p =  GetIndex( passprefix, PassLength );

          if ( p == pass )
          {
            prefix = kMer.substr( PassLength, kMerSizePrefix);
            index = GetIndex( prefix, kMerSizePrefix );
            suffix = kMer.substr( PassLength+kMerSizePrefix, kMerSizeSuffix);
            KMerCounter[index][suffix]++;
          }
        }
      }

      std::getline(inFile, line);

      if( !inFile.eof() )
      {
        std::getline(inFile, line);
      }

      readCounter++;
    }

    // Flip KMerCounters into a multimap
    // Merge with KMerCounterAccumuator
    for( unsigned int i = 0; i < counterSize; i++ )
    {
      AccumulatorMapType fKMerCounter = flip_map<std::string, size_t>( KMerCounter[i] );
      KMerCounter[i].clear();

      // Retain iTopCount elements;
      AccumulatorMapIteratorType it = fKMerCounter.begin();
      unsigned countLimit = iTopCount;

      if ( countLimit > fKMerCounter.size() )
      {
        countLimit = fKMerCounter.size();
      }

      for( unsigned int count = 0; count < countLimit; count++ )
      {
        ++it;
      }
     fKMerCounter.erase( it, fKMerCounter.end() );


     // Add to KMerCounterAccumulator and delete kMerCounter
      for ( it = fKMerCounter.begin(); it != fKMerCounter.end(); ++it)
      {
        kMer = GetString( pass, PassLength ) +
                           GetString( i, kMerSizePrefix ) +
                           (*it).second;
        KMerCounterAccumulator.insert( std::make_pair( (*it).first, kMer ) );
      }
      fKMerCounter.clear();
    }

    delete[] KMerCounter;

    // Accumulate only the iTopCount values
    unsigned int countLimit = iTopCount;

    if ( countLimit > KMerCounterAccumulator.size() )
    {
      countLimit = KMerCounterAccumulator.size();
    }

    AccumulatorMapIteratorType it = KMerCounterAccumulator.begin();
    for( unsigned int count = 0; count < countLimit; count++ )
    {
      ++it;
    }
    KMerCounterAccumulator.erase( it, KMerCounterAccumulator.end() );

    // Reset file stream to beginning
    inFile.clear();
    inFile.seekg(0, std::ios::beg);
  }
  inFile.close();


  // Write out the top 25 k-mers
  std::cout << "Print accumulator size " << std::endl;
  AccumulatorMapIteratorType it = KMerCounterAccumulator.begin();
  while ( it != KMerCounterAccumulator.end() )
  {
    std::cout << (*it).first << ' ' << (*it).second << std::endl;
    ++it;
  }

  cputimer.Stop();
  std::cout << "K-mer extraction took " << cputimer.GetMean() << " seconds" << std::endl;

  return EXIT_SUCCESS;
}
