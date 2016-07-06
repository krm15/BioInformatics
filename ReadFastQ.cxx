#include <itkVector.h>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "itkTimeProbe.h"


std::string GetString( size_t val, unsigned int iKMerSizep )
{
  std::string s;
  s.reserve( iKMerSizep );

  unsigned int rem;

  for( unsigned int i = 0; i < iKMerSizep; i++ )
  {
    rem = val%4;
    val = val/4;

    switch( rem )
    {
      case 0 :
        s.insert( iKMerSizep-1-i, 1, 'A' );
        break;
      case 1 :
        s.insert( iKMerSizep-1-i, 1, 'C' );
        break;
      case 2 :
        s.insert( iKMerSizep-1-i, 1, 'G' );
        break;
      case 3 :
        s.insert( iKMerSizep-1-i, 1, 'T' );
        break;
    }
  }

  return s;
}

size_t GetIndex( std::string kmer, unsigned int iKMerSizep )
{
  size_t s = 0;
  unsigned int val;
  for( unsigned int i = 0; i < iKMerSizep; i++ )
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

    s += size_t(pow(4, iKMerSizep-1-i)) * val;
  }

  return s;
}


template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
    return std::pair<B,A>(p.second, p.first);
}

template<typename A, typename B>
std::multimap<B,A> flip_map(const std::map<A,B> &src)
{
    std::multimap<B,A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()),
                   flip_pair<A,B>);
    return dst;
}


int main ( int argc, char* argv[] )
{
  if ( argc < 4 )
  {
    std::cout << "Usage: " << std::endl;
    std::cout << argv[0] << " iFilename iKMerSize iTopCount <iLargeReadFile>" << std::endl;
    return EXIT_FAILURE;
  }

  typedef std::map<std::string, unsigned int> KMerMapType;
  typedef KMerMapType::iterator KMerMapIteratorType;
  typedef std::multimap<unsigned int, std::string> AccumulatorMapType;
  typedef AccumulatorMapType::iterator AccumulatorMapIteratorType;

  unsigned int iKMerSize = atoi( argv[2] );
  unsigned int iTopCount = atoi( argv[3] );

  bool iLargeRead;
  unsigned int PrefixLengthInPass = 2;

  // Determine this from computer memory and file size
  iLargeRead = false;

  if ( argc > 4 )
  {
    iLargeRead = atoi( argv[4] );
  }

  if ( iKMerSize == 0 )
  {
    std::cout << "K-Mer size needs to be larger than 0" << std::endl;
  }

  if ( iTopCount == 0 )
  {
    std::cout << "iTopCount size needs to be larger than 0" << std::endl;
  }

  size_t index;
  unsigned int iKMerSizep = 0.5*(iKMerSize - PrefixLengthInPass);
  unsigned int iKMerSizes = iKMerSize - PrefixLengthInPass - iKMerSizep;
  size_t counterSize = size_t( pow(4, iKMerSizep) );

  std::string line, kmer;
  std::fstream inFile( argv[1], std::ios::in );

  itk::TimeProbe cputimer; 
  cputimer.Start();

  if( !inFile )
  {
    std::cout << "File not opened..." << std::endl;
    return EXIT_FAILURE;
  }

  AccumulatorMapType KMerCounterAccumulator;

  unsigned int len;
  unsigned int numOfPasses =  int( pow(4, PrefixLengthInPass) );

  for( unsigned int pass = 0; pass < numOfPasses; pass++ )
  {
    KMerMapType *KMerCounter;
    KMerCounter = new KMerMapType[ counterSize ];
    unsigned int cc = 0;
    while( !inFile.eof() )
    {
      if ( cc%100000 == 0 )
      {
        std::cout << cc/100000 << ' ' << KMerCounter[0].size() << std::endl;
      }

      std::getline(inFile, line);
      //std::cout << line << std::endl;

      std::getline(inFile, line);
      //std::cout << line << std::endl;

      len = std::strlen( line.c_str() );
      //std::cout << len << std::endl;
      if ( len >= iKMerSize )
      {
        for( unsigned int i = 0; i < len - iKMerSize; i++ )
        {
          //KMer km( line.substr( i, iKMerSize) );
          if ( GetIndex(line.substr( i, PrefixLengthInPass), PrefixLengthInPass) == pass )
          {
            index = GetIndex( line.substr( i+PrefixLengthInPass, iKMerSizep), iKMerSizep );
            KMerCounter[index][ line.substr( i+PrefixLengthInPass+iKMerSizep, iKMerSizes) ]++;//km
            //std::cout << km.GetNTide() << ' ' << KMerCounter.size() << std::endl;
          }
        }
      }

      std::getline(inFile, line);

      if( !inFile.eof() )
      {
        std::getline(inFile, line);
      }
      cc++;
    }

    // Merge KMerCounter and KMerCounterAccumuator
    for(unsigned int i = 0; i < counterSize; i++)
    {
      AccumulatorMapType fKMerCounter = flip_map<std::string, unsigned int>( KMerCounter[i] );
      KMerCounter[i].clear();

      // Retain iTopCount elements;
      AccumulatorMapIteratorType it = fKMerCounter.begin();
      unsigned int count = 0;
      unsigned countLimit = iTopCount;

      if ( countLimit > fKMerCounter.size() )
      {
        countLimit = fKMerCounter.size();
      }

      while( count < countLimit )
      {
        count++;
        ++it;
      }
     fKMerCounter.erase( it, fKMerCounter.end() );

     std::string kmer = GetString( pass, PrefixLengthInPass ) +
                        GetString( counterSize, iKMerSizep );

     // Add to KMerCounterAccumulator
      for (it=fKMerCounter.begin(); it!=fKMerCounter.end(); ++it)
      {
        kmer += (*it).second;
        KMerCounterAccumulator.insert( std::make_pair( (*it).first, kmer ) );
      }
      fKMerCounter.clear();
    }

    AccumulatorMapIteratorType it = KMerCounterAccumulator.begin();
    unsigned int count = 0;
    unsigned int countLimit = iTopCount;

    if ( countLimit > KMerCounterAccumulator.size() )
    {
      countLimit = KMerCounterAccumulator.size();
    }

    while( count < countLimit )
    {
      count++;
      ++it;
    }
    KMerCounterAccumulator.erase( it, KMerCounterAccumulator.end() );

    delete[] KMerCounter;

    inFile.clear();
    inFile.seekg(0, std::ios::beg);
  }

  inFile.close();
  std::cout << "Finished reading file " << std::endl;

  // Write out the top 25 k-mers

  cputimer.Stop();
  std::cout << "Reading took " << cputimer.GetMean() << " seconds" << std::endl;

  return EXIT_SUCCESS;
}
