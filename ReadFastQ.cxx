#include <itkVector.h>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "itkTimeProbe.h"


class KMer
{
public:
  KMer(int id)
  {
    iKMerSize = id;
    AC_GT = new bool[iKMerSize];
    AG_CT = new bool[iKMerSize];

    //AC_GT.resize(iKMerSize);
    //AG_CT.resize(iKMerSize);
  }

  KMer(std::string p)
  {
    iKMerSize = p.size();
    AC_GT = new bool[iKMerSize];
    AG_CT = new bool[iKMerSize];

    //AC_GT.resize(iKMerSize);
    //AG_CT.resize(iKMerSize);

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

  ~KMer()
  {
    //AC_GT.clear();
    //AG_CT.clear();

    //delete[] AC_GT;
    //delete[] AG_CT;
  }

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

  bool operator <(const KMer& rhs) const
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
  bool *AC_GT; // A (0) or C (1)
  bool *AG_CT; // G (2) or T (3)
};


int main ( int argc, char* argv[] )
{
  if ( argc < 4 )
  {
    std::cout << "Usage: " << std::endl;
    std::cout << argv[0] << " iFilename iKMerSize iTopCount <iLargeReadFile>" << std::endl;
    return EXIT_FAILURE;
  }

  typedef std::map<KMer, unsigned int> KMerMapType;
  typedef KMerMapType::iterator KMerMapIteratorType;

  unsigned int iKMerSize = atoi( argv[2] );
  unsigned int iTopCount = atoi( argv[3] );

  bool iLargeRead;

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

  KMerMapType KMerCounter;
  std::string line, kmer;
  std::fstream inFile( argv[1], std::ios::in );

  itk::TimeProbe cputimer;
  cputimer.Start();

  if( !inFile )
  {
    std::cout << "File not opened..." << std::endl;
    return EXIT_FAILURE;
  }

  unsigned int cc = 0;

  unsigned int len;

  if ( ! iLargeRead )
  {
    // Use unordered maps for small reads
    while( !inFile.eof() )
    {
      if ( cc%100000 == 0 )
      {
        std::cout << cc/100000 << std::endl;
      }
      cc++;

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
          KMer km( line.substr( i, iKMerSize) );
          KMerCounter[km]++;
          //std::cout << km.GetNTide() << ' ' << KMerCounter.size() << std::endl;
        }
      }

      std::getline(inFile, line);

      if( !inFile.eof() )
      {
        std::getline(inFile, line);
      }
    }
  }
  else // For large reads
  {
    char ch;
    kmer.reserve( iKMerSize );

    while( !inFile.eof() )
    {
      // Read first line -- usually small
      std::getline(inFile, line);

      // Streaming read till end-of-line is encountered
      for( unsigned int i = 0; i < iKMerSize-1; i++ )
      {
        if ( !inFile.eof() )
        {
          inFile >> ch;
          kmer.insert( i, 1, ch );
        }
      }

      if ( !inFile.eof() )
      {
        inFile >> ch;
        kmer.insert( iKMerSize - 1, 1, ch );
      }

      while ( ( ch != '\n' ) && ( !inFile.eof() ) )
      {
        KMer km(kmer);
        KMerCounter[km]++;

        // Read next character
        ch = inFile.get();

        // Move everything in the string up by one character
        // Append last character
        std::string ss(1, ch);
        kmer = kmer.substr(1, iKMerSize-1) + ss;
      }

      // Read third line -- usually small
      std::getline(inFile, line);

      ch = '0';
      // Read fourth line - streaming
      while ( ( ch != '\n' ) && !inFile.eof() )
      {
        ch = inFile.get();
      }

      kmer.clear();
    }
  }

  inFile.close();
  std::cout << "Finished reading file " << cc << std::endl;

  cputimer.Stop();
  std::cout << "Reading took " << cputimer.GetMean() << " seconds" << std::endl;

  // Adjust iTopCount to the size of the map, only if its smaller
  unsigned int size = KMerCounter.size();
  if (size < iTopCount)
  {
    iTopCount = size;
  }

  // Print out the top counts
  unsigned int count = 0;
  std::cout << "Top " << iTopCount << " K-Mers are: " << std::endl;
  for (KMerMapIteratorType iter = KMerCounter.begin(); count < iTopCount; iter++)
  {
    KMer km = iter->first;
    std::cout << km.GetNTide() << ' ' << iter->second << std::endl;
    count++;
  }

  KMerMapIteratorType iter = KMerCounter.begin();
  while ( iter != KMerCounter.end() )
  {
    KMerMapIteratorType toErase = iter;
    ++iter;
    KMerCounter.erase(toErase);
  }

  return EXIT_SUCCESS;
}
