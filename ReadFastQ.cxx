#include <itkVector.h>
#include <fstream>
#include <string>
#include <sstream>

#include "itkTimeProbe.h"

int main ( int argc, char* argv[] )
{
  if ( argc < 4 )
  {
    std::cout << "Usage: " << std::endl;
    std::cout << argv[0] << " iFilename iKMerSize iTopCount <iLargeReadFile>" << std::endl;
    return EXIT_FAILURE;
  }

  typedef std::pair<std::string,unsigned int> PairType;
  typedef std::map<std::string, unsigned int> MapType;
  typedef MapType::const_iterator MapIteratorType;
  typedef std::pair<MapIteratorType,bool> MapInsertReturnType;

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

  MapType KMerCounter;
  std::string line, kmer;
  std::fstream inFile( argv[1], std::ios::in );

  itk::TimeProbe cputimer;
  cputimer.Start();

  if( !inFile )
  {
    std::cout << "File not opened..." << std::endl;
    return EXIT_FAILURE;
  }

  unsigned int len;

  if ( ! iLargeRead )
  {
    // Use unordered maps for small reads
    while( !inFile.eof() )
    {
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
          kmer = line.substr( i, iKMerSize);
          KMerCounter[kmer]++;
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
        KMerCounter[kmer]++;

        // Read next character
        ch = inFile.get();

        // Move everything in the string up by one character
        // Append last character
        std::string ss(1, ch);
        kmer = kmer.substr(1,iKMerSize-1) + ss;
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
  std::cout << "Finished reading file" << std::endl;

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
  for (MapIteratorType iter = KMerCounter.begin(); count < iTopCount; iter++)
  {
    std::cout << iter->first << ' ' << iter->second << std::endl;
    count++;
  }

  return EXIT_SUCCESS;
}
