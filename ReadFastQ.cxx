#include <itkVector.h>
#include <fstream>
#include <string>
#include <sstream>

int main ( int argc, char* argv[] )
{
  if ( argc < 3 )
  {
    std::cout << "Usage: " << std::endl;
    std::cout << argv[0] << " iFilename iKMerSize iTopCount" << std::endl;
    return EXIT_FAILURE;
  }

  typedef std::map<std::string, unsigned int> MapType;
  typedef MapType::const_iterator MapIterator;

  unsigned int kmersize = atoi( argv[2] );
  unsigned int topcount = atoi( argv[3] );

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

  if( !inFile )
  {
    std::cout << "File not opened..." << std::endl;
    return EXIT_FAILURE;
  }

  while( !inFile.eof() )
  {
    std::getline(inFile, line);

    std::getline(inFile, line);
    unsigned len = line.length();

    for( unsigned int i = iMerSize -1; i <  )
    {
      MapType::iterator it = KMerCounter.find( kmer );
      if ( it != KMerCounter.end() )
      {
        KMerCounter[kmer]++;
      }
      else
      {
        KMerCounter[kmer] = 1;
      }
    }

    std::getline(inFile, line);

    std::getline(inFile, line);
  }

  inFile.close();


  // Print out the top counts
  unsigned int count = 0;
  for (MapIterator iter = KMerCounter.begin(); count < topcount; iter++)
  {
    std::cout << iter->first << ' ' << iter->second << std::endl;
    count++;
  }

  return EXIT_SUCCESS;
}

