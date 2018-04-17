#ifndef RW_BINARY_FILE_HH
#define RW_BINARY_FILE_HH

#include <string>
#include <fstream>

template< class a_type >
int read_from_raw_file( std::string filename, a_type* data, int length ) {
  std::ifstream raw_input( filename.c_str(), std::ios::binary );
  //check if file is open
  if ( raw_input.fail() )
    {
      std::cerr << "something is wrong with the input file." << std::endl;
      //exit(-1);
      return 0;
    };
  int i(0);
  while ( raw_input.good() )
    { 
      //check if data is full
      if ( i == length )
	{
	  std::cerr << "input file was to long." << std::endl;
	  exit(-1);
	};
      raw_input.read( reinterpret_cast< char* >( &data[i] ), sizeof data[i] );
      //check if there is more data or end of file
      raw_input.peek();
      i++;
    };
  //check if file contained enough data
  if ( i != length )
    {
      std::cerr << "input file was to short." << std::endl;
      exit(-1);
    };
  raw_input.close();
  return 1;
};

template< class a_type >
void write_to_raw_file( std::string filename, const a_type* data, int length ) {
  std::ofstream output( filename.c_str(), std::ios::binary );
  //check if file is open
  if ( output.fail() )
    {
      std::cerr << "something is wrong with the output file." << std::endl;
      exit(-1);
    };
  for ( int i(0); i < length; i++ )
    output.write( reinterpret_cast< const char* >( &data[i] ), sizeof data[i] );
};

#endif
