#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <complex>
#include <numeric>
#include <type_traits>
#include <boost/align/aligned_allocator.hpp>

#include <fftw3.h> /* must be included after <complex> */

#include "helpers.hh"

using namespace std;

using Signal = vector<complex<float>, boost::alignment::aligned_allocator<complex<float>, 64>>;

Signal read( const DAT & dat_file )
{
    Signal ret;
    ret.reserve( dat_file.IQ_sample_count() );
    for ( unsigned int i = 0; i < dat_file.IQ_sample_count(); i++ ) {
        ret.emplace_back( dat_file.I( i ), dat_file.Q( i ) );
    }

    return ret;
}

const float sample_rate = 15.36 * 1.0e6;

void program_body( const string & reference_filename, const string & data_filename )
{
    Signal reference = read( reference_filename ), data = read( data_filename );

    if ( reference.size() > data.size() ) {
        throw runtime_error( "reference length is longer than received data" );
    }

    /* pad the reference */
    reference.resize( data.size() );

    const fftw_plan forward = fftw_plan_dft_1d( data.size(), in, out, FFTW_FORWARD );
    const fftw_plan reverse = fftw_plan_dft_1d( data.size(), in, out, FFTW_BACKWARD );
}

int main( const int argc, const char * argv[] )
{
    if ( argc < 0 ) { abort(); }

    if ( argc != 3 ) {
        cerr << "Usage: " << argv[ 0 ] << " reference data\n";
        return EXIT_FAILURE;
    }

    try {
        program_body( argv[ 1 ], argv[ 2 ] );
    } catch ( const exception & e ) {
        cerr << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
