#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <complex>
#include <numeric>
#include <type_traits>
#include <thread>
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

void program_body( const string & reference_filename, const string & data_filename )
{
    Signal reference = read( reference_filename ), data = read( data_filename );

    if ( reference.size() > data.size() ) {
        throw runtime_error( "reference length is longer than received data" );
    }

    /* zero-pad data out to 2x its length */
    data.resize( data.size() * 10 );

    /* pad the reference */
    reference.resize( data.size() );

    Signal reference_fft( reference.size() );
    Signal data_fft( data.size() );
    Signal crosscorrelation( reference.size() );

    const fftwf_plan ref_forward = fftwf_plan_dft_1d( reference.size(),
                                                      reinterpret_cast<fftwf_complex *>( reference.data() ),
                                                      reinterpret_cast<fftwf_complex *>( reference_fft.data() ),
                                                      FFTW_FORWARD,
                                                      FFTW_ESTIMATE );

    const fftwf_plan data_forward = fftwf_plan_dft_1d( reference.size(),
                                                       reinterpret_cast<fftwf_complex *>( data.data() ),
                                                       reinterpret_cast<fftwf_complex *>( data_fft.data() ),
                                                       FFTW_FORWARD,
                                                       FFTW_ESTIMATE );

    const fftwf_plan reverse = fftwf_plan_dft_1d( reference.size(),
                                                  reinterpret_cast<fftwf_complex *>( reference_fft.data() ),
                                                  reinterpret_cast<fftwf_complex *>( crosscorrelation.data() ),
                                                  FFTW_BACKWARD,
                                                  FFTW_ESTIMATE );

    for ( auto & x : reference ) {
        x = 0;
    }

    for ( auto & x : data ) {
        x = 0;
    }

    reference.at( 24 ) = 1.0;

    data.at( 23 ) = 1.0;

    thread fft1 { [&ref_forward] { fftwf_execute( ref_forward ); } };
    thread fft2 { [&data_forward] { fftwf_execute( data_forward ); } };

    fft1.join();
    fft2.join();

    /*
    for ( unsigned int i = 0; i < reference_fft.size(); i++ ) {
        cout << i << " " << reference_fft.at( i ).real() << " " << reference_fft.at( i ).imag() << "\n";
    }

    return;
    */

    /* multiply */
    float reference_power = 0;
    for ( unsigned int i = 0; i < reference_fft.size(); i++ ) {
        reference_power += norm( reference_fft[ i ] );
        reference_fft[ i ] *= conj( data_fft[ i ] );
    }

    fftwf_execute( reverse );

    /* print */
    const float sample_rate = 15.36 * 1.0e6;
    for ( unsigned int lag = 0; lag < crosscorrelation.size(); lag++ ) {
        cout << lag / (2 * sample_rate) << " " << abs( crosscorrelation[ lag ] ) / reference_power << "\n";
    }
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
