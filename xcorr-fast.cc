#include <cstdlib>
#include <iostream>
#include <thread>

#include "helpers.hh"

#include <fftw3.h> /* must be included last, to come after <complex> */

using namespace std;

void program_body( const string & reference_filename, const string & data_filename )
{
    DAT reference_dat { reference_filename }, data_dat { data_filename };
    
    if ( reference_dat.IQ_sample_count() > data_dat.IQ_sample_count() ) {
        throw runtime_error( "reference length is longer than received data" );
    }

    /* set up signals */
    Signal data( data_dat.IQ_sample_count() * 2 );
    Signal reference( data.size() );

    Signal reference_fft( data.size() );
    Signal data_fft( data.size() );
    Signal crosscorrelation( data.size() );

    const fftwf_plan ref_forward = fftwf_plan_dft_1d( data.size(),
                                                      reinterpret_cast<fftwf_complex *>( reference.data() ),
                                                      reinterpret_cast<fftwf_complex *>( reference_fft.data() ),
                                                      FFTW_FORWARD,
                                                      FFTW_ESTIMATE );

    const fftwf_plan data_forward = fftwf_plan_dft_1d( data.size(),
                                                       reinterpret_cast<fftwf_complex *>( data.data() ),
                                                       reinterpret_cast<fftwf_complex *>( data_fft.data() ),
                                                       FFTW_FORWARD,
                                                       FFTW_ESTIMATE );

    const fftwf_plan reverse = fftwf_plan_dft_1d( data.size(),
                                                  reinterpret_cast<fftwf_complex *>( data_fft.data() ),
                                                  reinterpret_cast<fftwf_complex *>( crosscorrelation.data() ),
                                                  FFTW_BACKWARD,
                                                  FFTW_ESTIMATE );

    reference_dat.read( 0, reference );
    data_dat.read( 0, data );

    thread fft1 { [&ref_forward] { fftwf_execute( ref_forward ); } };
    thread fft2 { [&data_forward] { fftwf_execute( data_forward ); } };

    fft1.join();
    fft2.join();

    /* multiply */
    float reference_power = 0;
    for ( unsigned int i = 0; i < data_fft.size(); i++ ) {
        reference_power += norm( reference_fft[ i ] );
        data_fft[ i ] *= conj( reference_fft[ i ] );
    }

    fftwf_execute( reverse );

    /* print */
    const float sample_rate = 15.36 * 1.0e6;
    for ( unsigned int lag = 0; lag < crosscorrelation.size() / 2; lag++ ) {
        cout << lag / sample_rate << " " << abs( crosscorrelation[ lag ] ) / reference_power << "\n";
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
