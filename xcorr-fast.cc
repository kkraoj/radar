#include <cstdlib>
#include <iostream>
#include <thread>

#include "helpers.hh"
#include "cross-correlation.hh"

using namespace std;

void output_line( const float x, const float y )
{
    const float line[ 2 ] = { x, y };
    if ( 2 != fwrite( line, sizeof( float ), 2, stdout ) ) {
        throw unix_error( "fwrite" );
    }
}

void program_body( const string & reference_filename, const string & data_filename )
{
    /* read in DAT files */
    DAT reference_dat { reference_filename };
    Signal reference( reference_dat.IQ_sample_count() );
    reference_dat.read( 0, reference );

    DAT data_dat { data_filename };
    Signal data( data_dat.IQ_sample_count() );
    data_dat.read( 0, data );

    if ( reference.size() > data.size() ) {
        throw runtime_error( "reference length is longer than received data" );
    }

    /* make output */
    vector<float> result( data.size() - reference.size() );

    CrossCorrelator cross_correlator( reference.size(), data.size(), 7680000 );
    cross_correlator.correlate_fast( reference, data, result );

    /* print */
    const float sample_rate = 15.36 * 1.0e6;
    bool printing = false;
    for ( unsigned int lag = 0; lag < result.size(); lag++ ) {
        if ( result[ lag ] > 0.001 ) {
            if ( not printing and lag > 0 ) {
                output_line( (lag-1) / sample_rate, 0 );
                printing = true;
            }

            output_line( lag / sample_rate, result[ lag ] );
        } else {
            if ( printing ) {
                output_line( lag / sample_rate, 0 );
                printing = false;
            }
        }
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
