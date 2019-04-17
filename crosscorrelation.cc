#include <cstdlib>
#include <cstring>
#include <iostream>

#include "helpers.hh"

using namespace std;

void cross_correlate( const string & reference_filename, const string & data_filename )
{
    DAT reference { reference_filename };
    DAT data { data_filename };

    if ( reference.IQ_sample_count() > data.IQ_sample_count() ) {
        throw runtime_error( "reference length is longer than received data" );
    }

    const double sample_rate = 15.36 * 1.0e6;

    for ( unsigned int i = 0; i < data.IQ_sample_count() - reference.IQ_sample_count(); i++ ) {
        double correlation = 0;

        for ( unsigned int j = 0; j < reference.IQ_sample_count(); j++ ) {
            correlation += data.I( i + j ) * reference.I( j ) + data.Q( i + j ) * reference.Q( j );
        }

        cout << i / sample_rate << " " << correlation << "\n";
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
        cross_correlate( argv[ 1 ], argv[ 2 ] );
    } catch ( const exception & e ) {
        cerr << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
