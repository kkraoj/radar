#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <complex>
#include "cross-correlation.hh"

using namespace std;

void program_body( const string & reference_filename, const string & data_filename )
{
    /* read in DAT files */
    DAT reference_dat { reference_filename };
    Signal reference( reference_dat.IQ_sample_count() );
    reference_dat.read( 0, reference );

    DAT data_dat { data_filename };
    Signal data( data_dat.IQ_sample_count() );
    data_dat.read( 0, data );

    /* make output */
    Signal crosscorrelation( data.size() - reference.size() );

    /* do the cross-correlation with a slow O(N^2) algorithm */
    cross_correlate_slow( reference, data, crosscorrelation );

    /* print */
    const float sample_rate = 15.36 * 1.0e6;
    for ( unsigned int lag = 0; lag < crosscorrelation.size(); lag++ ) {
        cout << lag / sample_rate << " " << abs( crosscorrelation[ lag ] ) << "\n";
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
