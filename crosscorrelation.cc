#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <complex>
#include <numeric>

#include "helpers.hh"

using namespace std;

vector<complex<float>> read( const DAT & dat_file )
{
    vector<complex<float>> ret;
    ret.reserve( dat_file.IQ_sample_count() );
    for ( unsigned int i = 0; i < dat_file.IQ_sample_count(); i++ ) {
        ret.emplace_back( dat_file.I( i ), dat_file.Q( i ) );
    }
    return ret;
}

void program_body( const string & reference_filename, const string & data_filename )
{
    vector<complex<float>> reference = read( reference_filename ), data = read( data_filename );

    if ( reference.size() > data.size() ) {
        throw runtime_error( "reference length is longer than received data" );
    }

    /* conjugate reference signal, and get total power of reference */
    float reference_power = 0;
    for ( auto & x : reference ) {
        x = conj( x );
        reference_power += norm( x );
    }

    const float sample_rate = 15.36 * 1.0e6;

    unsigned int last_percent_complete = -1;

    for ( unsigned int lag = 0; lag < data.size() - reference.size(); lag++ ) {
        /* report status */
        const unsigned int percent_complete = 10000 * lag / (data.size() - reference.size());

        if ( percent_complete != last_percent_complete ) {
            cerr << "\r" << percent_complete / 100.0 << "%                ";
            last_percent_complete = percent_complete;
        }

        /* calculate correlation */
        complex<float> correlation = 0;

        for ( unsigned int i = 0; i < reference.size(); i++ ) {

            const auto & data_sample = data[ lag + i ];
            const auto & reference_sample = reference[ i ];

            correlation += data_sample * reference_sample;
        }

        cout << lag / sample_rate << " " << abs( correlation ) / reference_power << "\n";
    }

    cerr << "\n";
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
