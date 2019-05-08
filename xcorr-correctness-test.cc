#include <cstdlib>
#include <iostream>
#include <random>

#include "helpers.hh"
#include "cross-correlation.hh"

using namespace std;

const unsigned int ITERATIONS = 100;
const unsigned int MAX_REFERENCE = 1000;
const unsigned int MAX_DATA_MINUS_REFERENCE = 20000;

void program_body()
{
    uniform_int_distribution<> reference_size { 32, MAX_REFERENCE };
    uniform_int_distribution<> data_minus_reference_size { 1, MAX_DATA_MINUS_REFERENCE };
    uniform_real_distribution<> sample { -32768.0, 32768.0 };

    default_random_engine rng( random_device{}() );

    for ( unsigned int iteration = 0; iteration < ITERATIONS; iteration++ ) {
        Signal reference( reference_size( rng ) );
        Signal data( reference.size() + data_minus_reference_size( rng ) );

        for ( auto & x : reference ) {
            x = sample( rng );
        }

        for ( auto & x : data ) {
            x = sample( rng );
        }

        if ( iteration % 2 ) {
            /* add the reference in there somewhere */
            const unsigned int offset = reference_size( rng ) + data_minus_reference_size( rng );
            for ( unsigned int index = 0;
                  index < reference.size() and offset + index < data.size();
                  index++ ) {
                data.at( offset + index ) += reference.at( index );
            }
        }

        /* make output */
        vector<float> result_slow( data.size() - reference.size() ),
            result_fast( data.size() - reference.size() );

        correlate_slow( reference, data, result_slow );

        CrossCorrelator cross_correlator( reference.size(), data.size() );
        cross_correlator.correlate_fast( reference, data, result_fast );

        for ( unsigned int lag = 0; lag < result_slow.size(); lag++ ) {
            if ( abs( result_fast[ lag ] - result_slow[ lag ] ) > 0.0001 ) {
                throw runtime_error( "test failure at reference_size=" + to_string( reference.size() )
                                     + " data_size=" + to_string( data.size() )
                                     + " lag=" + to_string( lag )
                                     + " slow=" + to_string( result_slow[ lag ] )
                                     + " vs. fast=" + to_string( result_fast[ lag ] ) );
            }
        }
    }
}

int main( const int argc, const char * argv[] )
{
    if ( argc < 0 ) { abort(); }

    if ( argc != 1 ) {
        cerr << "Usage: " << argv[ 0 ] << "\n";
        return EXIT_FAILURE;
    }

    try {
        program_body();
    } catch ( const exception & e ) {
        cerr << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
