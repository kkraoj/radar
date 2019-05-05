#include <cstdlib>
#include <cstring>
#include <iostream>

#include "helpers.hh"

using namespace std;

void program_body( const string & filename )
{
    DAT dat { filename };

    const double sample_rate = 15.36 * 1.0e6;
    const uint64_t interval = sample_rate / 1000;

    cout.precision( 10 );

    for ( uint64_t interval_start = 0; interval_start < dat.IQ_sample_count(); interval_start += interval ) {
        uint64_t interval_end = min( interval_start + interval, dat.IQ_sample_count() );

        double power = 0, I_power = 0, Q_power = 0;

        for ( unsigned int i = interval_start; i < interval_end; i++ ) {
            int16_t I_sample = dat.I( i );
            int16_t Q_sample = dat.Q( i );
            power += I_sample * I_sample + Q_sample * Q_sample;
            I_power += I_sample * I_sample;
            Q_power += Q_sample * Q_sample;
        }

        cout << ((interval_start + interval_end) / 2) / sample_rate << " " << power / interval << " " << I_power / interval << " " << Q_power / interval << "\n";
    }
}

int main( const int argc, const char * argv[] )
{
    if ( argc < 0 ) { abort(); }

    if ( argc != 2 ) {
        cerr << "Usage: " << argv[ 0 ] << " filename\n";
        return EXIT_FAILURE;
    }

    try {
        program_body( argv[ 1 ] );
    } catch ( const exception & e ) {
        cerr << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
