#include "cross-correlation.hh"

using namespace std;

void cross_correlate_slow( const Signal & reference, const Signal & data, Signal & output )
{
    if ( reference.size() > data.size() ) {
        throw runtime_error( "reference length is longer than received data" );
    }

    /* get power of reference signal */
    float reference_power = 0;
    for ( const auto & x : reference ) {
        reference_power += norm( x );
    }

    /* cross-correlate for all positive lags, up to size of provided output signal */
    for ( unsigned int lag = 0; lag < output.size(); lag++ ) {
        complex<float> correlation = 0;

        for ( unsigned int i = 0; i < reference.size(); i++ ) {
            const auto & reference_sample = reference[ i ];
            const auto & data_sample = data[ lag + i ];

            correlation += data_sample * conj( reference_sample );
        }

        output[ lag ] = correlation / reference_power;
    }
}
