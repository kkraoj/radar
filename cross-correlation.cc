#include "cross-correlation.hh"

#include <iostream>
#include <cstring>
#include <thread>

using namespace std;

void correlate_slow( const Signal & reference, const Signal & data, vector<float> & output )
{
    if ( reference.size() > data.size() ) {
        throw runtime_error( "reference length is longer than received data" );
    }

    if ( output.size() != data.size() - reference.size() ) {
        throw runtime_error( "invalid output length (must be data_length - reference_length)" );
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

        output[ lag ] = abs( correlation ) / reference_power;
    }
}

/* make sure that global FFTW state is cleaned up when program exits */
class FFTW
{
public:
    FFTW()
    {
        if ( not fftwf_init_threads() ) {
            throw unix_error( "fftwf_init_threads" );
        }

        if ( thread::hardware_concurrency() > 1 ) {
            fftwf_plan_with_nthreads( thread::hardware_concurrency() / 2 );
        }
    }

    ~FFTW() { fftwf_cleanup_threads(); }
};

FFTW global_fftw_state;

template <typename T>
T * notnull( const string & what, T * x ) { if ( x ) { return x; } throw runtime_error( what ); }

FFTPlan::FFTPlan( Signal & input, Signal & output,
                  const int sign, const int flags )
    : plan_( notnull( "fftwf_plan_dft_1d",
                      fftwf_plan_dft_1d( input.size(),
                                         reinterpret_cast<fftwf_complex *>( input.data() ),
                                         reinterpret_cast<fftwf_complex *>( output.data() ),
                                         sign, flags ) ) )
{}

FFTPlan::~FFTPlan()
{
    fftwf_destroy_plan( plan_ );
}

void FFTPlan::execute()
{
    fftwf_execute( plan_ );
}

CrossCorrelator::CrossCorrelator( const size_t reference_length,
                                  const size_t data_length,
                                  const size_t max_chunk_size )
    : reference_length_( reference_length ),
      data_length_( data_length ),
      reference_( max( reference_length * 2, min( data_length, size_t( max_chunk_size ) ) ) ),
      reference_fft_( reference_.size() ),
      data_( reference_.size() ),
      data_fft_( reference_.size() ),
      reference_plan_( reference_, reference_fft_, FFTW_FORWARD, FFTW_ESTIMATE ),
      data_plan_( data_, data_fft_, FFTW_FORWARD, FFTW_ESTIMATE ),
      inverse_plan_( data_fft_, data_, FFTW_BACKWARD, FFTW_ESTIMATE )
{
    if ( reference_length < 1 ) {
        throw runtime_error( "invalid reference_length" );
    }

    if ( data_length < 1 ) {
        throw runtime_error( "invalid data_length" );
    }

    if ( reference_length > data_length ) {
        throw runtime_error( "reference_length cannot be longer than data_length" );
    }

    if ( reference_.size() < data_.size() ) {
        if ( (reference_.size() - reference_length_) % 64 ) {
            throw runtime_error( "max_chunk_size leads to unaligned access" );
        }
    }
}

void CrossCorrelator::correlate_fast( const Signal & reference, const Signal & data,
                                      vector<float> & output )
{
    if ( reference.size() != reference_length_ ) {
        throw runtime_error( "invalid reference length" );
    }

    if ( data.size() != data_length_ ) {
        throw runtime_error( "invalid data length" );
    }

    if ( output.size() != data_length_ - reference_length_ ) {
        throw runtime_error( "invalid output length (must be data_length - reference_length)" );
    }

    const unsigned int interval = reference_.size() - reference_length_;

    memcpy( reference_.data(), reference.data(), reference_length_ * sizeof( complex<float> ) );
    reference_plan_.execute();

    float reference_power = 0;
    for ( unsigned int i = 0; i < reference_fft_.size(); i++ ) {
        reference_power += norm( reference_fft_[ i ] );
    }

    for ( unsigned int offset = 0; offset < data.size(); offset += interval ) {
        cerr << "\r" << int( 100.0 * offset / float( data.size() ) ) << "%             ";

        fill( data_.begin(), data_.end(), 0 );
        memcpy( data_.data(), data.data() + offset, min( data_.size(),
                                                         data.size() - offset ) * sizeof( complex<float> ) );
        data_plan_.execute();

        /* multiply data_fft_ in place by conjugate of reference */
        for ( unsigned int i = 0; i < data_fft_.size(); i++ ) {
            data_fft_[ i ] *= conj( reference_fft_[ i ] );
        }

        /* inverse FFT */
        inverse_plan_.execute();

        for ( unsigned int lag = 0; lag < interval and lag + offset < output.size(); lag++ ) {
            output[ offset + lag ] = abs( data_[ lag ] ) / reference_power;
        }
    }

    cerr << "\r                                    \n";
}
