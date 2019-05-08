#include "helpers.hh"

#include <fftw3.h> /* must be included last, after <complex> */

void correlate_slow( const Signal & reference, const Signal & data, std::vector<float> & output );

class FFTPlan
{
    Signal *input_, *output_;
    
    fftwf_plan plan_;

public:
    FFTPlan( Signal & input, Signal & output,
             const int sign, const int flags );
    ~FFTPlan();
    void execute();

    FFTPlan( const FFTPlan & other ) = delete;
    FFTPlan & operator=( const FFTPlan & other ) = delete;
};

class CrossCorrelator
{
    size_t reference_length_, data_length_;

    Signal reference_, reference_fft_;
    Signal data_, data_fft_;

    FFTPlan reference_plan_, data_plan_, inverse_plan_;

public:
    CrossCorrelator( const size_t reference_length,
                     const size_t data_length,
                     const size_t max_chunk_size );

    void correlate_fast( const Signal & reference, const Signal & data, std::vector<float> & output );
};
