AM_CPPFLAGS = $(CXX17_FLAGS) $(fftw3f_CFLAGS)
AM_CXXFLAGS = $(PICKY_CXXFLAGS)

SHARED = helpers.hh

bin_PROGRAMS = power parse crosscorrelation fftxcorr

power_SOURCES = power.cc $(SHARED)

parse_SOURCES = parse.cc $(SHARED)

crosscorrelation_SOURCES = crosscorrelation.cc $(SHARED)

fftxcorr_SOURCES = fftxcorr.cc $(SHARED)
fftxcorr_LDFLAGS = -pthread
fftxcorr_LDADD = $(fftw3f_LIBS)
