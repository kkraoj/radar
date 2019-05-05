AM_CPPFLAGS = $(CXX17_FLAGS) $(fftw3f_CFLAGS)
AM_CXXFLAGS = $(PICKY_CXXFLAGS)

SHARED = helpers.hh cross-correlation.hh cross-correlation.cc

bin_PROGRAMS = power parse xcorr-slow xcorr-fast

power_SOURCES = power.cc $(SHARED)

parse_SOURCES = parse.cc $(SHARED)

xcorr_slow_SOURCES = xcorr-slow.cc $(SHARED)

xcorr_fast_SOURCES = xcorr-fast.cc $(SHARED)
xcorr_fast_LDFLAGS = -pthread
xcorr_fast_LDADD = $(fftw3f_LIBS)
