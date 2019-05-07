AM_CPPFLAGS = $(CXX17_FLAGS) $(fftw3f_CFLAGS)
AM_CXXFLAGS = $(PICKY_CXXFLAGS)

LDFLAGS = -pthread
LDADD = $(fftw3f_LIBS)

SHARED = helpers.hh cross-correlation.hh cross-correlation.cc

bin_PROGRAMS = power parse xcorr-slow xcorr-fast

power_SOURCES = power.cc $(SHARED)

parse_SOURCES = parse.cc $(SHARED)

xcorr_slow_SOURCES = xcorr-slow.cc $(SHARED)

xcorr_fast_SOURCES = xcorr-fast.cc $(SHARED)

check_PROGRAMS = xcorr-correctness-test

xcorr_correctness_test_SOURCES = xcorr-correctness-test.cc $(SHARED)

TESTS = $(check_PROGRAMS)
