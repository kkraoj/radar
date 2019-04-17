AM_CPPFLAGS = $(CXX17_FLAGS)
AM_CXXFLAGS = $(PICKY_CXXFLAGS)

SHARED = helpers.hh

bin_PROGRAMS = power parse crosscorrelation

power_SOURCES = power.cc $(SHARED)
parse_SOURCES = parse.cc $(SHARED)
crosscorrelation_SOURCES = crosscorrelation.cc $(SHARED)
