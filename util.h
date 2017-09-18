#ifndef UTIL_H_
#define UTIL_H_

#include "datastructure.h"

#include <utility>
#include <vector>

#define ONE_GIG 1000000000ULL
#define EPS 1e-8

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__) " "

#ifdef DBG
#define DEBUG(...) PrintDebugMessage(AT, __VA_ARGS__)
#else
#define DEBUG(...)
#endif

extern const int kInfinity;

// Prints a message formatted using fmt_string. Prepends the location of the
// file before the message.
void PrintDebugMessage(const char* location, const char* fmt_string, ...);

// Returns the mean of the vector of data.
template <typename T>
double GetMean(const std::vector<T>& data);

// Returns the Nth percentile of the vector of data. Nth percentile is
// calculated as the (ceil(N / 100) * data.size() - 1)-th element from an array
// obtained by sorting data.
template <typename T>
T GetNthPercentile(const std::vector<T>& data, int n);

// Returns the Cumalitive Distribution Frequence of the data items sotred in
// data. If the data items sotred in data are of type double, then a precision
// of 3 digits after the decimal points is used.
template <typename T>
std::vector<std::pair<T, double> > GetCDF(const std::vector<T>& data);

#endif  // UTIL_H_
