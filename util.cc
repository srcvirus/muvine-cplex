#include "util.h"

#include <stdarg.h>
#include <stdio.h>
#include <algorithm>
#include <map>
#include <string>

const int kInfinity = 100000000;

void PrintDebugMessage(const char* location, const char* fmt_string, ...) {
  va_list args;
  va_start(args, fmt_string);
  std::string str = location;
  str += fmt_string;
  vprintf(str.c_str(), args);
  fflush(stdout);
  va_end(args);
}

template <class T>
double GetMean(const std::vector<T>& data) {
  T sum = T(0);
  const size_t kNumElements = data.size();
  for (int i = 0; i < data.size(); ++i) {
    sum += data[i];
  }
  return sum / static_cast<T>(kNumElements);
}

template <class T>
T GetNthPercentile(const std::vector<T>& data, int n) {
  std::vector<T> temp_data_buffer = data;
  sort(temp_data_buffer.begin(), temp_data_buffer.end());
  const size_t kNumElements = data.size();
  int rank = n * kNumElements;
  if (rank % 100) {
    rank = (rank / 100) + 1;
  } else
    rank /= 100;
  --rank;
  return temp_data_buffer[rank];
}

template <class T>
std::vector<std::pair<T, double> > GetCDF(const std::vector<T>& data) {
  int precision = 1;
  std::vector<T> temp_data_buffer = data;
  if (typeid(temp_data_buffer[0]) == typeid(double) ||
      typeid(temp_data_buffer[0]) == typeid(float)) {
    precision = 1000;
  }
  std::map<int, int> cdf;
  for (int i = 0; i < temp_data_buffer.size(); ++i) {
    int bucket_index = temp_data_buffer[i] * precision;
    if (cdf[bucket_index])
      cdf[bucket_index]++;
    else
      cdf[bucket_index] = 1;
  }
  std::map<int, int>::iterator prev = cdf.begin(), current = cdf.begin();
  current++;
  for (; current != cdf.end(); current++, prev++) {
    current->second += prev->second;
  }
  int total = temp_data_buffer.size();
  std::vector<std::pair<T, double> > ret;
  for (current = cdf.begin(); current != cdf.end(); ++current) {
    T first = static_cast<T>(current->first) / static_cast<T>(precision);
    double second =
        static_cast<double>(current->second) / static_cast<double>(total);
    ret.push_back(std::make_pair(first, second));
  }
  return ret;
}
