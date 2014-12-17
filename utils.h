#ifndef kuramoto_utils_h
#define kuramoto_utils_h

#include <iostream>
#include <iterator>

#define PO(e) std::cout << #e << ": " << std::setprecision(10) << e << std::endl
#define PIO(v) std::cout << #v << ": "; copy(v.begin(), v.end(), ostream_iterator<double>(cout, " ")); std::cout << " (size: " << v.size() << ")" << std::endl

#endif
