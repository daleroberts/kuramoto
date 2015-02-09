#pragma once
#include <algorithm>

namespace sa
{
    template <typename Container>
    void sort(Container &c) {
        std::sort(c.begin(), c.end());
    }
}
