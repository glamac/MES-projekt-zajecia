#ifndef _MES_TEST_H
#define _MES_TEST_H
#include "MESMatrix.hpp"

namespace MES::Test {
    inline void testMatrix() {
        constexpr Matrix<2, 3> m1 = {
            2, 2, 2,
            2, 2, 2
        };

        constexpr Matrix<3, 2> m2 = {
            3, 3,
            3, 3,
            3, 3
        };

        constexpr Matrix<3, 3> expected = {
            12, 12, 12,
            12, 12, 12,
            12, 12, 12,
        };

        constexpr auto m3 = m2 * m1;
        static_assert(m3 == expected, "Either matrix multiplication broke or the initializer list thing did");
    }
}
#endif // _MES_TEST_H