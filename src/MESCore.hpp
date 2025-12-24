#ifndef _MES_CORE_H
#define _MES_CORE_H
#include "MESMatrix.hpp"
#include <memory>
#define INTEGRATION_POINTS 2
#include <array>


namespace MES {
	typedef double real_t;

	struct node {
		union {
			struct {real_t x, y;};
			struct {real_t ksi, eta;};
		};
		bool BC = false;
	};
	struct element {
		std::array<int, 4> ID;
		std::array<MES::node*, 4> Node;
		MES::Matrix<4, 4> H;
		MES::Matrix<4, 4> C;
		bool operator==(MES::element rhs) {
            for(int i = 0; i < ID.size(); i++) {
                if (rhs.ID[i] != ID[i]) return false;
            }
            return true;
        }
	};
}



#endif // _MES_CORE_H
