#ifndef _MES_CORE_H
#define _MES_CORE_H
#include "MESMatrix.hpp"
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
		MES::Matrix H = MES::Matrix(4, 4);
		MES::Matrix C = MES::Matrix(4, 4);
		bool operator==(MES::element rhs) {
            for(int i = 0; i < ID.size(); i++) {
                if (rhs.ID[i] != ID[i]) return false;
            }
            return true;
        }
        element& calculateMatrices();
        // H_BC before summing into H, for debugging purposes.
        // Matrix H has this added to it after calculateBoundaryCondution().
        MES::Matrix H_BC = MES::Matrix(4, 4);
        MES::Matrix vP = MES::Matrix(4, 1);
        element& calculateBoundaryCondition();
        element& aggregateToGlobal();
	};
}



#endif // _MES_CORE_H
