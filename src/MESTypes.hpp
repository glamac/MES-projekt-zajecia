#ifndef _MES_TYPES_H
#define _MES_TYPES_H
#include <string>
#include <vector>
#include "MESCore.hpp"
#include "MESData.h"

namespace MES {

	struct grid {
		int numNodes;
		int numElements;
        std::vector<element> Elements;
        std::vector<node> Nodes;
        void setNodes(int nN) {numNodes = nN;};
        void setElems(int nE) {numElements = nE;};
        simulationSharedData simData;
        void clearGlobalMatrices() {
	       	simData.H = {0};
	        simData.C = {0};
	        simData.vP = {0};
        }
	};

	void parseTextFile(grid& SimulationGrid, std::string Filename);
}
#endif // _MES_TYPES_H
