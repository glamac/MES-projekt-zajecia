#include <string>
#include <vector>
#include "MESCore.hpp"

namespace MES {

	struct grid {
		int numNodes;
		int numElements;
        std::vector<element> Element;
        std::vector<node> Node;
		grid(const int nN, const int nE) : numNodes(nN), numElements(nE) {
            Element.resize(nE);
            Node.resize(nN);
		};
        grid() : grid(0,0) {};
        void setNodes(int nN) {numNodes = nN;};
        void setElems(int nE) {numElements = nE;};
		
	};
	
	void parseTextFile(grid& SimulationGrid, std::string Filename);

	struct globalData {
		real_t SimulationTime;
		real_t SimulationStepTime;
		real_t Conductivity;
		real_t Alpha;
		real_t Tot;
		real_t InitialTemp;
		real_t Density;
		real_t SpecificHeat;
		int numNodes;
		int numElements;
	};
}

extern MES::globalData GLOB;

