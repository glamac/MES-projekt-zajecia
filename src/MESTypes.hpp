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
	void globMatrixInit();

	struct globalData {
		// Time to run the simulation.
		real_t SimulationTime;
		// Time to step the simulation by.
		real_t SimulationStepTime;
		// Conductivity of material.
		real_t Conductivity;
		real_t Alpha;
		// Temperature on the boundaries.
		real_t Tot;
		// Initial temperature.
		real_t InitialTemp;
		real_t Density;
		real_t SpecificHeat;
		int numNodes;
		int numElements;

		// Global Matrices; initialized
		// when file is read.

		// Global H Matrix
		Matrix H;
		// Global C Matrix
		Matrix C;
		// Global P Vector
		Matrix vP;

		// Number of integration points to use when
		// calculating H and C matrices.
		unsigned int integrationPoints_Surface = 2;
		// Number of integration points to use when
		// calculating the boundary condition math.
		unsigned int integrationPoints_Boundary = 2;

	};
}

extern MES::globalData GLOB;
