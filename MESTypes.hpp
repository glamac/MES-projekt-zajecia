#include <string>
namespace MES {
	typedef double real_t;
	struct node {
		real_t x;
		real_t y;
	};
	struct element {
		int ID[4];
	};
	struct grid {
		const int numNodes;
		const int numElements;
		element *Element;
		node *Node;
		grid(const int nN, const int nE) : numNodes(nN), numElements(nE) {
			Element = new element[nE];
			Node = new node[nN];
		};
	};
	
	bool parseTextFile(grid *Grid, std::string Filename);
	void parserHandeLine(grid *Grid, std::string Line);

	struct globalData {
		real_t SimulationTime;
		real_t SimulationStepTime;
		real_t Conductivity;
		real_t Alpha;
		real_t Tot;
		real_t InitialTemp;
		real_t Density;
		real_t SpecificHeat;
		const int numNodes;
		const int numElements;
	};

}
