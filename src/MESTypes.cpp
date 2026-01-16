#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include "MESTypes.hpp"

void MES::parseTextFile(MES::grid& SimulationGrid, std::string Filename) {
	std::ifstream DataFile(Filename);
	if (!DataFile.is_open()) {
		std::cout << "Error opening file " << Filename << "\n";
		return;
	}

	simulationSharedData& simData = SimulationGrid.simData;

	std::string currline, junk;
	// parsing logic
	while(std::getline(DataFile, currline)){
        std::istringstream lineContent(currline);
        if(currline.starts_with("SimulationTime ")){
            lineContent >> junk >> simData.SimulationTime;
        }
        if(currline.starts_with("SimulationStepTime ")){
            lineContent >> junk >> simData.SimulationStepTime;
        }
        if(currline.starts_with("Conductivity ")){
            lineContent >> junk >> simData.Conductivity;
        }
        if(currline.starts_with("Alfa ")){
            lineContent >> junk >> simData.Alpha;
        }
        if(currline.starts_with("Tot ")){
            lineContent >> junk >> simData.Tot;
        }
        if(currline.starts_with("InitialTemp ")){
            lineContent >> junk >> simData.InitialTemp;
        }
        if(currline.starts_with("Density ")){
            lineContent >> junk >> simData.Density;
        }
        if(currline.starts_with("SpecificHeat ")){
            lineContent >> junk >> simData.SpecificHeat;
        }
        if(currline.starts_with("Nodes number ")){
            lineContent >> junk >> junk >> SimulationGrid.numNodes;
            simData.numNodes = SimulationGrid.numNodes;
        }
        if(currline.starts_with("Elements number ")){
            lineContent >> junk >> junk >> SimulationGrid.numElements;
            simData.numElements = SimulationGrid.numElements;
        }
        if(currline.starts_with("*Node")){
            int nodeIndex;
            SimulationGrid.Nodes.resize(SimulationGrid.numNodes);
            real_t x, y;
            for(int i = 0; i < SimulationGrid.numNodes; i++) {
                std::getline(DataFile, currline);
                std::istringstream nodeLine(currline);
                nodeLine >> nodeIndex >> junk;
                nodeLine >> x >> junk;
                nodeLine >> y;
                SimulationGrid.Nodes[i] = {x, y};
            }
        }
        if(currline.starts_with("*Element")){
            int elemIndex;
            SimulationGrid.Elements.reserve(SimulationGrid.numElements);
            for(int i = 0; i < SimulationGrid.numElements; i++) {
                std::getline(DataFile, currline);
                std::istringstream elemLine(currline);
                auto& currElem = SimulationGrid.Elements.emplace_back(SimulationGrid.simData);
                elemLine >> elemIndex >> junk;
                elemLine >> currElem.ID[0] >> junk;
                elemLine >> currElem.ID[1] >> junk;
                elemLine >> currElem.ID[2] >> junk;
                elemLine >> currElem.ID[3];
            }
            for(int i = 0; i < SimulationGrid.numElements; i++) {
                for(int j = 0; j < 4; j++) {
                    auto& currElem = SimulationGrid.Elements[i];
                    currElem.ID[j] -= 1;
                    currElem.Node[j] = &SimulationGrid.Nodes[currElem.ID[j]];
                }
            }
        }
        if(currline.starts_with("*BC")){
            std::getline(DataFile, currline);
            int nodeID;
            std::istringstream bcLine(currline);
            while(bcLine) {
                bcLine >> nodeID;
                SimulationGrid.Nodes[nodeID-1].BC = true;
                bcLine >> junk;
            }
        }
    }

	const auto globMatSize = simData.numNodes;
	simData.H = Matrix(globMatSize, globMatSize);
	simData.C = Matrix(globMatSize, globMatSize);
	simData.vP = Matrix(globMatSize, 1);
}
