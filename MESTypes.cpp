#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include "MESTypes.hpp"

MES::globalData GLOB;

void MES::parseTextFile(MES::grid& SimulationGrid, std::string Filename) {
	std::ifstream DataFile(Filename);
	if (!DataFile.is_open()) {
		std::cout << "Error opening file " << Filename << "\n";
		return;
	}

	std::string currline, junk;

	while(std::getline(DataFile, currline)){

        std::istringstream lineContent(currline);
        if(currline.starts_with("SimulationTime ")){
            lineContent >> junk >> GLOB.SimulationTime;
        }
        if(currline.starts_with("SimulationStepTime ")){
            lineContent >> junk >> GLOB.SimulationStepTime;
        }
        if(currline.starts_with("Conductivity ")){
            lineContent >> junk >> GLOB.Conductivity;
        }
        if(currline.starts_with("Alfa ")){
            lineContent >> junk >> GLOB.Alpha;
        }
        if(currline.starts_with("Tot ")){
            lineContent >> junk >> GLOB.Tot;
        }
        if(currline.starts_with("InitialTemp ")){
            lineContent >> junk >> GLOB.InitialTemp;
        }
        if(currline.starts_with("Density ")){
            lineContent >> junk >> GLOB.Density;
        }
        if(currline.starts_with("SpecificHeat ")){
            lineContent >> junk >> GLOB.SpecificHeat;
        }
        if(currline.starts_with("Nodes number ")){
            lineContent >> junk >> junk >> SimulationGrid.numNodes;
            GLOB.numNodes = SimulationGrid.numNodes;
        }
        if(currline.starts_with("Elements number ")){
            lineContent >> junk >> junk >> SimulationGrid.numElements;
            GLOB.numElements = SimulationGrid.numElements;
        }
        if(currline.starts_with("*Node")){
            int nodeIndex;
            SimulationGrid.Node.resize(SimulationGrid.numNodes);
            real_t x, y;
            for(int i = 0; i < SimulationGrid.numNodes; i++) {
                std::getline(DataFile, currline);
                std::istringstream nodeLine(currline);
                nodeLine >> nodeIndex >> junk;
                nodeLine >> x >> junk;
                nodeLine >> y;
                SimulationGrid.Node[i] = {x, y};
            }
        }
        if(currline.starts_with("*Element")){
            int elemIndex;
            SimulationGrid.Element.resize(SimulationGrid.numElements);
            for(int i = 0; i < SimulationGrid.numElements; i++) {
                std::getline(DataFile, currline);
                std::istringstream elemLine(currline);
                auto& currElem = SimulationGrid.Element[i];
                elemLine >> elemIndex >> junk;
                elemLine >> currElem.ID[0] >> junk;
                elemLine >> currElem.ID[1] >> junk;
                elemLine >> currElem.ID[2] >> junk;
                elemLine >> currElem.ID[3];
            }
            for(int i = 0; i < SimulationGrid.numElements; i++) {
                for(int j = 0; j < 4; j++) {
                    auto& currElem = SimulationGrid.Element[i];
                    currElem.ID[j] -= 1;
                    currElem.Node[j] = &SimulationGrid.Node[currElem.ID[j]];
                }
            }
        }
        if(currline.starts_with("*BC")){
            std::getline(DataFile, currline);
            int nodeID;
            std::istringstream bcLine(currline);
            while(bcLine) {
                bcLine >> nodeID;
                SimulationGrid.Node[nodeID-1].BC = true;
                bcLine >> junk;
            }
        }

    }
}

