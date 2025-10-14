#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include "MESTypes.hpp"

MES::globalData GLOB;

MES::grid MES::parseTextFile(std::string Filename) {
	std::ifstream DataFile(Filename);
	if (!DataFile.is_open()) {
		std::cout << "Error opening file " << Filename << "\n";
		return grid();
	}

	std::string currline, junk;

    grid SimulationGrid;

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
        }
        if(currline.starts_with("Elements number ")){
            lineContent >> junk >> junk >> SimulationGrid.numElements;
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
                elemLine >> elemIndex >> junk;
                elemLine >> SimulationGrid.Element[i].ID[0] >> junk;
                elemLine >> SimulationGrid.Element[i].ID[1] >> junk;
                elemLine >> SimulationGrid.Element[i].ID[2] >> junk;
                elemLine >> SimulationGrid.Element[i].ID[3];
            }
        }

    }
    return SimulationGrid;
}

