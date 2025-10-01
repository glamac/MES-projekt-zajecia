#include <string>
#include <fstream>
#include "MESTypes.hpp"

bool MES::parseTextFile(grid *Grid, std::string Filename) {
	std::ifstream DataFile(Filename);
	if (!DataFile.is_open()) {
		std::cout << "Error opening file " << Filename << "\n";
		return false;
	}

	std::string currline

	while(DataFile.get

}

void MES::parserHandleLine(grid *Grid, std::string Line);
