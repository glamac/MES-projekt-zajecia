#include <cmath>
#include <fstream>
#include <ios>
#include <iostream>

#include "MESCore.hpp"
#include "MESTypes.hpp"
#include "MESMatrix.hpp"
#include <getopt.h>
#include <ostream>

void printUsage() {
	std::cout << "Usage: messim [-ip 2|3|4] [-o <outFileName>] <gridFile>";
	exit(-1);
}

int main(int argc, char* argv[]) {
    MES::grid SimulationGrid;

    if (argc <= 0) {
   		printUsage();
    }
    bool fileIO;
    std::string fileIOname;

    for(int i = 1; i < argc; i++) {
    	std::string arg = argv[i];
     	if(arg.starts_with("-ip") || arg.starts_with("--integration-points")) {
      		int integration_points = std::stoi(argv[++i]);
       		SimulationGrid.simData.integrationPoints_Boundary = integration_points;
       		SimulationGrid.simData.integrationPoints_Surface = integration_points;
    	}
      	if(arg.starts_with("-o") || arg.starts_with("--output")) {
       		fileIO = true;
        	fileIOname = argv[++i];
       	}
      	if(!arg.starts_with("-")) {
     		std::cout << "Running simulation; file name '" << arg << "', " << SimulationGrid.simData.integrationPoints_Surface << " Integration Points\n";
     		MES::parseTextFile(SimulationGrid, arg);
       		goto filenameFound;
    	}
    }
    // wrong usage; exit
    printUsage();
    filenameFound: // right usage; continue

    const auto& GLOB = SimulationGrid.simData;
    const auto stepTime = GLOB.SimulationStepTime;

    MES::Matrix t0(GLOB.numNodes, 1, GLOB.InitialTemp);

    std::ofstream fileOut;
    std::ostream* output = &std::cout;
    if(fileIO) {
   		fileOut.open(fileIOname.c_str());
    	output = &fileOut;
    }

    for(MES::real_t time = GLOB.SimulationStepTime; time <= GLOB.SimulationTime; time += stepTime) {

  		SimulationGrid.clearGlobalMatrices();
	    for(auto& elem : SimulationGrid.Elements) {
	    	elem
	     		.calculateMatrices()
	     		.calculateBoundaryCondition()
	     		.aggregateToGlobal();
	    }

        auto Hg_Cg_dt = (GLOB.H + GLOB.C / stepTime);
        auto Pg_Cg_dt_t0 = -1 * (-(GLOB.C / stepTime) * t0 - GLOB.vP );
        auto t1 = MES::gauss(Hg_Cg_dt, Pg_Cg_dt_t0);
        // auto t1 = MES::gauss(a, b);
        (*output) << "t: " << time << "s\t";
        (*output) << "min: " << t1.min() << " max: " << t1.max() << '\n';
        // std::cout << "[H]+[C]/d: " << a << '\n';
        // std::cout << "{P}+{[C]/dT}*{T0}: " << b.transpose() << '\n';
        // std::cout << " temperatures: " << t1.transpose() << std::endl;
        t0 = t1;
    }

    if(fileIO) {
    	fileOut.close();
    }


}
