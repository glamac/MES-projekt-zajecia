#include "paraOutputWriter.hpp"
#include "MESData.h"
#include <filesystem> // C++17
#include <fstream>
#include <ostream>
#include <string>

namespace MES {

namespace fs = std::filesystem;
	void paraOutputWriter::initialize(std::string name) {
		simulationName = name;
		std::error_code ec;
		if (!fs::create_directory(simulationName, ec) && !fs::exists(simulationName)) {
        throw std::runtime_error("Failed to create directory: " +
                                 simulationName + " - " + ec.message());
    }
		outputPath = simulationName + "/";
    std::string fullFilePath = fs::current_path() / outputPath / (simulationName + ".pvd");

    outPVDFile.open(fullFilePath);

    if(!outPVDFile.good()) {
        throw std::runtime_error("Error: can't write to file: " + fullFilePath);
    }
		outPVDFile << "<?xml version=\"1.0\"?>\n"
		<< "<VTKFile type=\"Collection\" version=\"0.1\">\n"
		<< "<Collection>\n";
		std::cout << "writing!\n";
	}
	void paraOutputWriter::writeTimeStep(MES::grid& simGrid, MES::Matrix temps, MES::real_t time) {
		std::string fileName = simulationName + "_" + std::to_string(steps++) + ".vtu";
		std::ofstream vtuFile(
			outputPath + fileName
		);
		if(!vtuFile.is_open() || !outPVDFile.is_open()) throw ("Error: can't write to file.");
		outPVDFile << "<DataSet timestep=\"" << time << "\" file=\"" << fileName << "\"/>\n";
		vtuFile << "<?xml version=\"1.0\"?>\n"
		<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n"
		<< "<UnstructuredGrid>\n"
		<< "<Piece NumberOfPoints=\"" << simGrid.numNodes
		<< "\" NumberOfCells=\"" << simGrid.numElements << "\">\n"
		<< "<Points>\n<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		for(auto& node : simGrid.Nodes) {
			vtuFile << node.x << " " << node.y << " 0.0\n";
		}
		vtuFile << "</DataArray>\n</Points>\n"
		<< "<Cells>\n"
		<< "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
		for(auto& el : simGrid.Elements) {
			vtuFile << el.ID[1] << " " << el.ID[0] << " " << el.ID[3] << " " << el.ID[2] << '\n';
		}
		vtuFile << "</DataArray>\n"
		<< "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
		for(int i = 0; i < simGrid.numElements; i++) {
			vtuFile << 4*(i+1) << " ";
		}
		vtuFile << "\n</DataArray>\n"
		<< "<DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
		for(int i = 0; i < simGrid.numElements; i++) {
			vtuFile << "9" << " "; // 9 to quad, chyba.
		}
		vtuFile << "\n</DataArray>\n"
		<< "</Cells>\n"
		<< "<PointData Scalars=\"Temperature\">\n"
		<< "<DataArray type=\"Float64\" Name=\"Temperature\" NumberOfComponents=\"1\" format=\"ascii\">\n";
		for(int i = 0; i < simGrid.numNodes; i++) {
			vtuFile << temps[i, 0] << " ";
		}
		vtuFile << "\n</DataArray>\n</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
		vtuFile.close();
	}

	void paraOutputWriter::close() {
		if(!outPVDFile.is_open()) throw ("Error: can't write to file.");
		outPVDFile << "</Collection>\n</VTKFile>\n";
		outPVDFile.close();
	}
}
