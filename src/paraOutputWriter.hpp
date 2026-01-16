#ifndef _PARA_OUTPUT_WRITER_H
#define _PARA_OUTPUT_WRITER_H
#include <fstream>
#include "MESTypes.hpp"
#include "MESMatrix.hpp"

namespace MES {
	class paraOutputWriter {
		std::string simulationName;
		std::string outputPath;
		int steps = 0;

		std::ofstream outPVDFile;
	public:
		void initialize(std::string name);
		void writeTimeStep(MES::grid& grid, MES::Matrix temperatures, MES::real_t time);
		void close();

	};
}

#endif // _PARA_OUTPUT_WRITER_H
