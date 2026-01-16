#ifndef _MES_DATA_H
#define _MES_DATA_H
#include "MESMatrix.hpp"
namespace MES {
	typedef double real_t;
	struct simulationSharedData {
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
#endif // _MES_DATA_H
