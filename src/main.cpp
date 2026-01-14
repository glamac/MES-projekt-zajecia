#include <cmath>
#include <iostream>

#include "MESCore.hpp"
#include "MESTypes.hpp"
#include "MESMatrix.hpp"
#include <getopt.h>

int main(int argc, char* argv[]) {
    MES::grid A;
    //MES::parseTextFile(A, "Test2_4_4_Mixgrid.txt");
    MES::parseTextFile(A, "Test1_4_4.txt");
    printf("numElements: %u", GLOB.numElements);

    for(auto& elem : A.Element) {
    	elem
     		.calculateMatrices()
     		.calculateBoundaryCondition()
     		.aggregateToGlobal();
    }

    auto tauInf = MES::gauss(GLOB.H, GLOB.vP).transpose();

    MES::Matrix t0(16, 1, GLOB.InitialTemp);
    const auto stepTime = GLOB.SimulationStepTime;
    for(MES::real_t deltaTau = 50; deltaTau <= GLOB.SimulationTime; deltaTau += GLOB.SimulationStepTime) {
        auto a = (GLOB.H + GLOB.C / stepTime);
        auto b = -1 * (-(GLOB.C / stepTime) * t0 - GLOB.vP );
        auto Hg_Cg_dt = GLOB.H + (GLOB.C / stepTime);
        auto Pg_Cg_dt_t0 = GLOB.vP + (GLOB.C/ stepTime)*t0;
        //auto t1 = MES::gauss(Hg_Cg_dt, -1*Pg_Cg_dt_t0);
        auto t1 = MES::gauss(a, b);
        std::cout << "Time: " << deltaTau << std::endl;
        std::cout << "min: " << t1.min() << " max: " << t1.max() << std::endl;
        // std::cout << "[H]+[C]/d: " << a << '\n';
        // std::cout << "{P}+{[C]/dT}*{T0}: " << b.transpose() << '\n';
        std::cout << " temperatures: " << t1.transpose() << std::endl;
        t0 = t1;
        break;
    }




}
