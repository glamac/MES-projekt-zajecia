#include "MESTypes.hpp"
#include "MESCalculus.hpp"
#include <iostream>

MES::real_t testfunc(MES::real_t x) {
    return 5*x*x + 3*x + 6;
}

MES::real_t testfunc2D(MES::real_t x, MES::real_t y) {
    return 5*x*x*y*y + 3*x*y + 6;
}

int main() {
    MES::grid A = MES::parseTextFile("Test1_4_4.txt");
    std::cout << "After Grid Get: " << A.numNodes << " " << A.numElements << "\n";
    for(int i = 0; i < A.numNodes; i++) {
        std::cout << i << ": " << A.Node[i].x << " " << A.Node[i].y << std::endl;
    }
    for(int i = 0; i < A.numElements; i++) {
        std::cout << i << ": " << A.Element[i].ID[0] << " " << A.Element[i].ID[1] << " "
            << A.Element[i].ID[2] << " " << A.Element[i].ID[3] << std::endl;
    }
    std::cout << "----";
    std::cout << "Integrating 2D func" << MES::CalculusHelper::integrate(testfunc2D, 3);
}
