#include "MESTypes.hpp"
#include <iostream>

int main() {
    MES::grid A = MES::parseTextFile("Test3_31_31_kwadrat.txt");
    std::cout << "After Grid Get: " << A.numNodes << " " << A.numElements << "\n";
    for(int i = 0; i < A.numNodes; i++) {
        std::cout << i << ": " << A.Node[i].x << " " << A.Node[i].y << std::endl;
    }
    for(int i = 0; i < A.numElements; i++) {
        std::cout << i << ": " << A.Element[i].ID[0] << " " << A.Element[i].ID[1] << " "
            << A.Element[i].ID[2] << " " << A.Element[i].ID[3] << std::endl;
    }
    std::cout << "----";
}
