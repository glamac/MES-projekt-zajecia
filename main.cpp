#include "MESCore.hpp"
#include "MESTypes.hpp"
#include "MESCalculus.hpp"
#include <iostream>

MES::real_t testfunc(MES::real_t x) {
    return 5*x*x + 3*x + 6;
}

MES::real_t testfunc2D(MES::real_t x, MES::real_t y) {
    return 5*x*x*y*y + 3*x*y + 6;
}
void jacobian(MES::grid A, MES::node p1) {
    using MSF = MES::ShapeFunctions;
    printf("p1: %.8lf, %.8lf\n", p1.x, p1.y);
    for(int i = 0; i < 4; i++) {
        auto myNode = A.Node[A.Element[0].ID[i]-1];
        printf("n%i: %.8lf, %.8lf\n", i, myNode.x, myNode.y);
    }
    MES::real_t J[2][2] = {0};
    printf("%.8lf, %.8lf\n",
            MSF::dKsi::N1(p1.eta) * A.Node[A.Element[1].ID[0]].x +
            MSF::dKsi::N2(p1.eta) * A.Node[A.Element[1].ID[1]].x +
            MSF::dKsi::N3(p1.eta) * A.Node[A.Element[1].ID[2]].x +
            MSF::dKsi::N4(p1.eta) * A.Node[A.Element[1].ID[3]].x,

            MSF::dKsi::N1(p1.eta) * A.Node[A.Element[1].ID[0]].y +
            MSF::dKsi::N2(p1.eta) * A.Node[A.Element[1].ID[1]].y +
            MSF::dKsi::N3(p1.eta) * A.Node[A.Element[1].ID[2]].y +
            MSF::dKsi::N4(p1.eta) * A.Node[A.Element[1].ID[3]].y
        );
    printf("%.8lf, %.8lf\n",
            MSF::dEta::N1(p1.ksi) * A.Node[A.Element[1].ID[0]].x +
            MSF::dEta::N2(p1.ksi) * A.Node[A.Element[1].ID[1]].x +
            MSF::dEta::N3(p1.ksi) * A.Node[A.Element[1].ID[2]].x +
            MSF::dEta::N4(p1.ksi) * A.Node[A.Element[1].ID[3]].x,

            MSF::dEta::N1(p1.ksi) * A.Node[A.Element[1].ID[0]].y +
            MSF::dEta::N2(p1.ksi) * A.Node[A.Element[1].ID[1]].y +
            MSF::dEta::N3(p1.ksi) * A.Node[A.Element[1].ID[2]].y +
            MSF::dEta::N4(p1.ksi) * A.Node[A.Element[1].ID[3]].y
        );
}

int main() {
    MES::grid A = MES::parseTextFile("Test2_4_4.txt");
    std::cout << "After Grid Get: " << A.numNodes << " " << A.numElements << "\n";
    for(int i = 0; i < A.numNodes; i++) {
        std::cout << i << ": " << A.Node[i].x << " " << A.Node[i].y << std::endl;
    }
    for(int i = 0; i < A.numElements; i++) {
        std::cout << i << ": " << A.Element[i].ID[0] << " " << A.Element[i].ID[1] << " "
            << A.Element[i].ID[2] << " " << A.Element[i].ID[3] << std::endl;
    }
    std::cout << "----\n";
    for(int i = 0; i < 4; i++) {
        auto intP = MES::getIntegrationPoint(i);
        printf("Integration point %i: %lf, %lf\n", i, intP.x, intP.y);
    }
    std::cout << "---- dN/dKsi\n";
    for(int i = 0; i < INTEGRATION_POINTS*INTEGRATION_POINTS; i++) {
        printf("pc%i:  %.5lf  %.5lf  %.5lf  %.5lf\n",
            i, 
            MES::ShapeFunctions::dKsi::N1(MES::getIntegrationPoint(i).eta),
            MES::ShapeFunctions::dKsi::N2(MES::getIntegrationPoint(i).eta),
            MES::ShapeFunctions::dKsi::N3(MES::getIntegrationPoint(i).eta),
            MES::ShapeFunctions::dKsi::N4(MES::getIntegrationPoint(i).eta)
        );
    }
    std::cout << "---- dN/dEta\n";
    for(int i = 0; i < INTEGRATION_POINTS*INTEGRATION_POINTS; i++) {
        printf("pc%i:  %.5lf  %.5lf  %.5lf  %.5lf\n",
            i, 
            MES::ShapeFunctions::dEta::N1(MES::getIntegrationPoint(i).ksi),
            MES::ShapeFunctions::dEta::N2(MES::getIntegrationPoint(i).ksi),
            MES::ShapeFunctions::dEta::N3(MES::getIntegrationPoint(i).ksi),
            MES::ShapeFunctions::dEta::N4(MES::getIntegrationPoint(i).ksi)
        );
    }
    for(int i = 0; i < INTEGRATION_POINTS*INTEGRATION_POINTS; i++) {
        jacobian(A, MES::getIntegrationPoint(i));
        std::cout << std::endl;
    }
}
