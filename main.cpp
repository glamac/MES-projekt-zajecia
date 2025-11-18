#include <iostream>

#include "MESCore.hpp"
#include "MESTests.hpp"
#include "MESTypes.hpp"
#include "MESCalculus.hpp"
#include "MESMatrix.hpp"

MES::real_t testfunc(MES::real_t x) {
    return 5*x*x + 3*x + 6;
}

MES::real_t testfunc2D(MES::real_t x, MES::real_t y) {
    return 5*x*x*y*y + 3*x*y + 6;
}

MES::Jakobian jacobian(MES::grid& A, MES::node p1, int element = 0) {
    using MSF = MES::ShapeFunctions;

    MES::Jakobian j;

    j.J = {
        MSF::dKsi::N1(p1.eta) * A.Node[A.Element[element].ID[0]].x +
        MSF::dKsi::N2(p1.eta) * A.Node[A.Element[element].ID[1]].x +
        MSF::dKsi::N3(p1.eta) * A.Node[A.Element[element].ID[2]].x +
        MSF::dKsi::N4(p1.eta) * A.Node[A.Element[element].ID[3]].x,
        MSF::dKsi::N1(p1.eta) * A.Node[A.Element[element].ID[0]].y +
        MSF::dKsi::N2(p1.eta) * A.Node[A.Element[element].ID[1]].y +
        MSF::dKsi::N3(p1.eta) * A.Node[A.Element[element].ID[2]].y +
        MSF::dKsi::N4(p1.eta) * A.Node[A.Element[element].ID[3]].y,
        MSF::dEta::N1(p1.ksi) * A.Node[A.Element[element].ID[0]].x +
        MSF::dEta::N2(p1.ksi) * A.Node[A.Element[element].ID[1]].x +
        MSF::dEta::N3(p1.ksi) * A.Node[A.Element[element].ID[2]].x +
        MSF::dEta::N4(p1.ksi) * A.Node[A.Element[element].ID[3]].x,
        MSF::dEta::N1(p1.ksi) * A.Node[A.Element[element].ID[0]].y +
        MSF::dEta::N2(p1.ksi) * A.Node[A.Element[element].ID[1]].y +
        MSF::dEta::N3(p1.ksi) * A.Node[A.Element[element].ID[2]].y +
        MSF::dEta::N4(p1.ksi) * A.Node[A.Element[element].ID[3]].y
    };
    
    j.detJ = j.J[0, 0] * j.J[1, 1] - j.J[0, 1]*j.J[1, 0];

    j.J1[0, 0] = j.J[1, 1];
    j.J1[0, 1] = -j.J[0, 1];
    j.J1[1, 0] = -j.J[1, 0];
    j.J1[1, 1] = j.J[0, 0];

    j.J1 *= 1/j.detJ;
    
    return j;

}

MES::real_t dNdKsi (int N, MES::node pc) {
    using MSF = MES::ShapeFunctions;
    std::function<MES::real_t(MES::real_t)> NdKsi[4] = {
        MSF::dKsi::N1,
        MSF::dKsi::N2,
        MSF::dKsi::N3,
        MSF::dKsi::N4
    };
    return NdKsi[N](pc.eta);
}

MES::real_t dNdEta (int N, MES::node pc) {
    using MSF = MES::ShapeFunctions;
    std::function<MES::real_t(MES::real_t)> NdEta[4] = {
        MSF::dEta::N1,
        MSF::dEta::N2,
        MSF::dEta::N3,
        MSF::dEta::N4
    };
    return NdEta[N](pc.ksi);
}

std::array<MES::Matrix<2, 1>, 4> calc_dN_dxdy(MES::Jakobian Jak, MES::node pc) {
    using MSF = MES::ShapeFunctions;
    std::array<MES::Matrix<2, 1>, 4> ret = {};
    for(int i = 0; i < 4; i++) {
        MES::Matrix<2, 1> a;
        a[0, 0] =  dNdKsi(i, pc);
        a[1, 0] =  dNdEta(i, pc);
        ret[i] = Jak.J1 * a;
    }
    return ret;
}

int main() {
    MES::Test::testMatrix();
    MES::grid A;
    //MES::parseTextFile(A, "Test2_4_4_Mixgrid.txt");
    MES::parseTextFile(A, "Test1_4_4.txt");
    std::cout << "After Grid Get: " << A.numNodes << " " << A.numElements << "\n";
    /*
    for(int i = 0; i < A.numNodes; i++) {
        std::cout << i << ": " << A.Node[i].x << " " << A.Node[i].y << " | BC: " << A.Node[i].BC << std::endl;
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
        */
    printf("numElements: %u", GLOB.numElements);

    MES::Matrix<16, 16> H_glob;
    
    for(int j = 0; j < GLOB.numElements; j++) {
        MES::Matrix<4, 4> H;
        for(int i = 0; i < INTEGRATION_POINTS*INTEGRATION_POINTS; i++) {
            auto point = MES::getIntegrationPoint(i);
            auto weight = MES::getIntegrationPointWeight(i);
            auto k = jacobian(A, point, j);
            auto dN_dxdy = calc_dN_dxdy(k, point);
            MES::Matrix<4, 1> dNdx = {
                dN_dxdy[0][0, 0],
                dN_dxdy[1][0, 0],
                dN_dxdy[2][0, 0],
                dN_dxdy[3][0, 0]
            };
            MES::Matrix<4, 1> dNdy = {
                dN_dxdy[0][1, 0],
                dN_dxdy[1][1, 0],
                dN_dxdy[2][1, 0],
                dN_dxdy[3][1, 0]
            };

            auto H_i = 
                (dNdx * dNdx.transpose() + dNdy * dNdy.transpose());
            H_i *= (k.detJ * GLOB.Conductivity);
            H = H + H_i * (weight.x * weight.y);
        }
        for(int a = 0; a < 4; a++) {
            auto currElemID = A.Element[j].ID;
            for(int b = 0; b < 4; b++) {
                H_glob[currElemID[a], currElemID[b]] += H[a, b];
            }
        }
            std::cout << H_glob;
        std::cout << "\nMacierz H:\n" << H << std::endl;
        std::cout << "\nNr. Elem: " << j << std::endl;
    }

    std::cout << H_glob;

    // H_bc

    auto& el = A.Element[0];
    auto& nodeVec = el.ID;
    for(int i = 0, j = 1; i < 4; ++i, j=(j+1)%4) {
        auto n1 = nodeVec[(i+2)%4], n2 = nodeVec[(j+2)%4];
        printf("Node-y boku %u: %i, %i\n", i, n1+1, n2+1);
        auto e = MES::getIntegrationPoint(i);
        std::cout << e.ksi << " " <<  e.eta << std::endl;
    }

    

}
