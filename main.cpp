#include <cmath>
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
     {
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
     }}
     */
    printf("numElements: %u", GLOB.numElements);

    MES::Matrix<16, 16> H_glob;
    MES::Matrix<16, 16> C_glob;

    for(int j = 0; j < GLOB.numElements; j++) {
        MES::Matrix<4, 4> H;
        MES::Matrix<4, 4> C;
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
            MES::Matrix<4, 1> Nmat = {
                MES::ShapeFunctions::N1(point.ksi, point.eta),
                MES::ShapeFunctions::N2(point.ksi, point.eta),
                MES::ShapeFunctions::N3(point.ksi, point.eta),
                MES::ShapeFunctions::N4(point.ksi, point.eta)
            };

            auto H_i =
                (dNdx * dNdx.transpose() + dNdy * dNdy.transpose());
            H_i *= (k.detJ * GLOB.Conductivity);
            auto C_i = Nmat * Nmat.transpose();
            C_i *= GLOB.Density * GLOB.SpecificHeat * k.detJ;
            H = H + H_i * (weight.x * weight.y);
            C = C + C_i * (weight.x * weight.y);
        }
        for(int a = 0; a < 4; a++) {
            auto currElemID = A.Element[j].ID;
            for(int b = 0; b < 4; b++) {
                H_glob[currElemID[a], currElemID[b]] += H[a, b];
                C_glob[currElemID[a], currElemID[b]] += C[a, b];
            }
        }
            std::cout << H_glob;
        std::cout << "\nMacierz H:\n" << H << std::endl;
        std::cout << "\nNr. Elem: " << j << std::endl;
    }

    std::cout << H_glob;

    // H_bc

    auto& el = A.Element[0];
    MES::Matrix<16, 1> vP_glob;

    for(auto& el: A.Element) {
        auto elemID = find(A.Element.begin(), A.Element.end(), el) - A.Element.begin();
        auto& nodeVec = el.ID;
        using MSF = MES::ShapeFunctions;
        MES::Matrix<4, 4> H_BC;
        MES::Matrix<4, 1> vP;


        for(int i = 0; i < 4; i++) {
            auto index1 = (i+2)%4, index2 = (i+3)%4;
            auto n1 = nodeVec[index1], n2 = nodeVec[index2];
            // no boundary condition
            if (!el.Node[index1]->BC || !el.Node[index2]->BC)
                continue;
            printf("Node-y boku %u elem. %u: %i, %i\n", i, elemID, n1+1, n2+1);
            printf("%i, %i\n", el.Node[index1]->BC, el.Node[index2]->BC);
            MES::Matrix<4, 4> H_BC_L;
            MES::Matrix<4, 1> Vec_P;
            for(int j = 0; j < INTEGRATION_POINTS; j++) {
                auto pc = MES::getIntegrationPointForSide(j, i);
                std::cout << pc.x << " | " << pc.y << " | " << j << std::endl;
                MES::Matrix<4, 1> N = {
                    MSF::N1(pc.ksi, pc.eta),
                    MSF::N2(pc.ksi, pc.eta),
                    MSF::N3(pc.ksi, pc.eta),
                    MSF::N4(pc.ksi, pc.eta),
                };
                std::cout << "N:" << N << std::endl;
                auto weight = MES::GaussQuadratureTables::point_weight[INTEGRATION_POINTS-1][j];
                H_BC_L = H_BC_L + weight * (N * N.transpose());
                Vec_P = Vec_P + weight * N;
            }
            auto p1 = el.Node[index1], p2 = el.Node[index2];
            // multiply by determinant (which is just length over 2)
            auto detJ = (std::sqrt(
                std::pow(p1->x - p2->x, 2) +
                std::pow(p1->y - p2->y, 2)
            ) )/ 2.0;
            H_BC_L *= GLOB.Alpha * detJ;
            Vec_P *= GLOB.Alpha * GLOB.Tot * detJ;
            std::cout << detJ << std::endl;
            H_BC = H_BC + H_BC_L;
            vP = vP + Vec_P;
            std::cout << H_BC_L << std::endl;
            std::cout << Vec_P;
            for(int i = 0; i < el.ID.size(); i++) {
                vP_glob[el.ID[i], 0] += vP[i, 0];
            }

            //aggregate h glob
            for(int a = 0; a < 4; a++) {
                for(int b = 0; b < 4; b++) {
                    auto ID_a = el.ID[a];
                    auto ID_b = el.ID[b];
                    H_glob[ID_a, ID_b] += H_BC[a, b];
                }
            }

        }
        std::cout << H_BC << std::endl;
        std::cout << vP << std::endl;
    }

    std::cout << vP_glob << std::endl;

    std::cout << H_glob;

    auto tauInf = MES::gauss(H_glob, vP_glob).transpose();
    std::cout << tauInf;

    MES::real_t deltaTau = 0;

    MES::Matrix<16, 1> t0(GLOB.InitialTemp);

    for(MES::real_t i = 0; i <= GLOB.SimulationTime; i += GLOB.SimulationStepTime) {
        deltaTau = i + 0.000001; // bo inaczej dziele przez 0 na początku
        auto a = (H_glob + C_glob / deltaTau);
        auto b = -1 * (-(C_glob / deltaTau) * t0 - vP_glob );
        //std::cout << "Matrix A: " << a << std::endl;
        //std::cout << "Matrix B: " << b << std::endl;
        //std::cout << "t0: " << t0 << std::endl << GLOB.InitialTemp << std::endl;
        auto t1 = MES::gauss(a, b);
        std::cout << "Time: " << deltaTau << std::endl;
        std::cout << "min: " << t1.min() << " max: " << t1.max() << std::endl;
        std::cout << " temperatures: " << t1.transpose() << std::endl;
        t0 = t1;
    }




}
