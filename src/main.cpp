#include <cmath>
#include <iostream>

#include "MESCore.hpp"
#include "MESShapeFunctions.hpp"
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

std::tuple<MES::Matrix<4, 1>, MES::Matrix<4, 1>> calc_dN_dxdy(MES::Jakobian Jak, MES::node pc) {
    using MSF = MES::ShapeFunctions;
    std::array<MES::Matrix<2, 1>, 4> ret = {};
    auto dNdKsi = MSF::dKsi::NMatrix(pc);
    auto dNdEta = MSF::dEta::NMatrix(pc);
    for(int i = 0; i < 4; i++) {
        MES::Matrix<2, 1> a;
        a[0, 0] =  dNdKsi[i];
        a[1, 0] =  dNdEta[i];
        ret[i] = Jak.J1 * a;
    }
    return {
        MES::Matrix<4, 1>{
            ret[0][0, 0],
            ret[1][0, 0],
            ret[2][0, 0],
            ret[3][0, 0]
        },
        MES::Matrix<4, 1>{
            ret[0][1, 0],
            ret[1][1, 0],
            ret[2][1, 0],
            ret[3][1, 0]
        },
    };
}

int main() {
    MES::Test::testMatrix();
    MES::grid A;
    //MES::parseTextFile(A, "Test2_4_4_Mixgrid.txt");
    MES::parseTextFile(A, "Test1_4_4.txt");
    //std::cout << "After Grid Get: " << A.numNodes << " " << A.numElements << "\n";
    /*
     {
     for(int i = 0; i < A.numNodes; i++) {
         //std::cout << i << ": " << A.Node[i].x << " " << A.Node[i].y << " | BC: " << A.Node[i].BC << std::endl;
         }
     for(int i = 0; i < A.numElements; i++) {
         //std::cout << i << ": " << A.Element[i].ID[0] << " " << A.Element[i].ID[1] << " "
             << A.Element[i].ID[2] << " " << A.Element[i].ID[3] << std::endl;
     }
     //std::cout << "----\n";
     for(int i = 0; i < 4; i++) {
         auto intP = MES::getIntegrationPoint(i);
         printf("Integration point %i: %lf, %lf\n", i, intP.x, intP.y);
     }
     //std::cout << "---- dN/dKsi\n";
     for(int i = 0; i < INTEGRATION_POINTS*INTEGRATION_POINTS; i++) {
         printf("pc%i:  %.5lf  %.5lf  %.5lf  %.5lf\n",
             i,
             MES::ShapeFunctions::dKsi::N1(MES::getIntegrationPoint(i).eta),
             MES::ShapeFunctions::dKsi::N2(MES::getIntegrationPoint(i).eta),
             MES::ShapeFunctions::dKsi::N3(MES::getIntegrationPoint(i).eta),
             MES::ShapeFunctions::dKsi::N4(MES::getIntegrationPoint(i).eta)
         );
     }
     //std::cout << "---- dN/dEta\n";
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

    for(auto& el : A.Element) {
        MES::Matrix<4, 4> H;
        MES::Matrix<4, 4> C;
        for(int i = 0; i < INTEGRATION_POINTS*INTEGRATION_POINTS; i++) {
            auto point = MES::IntegrationPoint::get(i);
            auto weight = MES::IntegrationPoint::getWeight(i);
            auto k = MES::Jakobian(point, el);
            auto [dNdx, dNdy] = calc_dN_dxdy(k, point);
            MES::Matrix<4, 1> Nmat = MES::ShapeFunctions::NMatrix(point);

            auto H_i =
                (dNdx * dNdx.transpose() + dNdy * dNdy.transpose());
            H_i *= (k.detJ * GLOB.Conductivity);
            auto C_i = Nmat * Nmat.transpose();
            C_i *= GLOB.Density * GLOB.SpecificHeat * k.detJ;
            H = H + H_i * (weight.x * weight.y);
            C = C + C_i * (weight.x * weight.y);
            el.H = H;
            el.C = C;
        }
    }

    //std::cout << H_glob;
    std::cout << "H glob:\n" <<  H_glob;

    // H_bc
    MES::Matrix<16, 1> vP_glob;

    for(auto& el: A.Element) {
        auto elemID = find(A.Element.begin(), A.Element.end(), el) - A.Element.begin();
        auto& nodeVec = el.ID;
        using MSF = MES::ShapeFunctions;
        MES::Matrix<4, 4> H_BC;
        MES::Matrix<4, 1> vP;



        for(int i = 0; i < 4; i++) { // iterate each side
            auto index1 = (i+0)%4, index2 = (i+1)%4;
            auto n1 = nodeVec[index1], n2 = nodeVec[index2];
            // no boundary condition
            if (!el.Node[index1]->BC || !el.Node[index2]->BC)
                continue;
            printf("Node-y boku %i elem. %u: %i, %i\n", i, elemID, n1+1, n2+1);
            printf("%i, %i\n", el.Node[index1]->BC, el.Node[index2]->BC);
            MES::Matrix<4, 4> H_BC_L;
            MES::Matrix<4, 1> Vec_P;
            for(int j = 0; j < INTEGRATION_POINTS; j++) {
                auto pc = MES::IntegrationPoint::getForSide(j, i);
                MES::Matrix<4, 1> N = MSF::NMatrix(pc);
                auto weight = MES::IntegrationPoint::getWeightForSide(j);
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
            //std::cout << detJ << std::endl;
            H_BC = H_BC + H_BC_L;
            vP = vP + Vec_P;
            //std::cout << H_BC_L << std::endl;
            //std::cout << Vec_P;
            for(int i = 0; i < el.ID.size(); i++) {
                vP_glob[el.ID[i], 0] += vP[i, 0];
            }
            std::cout << "BBBB\n";
        }
        std::cout << "Element " << elemID + 1 << std::endl;
        std::cout << H_BC << std::endl;
        std::cout << "H elem: " << el.H << '\n';
        std::cout << "H elem + HBC: " << (el.H + H_BC) << '\n';
        el.H += H_BC;

        //std::cout << vP << std::endl;
    }


    int elemNum = 0;
    for(auto& el : A.Element) {
        std::cout << "Element id " << ++elemNum << " H: " << el.H << '\n';
        for(int a = 0; a < 4; a++) {
            auto currElemID = el.ID;
            for(int b = 0; b < 4; b++) {
                H_glob[currElemID[a], currElemID[b]] += el.H[a, b];
                C_glob[currElemID[a], currElemID[b]] += el.C[a, b];
            }
        }
            //std::cout << H_glob;
    }


    std::cout << "H glob + Hbc:\n" <<  H_glob << '\n';

    std::cout << "C glob: " << C_glob << '\n';

    std::cout << "Wektor Pglob:\n" << vP_glob.transpose() << std::endl;

    auto tauInf = MES::gauss(H_glob, vP_glob).transpose();
    std::cout << "temp. koncowa: " << tauInf << '\n';



    MES::Matrix<16, 1> t0(GLOB.InitialTemp);

    for(MES::real_t deltaTau = 50; deltaTau <= GLOB.SimulationTime; deltaTau += GLOB.SimulationStepTime) {
        auto a = (H_glob + C_glob / deltaTau);
        auto b = -1 * (-(C_glob / deltaTau) * t0 - vP_glob );
        auto Hg_Cg_dt = H_glob + (C_glob / deltaTau);
        auto Pg_Cg_dt_t0 = vP_glob + (C_glob/deltaTau)*t0;
        //auto t1 = MES::gauss(Hg_Cg_dt, -1*Pg_Cg_dt_t0);
        auto t1 = MES::gauss(a, b);
        std::cout << "Time: " << deltaTau << std::endl;
        std::cout << "min: " << t1.min() << " max: " << t1.max() << std::endl;
        std::cout << "[H]+[C]/d: " << a << '\n';
        std::cout << "P}+{[C]/dT}*{T0}: " << b.transpose() << '\n';
        std::cout << " temperatures: " << t1.transpose() << std::endl;
        t0 = t1;
        break;
    }




}
