

#include "MESCalculus.hpp"
#include "MESCore.hpp"
#ifdef DEBUG
#include <iostream>
#include <ostream>
#endif
MES::real_t MES::CalculusHelper::_integrate1D(integrable1D ff, int numPoints) {
    real_t sum = 0.0;
    numPoints--;
    for(int i = 0; i <= numPoints; i++) {
        sum +=
            GaussQuadratureTables::point_weight[numPoints][i] *
            ff(GaussQuadratureTables::point_x[numPoints][i]);
        #ifdef DEBUG
        printf("\ni: %i, weight: %lf, x: %lf, w*ff(x): %lf",
            i, GaussQuadratureTables::point_weight[numPoints][i],
            GaussQuadratureTables::point_x[numPoints][i],
            GaussQuadratureTables::point_weight[numPoints][i]*ff(GaussQuadratureTables::point_x[numPoints][i])
        );
        #endif
    }
    #ifdef DEBUG
    printf("\n");
    #endif
    return sum;
}

MES::real_t MES::CalculusHelper::_integrate2D(integrable2D ff, int numPoints) {
    real_t sum = 0.0;
    numPoints--;
    for(int i = 0; i <= numPoints; i++) {
        for(int j = 0 ; j <= numPoints; j++) {
            sum +=
                GaussQuadratureTables::point_weight[numPoints][i] * GaussQuadratureTables::point_weight[numPoints][j] *
                ff(GaussQuadratureTables::point_x[numPoints][i], GaussQuadratureTables::point_x[numPoints][j]);
            #ifdef DEBUG
            printf("\ni: %i, j: %i, weight i: %lf, weight j: %lf, x: %lf, y: %lf w*ff(x): %lf",
                i, j, GaussQuadratureTables::point_weight[numPoints][i], GaussQuadratureTables::point_weight[numPoints][j],
                GaussQuadratureTables::point_x[numPoints][i], GaussQuadratureTables::point_x[numPoints][j],
                GaussQuadratureTables::point_weight[numPoints][i]*GaussQuadratureTables::point_weight[numPoints][j]*
                ff(GaussQuadratureTables::point_x[numPoints][i], GaussQuadratureTables::point_x[numPoints][j])
            );
            #endif
        }
    }
    #ifdef DEBUG
    printf("\n");
    #endif
    return sum;
}

MES::node MES::IntegrationPoint::get(int n, int numPoints) {
    int k = 0;
    for(int i = 0; i < numPoints; i++) {
        for(int j = 0; j < numPoints; j++) {
            if(k++ >= n) {
                return MES::node{
                    MES::GaussQuadratureTables::point_x[numPoints-1][i],
                    MES::GaussQuadratureTables::point_x[numPoints-1][j]
                };
            }
        }
    }
    return MES::node{0,0};
}

MES::node MES::IntegrationPoint::getForSide(int n, int side, int numPoints) {
    MES::node ret = MES::node{0, 0};
    ret = MES::node {
        MES::GaussQuadratureTables::point_x[numPoints-1][n],
        MES::GaussQuadratureTables::point_x[numPoints-1][n]
    };
    switch(side) {
        case SIDE_BOTTOM: // dolny bok
            ret.eta = -1;
        break;
        case SIDE_RIGHT: // prawy bok
            ret.ksi = 1;
        break;
        case SIDE_TOP: // górny bok
            ret.eta = 1;
        break;
        case SIDE_LEFT:
        default: // lewy bok
            ret.ksi = -1;
        break;
    }

    return ret;
}


MES::node MES::IntegrationPoint::getWeight(int n, int numPoints) {
    int k = 0;
    for(int i = 0; i < numPoints; i++) {
        for(int j = 0; j < numPoints; j++) {
            if(k++ >= n) {
                return MES::node{
                    MES::GaussQuadratureTables::point_weight[numPoints-1][i],
                    MES::GaussQuadratureTables::point_weight[numPoints-1][j]
                };
            }
        }
    }
    return MES::node{0,0};
}

MES::real_t MES::IntegrationPoint::getWeightForSide(int n, int numPoints) {
    return MES::GaussQuadratureTables::point_weight[numPoints-1][n];
}

MES::Jakobian::Jakobian(node integrationP, element E) {
    using MSF = MES::ShapeFunctions;
    J = {
        MSF::dKsi::N1(integrationP) * E.Node[0]->x +
        MSF::dKsi::N2(integrationP) * E.Node[1]->x +
        MSF::dKsi::N3(integrationP) * E.Node[2]->x +
        MSF::dKsi::N4(integrationP) * E.Node[3]->x,

        MSF::dKsi::N1(integrationP) * E.Node[0]->y +
        MSF::dKsi::N2(integrationP) * E.Node[1]->y +
        MSF::dKsi::N3(integrationP) * E.Node[2]->y +
        MSF::dKsi::N4(integrationP) * E.Node[3]->y,

        MSF::dEta::N1(integrationP) * E.Node[0]->x +
        MSF::dEta::N2(integrationP) * E.Node[1]->x +
        MSF::dEta::N3(integrationP) * E.Node[2]->x +
        MSF::dEta::N4(integrationP) * E.Node[3]->x,

        MSF::dEta::N1(integrationP) * E.Node[0]->y +
        MSF::dEta::N2(integrationP) * E.Node[1]->y +
        MSF::dEta::N3(integrationP) * E.Node[2]->y +
        MSF::dEta::N4(integrationP) * E.Node[3]->y
    };
    detJ =  J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0];

    J1[0, 0] = J[1, 1];
    J1[0, 1] = -J[0, 1];
    J1[1, 0] = -J[1, 0];
    J1[1, 1] = J[0, 0];
    J1 *= 1/detJ;
}
