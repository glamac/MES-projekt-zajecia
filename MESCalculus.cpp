

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

MES::node MES::getIntegrationPoint(int n) {
    int k = 0;
    for(int i = 0; i < INTEGRATION_POINTS; i++) {
        for(int j = 0; j < INTEGRATION_POINTS; j++) {
            if(k++ >= n) {
                return MES::node{
                    MES::GaussQuadratureTables::point_x[INTEGRATION_POINTS-1][i],
                    MES::GaussQuadratureTables::point_x[INTEGRATION_POINTS-1][j]
                };
            }
        }
    }
    return MES::node{0,0};
}

MES::node MES::getIntegrationPointForSide(int n, int side) {
    MES::node ret = MES::node{0, 0};
    ret = MES::node {
        MES::GaussQuadratureTables::point_x[INTEGRATION_POINTS-1][n],
        MES::GaussQuadratureTables::point_x[INTEGRATION_POINTS-1][n]
    };
    switch(side) {
        case 0: // dolny bok
            ret.eta = -1;
        break;
        case 1: // prawy bok
            ret.ksi = 1;
        break;
        case 2: // górny bok
            ret.eta = 1;
        break;
        default: // lewy bok
            ret.ksi = -1;
        break;
    }

    return ret;
}


MES::node MES::getIntegrationPointWeight(int n) {
    int k = 0;
    for(int i = 0; i < INTEGRATION_POINTS; i++) {
        for(int j = 0; j < INTEGRATION_POINTS; j++) {
            if(k++ >= n) {
                return MES::node{
                    MES::GaussQuadratureTables::point_weight[INTEGRATION_POINTS-1][i],
                    MES::GaussQuadratureTables::point_weight[INTEGRATION_POINTS-1][j]
                };
            }
        }
    }
    return MES::node{0,0};
}

MES::real_t MES::ShapeFunctions::N1(real_t ksi, real_t eta) {
    return 0.25 * (1 - ksi) * (1 - eta);
}
MES::real_t MES::ShapeFunctions::N2(real_t ksi, real_t eta) {
    return 0.25 * (1 + ksi) * (1 - eta);
}
MES::real_t MES::ShapeFunctions::N3(real_t ksi, real_t eta) {
    return 0.25 * (1 + ksi) * (1 + eta);
}
MES::real_t MES::ShapeFunctions::N4(real_t ksi, real_t eta) {
    return 0.25 * (1 - ksi) * (1 + eta);
}
MES::real_t MES::ShapeFunctions::dKsi::N1(real_t eta) {
    return -0.25 * (1 - eta);
};
MES::real_t MES::ShapeFunctions::dKsi::N2(real_t eta) {
    return 0.25 * (1 - eta);
};
MES::real_t MES::ShapeFunctions::dKsi::N3(real_t eta) {
    return 0.25 * (1 + eta);
};
MES::real_t MES::ShapeFunctions::dKsi::N4(real_t eta) {
    return -0.25 * (1 + eta);
};

MES::real_t MES::ShapeFunctions::dEta::N1(real_t ksi) {
    return -0.25 * (1 - ksi);
}
MES::real_t MES::ShapeFunctions::dEta::N2(real_t ksi) {
    return -0.25 * (1 + ksi);
}
MES::real_t MES::ShapeFunctions::dEta::N3(real_t ksi) {
    return 0.25 * (1 + ksi);
}
MES::real_t MES::ShapeFunctions::dEta::N4(real_t ksi) {
    return 0.25 * (1 - ksi);
}
