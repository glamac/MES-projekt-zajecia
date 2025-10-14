#include "MESCalculus.hpp"
#include <iostream>
MES::real_t MES::CalculusHelper::_integrate1D(integrable1D ff, int numPoints) {
    real_t sum = 0.0;
    numPoints--;
    for(int i = 0; i <= numPoints; i++) {
        sum += 
            GaussQuadratureTables::point_weight[numPoints][i] * 
            ff(GaussQuadratureTables::point_x[numPoints][i]);
        printf("\ni: %i, weight: %lf, x: %lf, w*ff(x): %lf",
            i, GaussQuadratureTables::point_weight[numPoints][i],
            GaussQuadratureTables::point_x[numPoints][i],
            GaussQuadratureTables::point_weight[numPoints][i]*ff(GaussQuadratureTables::point_x[numPoints][i])
        );
    }
    printf("\n");
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
            printf("\ni: %i, j: %i, weight i: %lf, weight j: %lf, x: %lf, y: %lf w*ff(x): %lf",
                i, j, GaussQuadratureTables::point_weight[numPoints][i], GaussQuadratureTables::point_weight[numPoints][j],
                GaussQuadratureTables::point_x[numPoints][i], GaussQuadratureTables::point_x[numPoints][j],
                GaussQuadratureTables::point_weight[numPoints][i]*GaussQuadratureTables::point_weight[numPoints][j]*
                ff(GaussQuadratureTables::point_x[numPoints][i], GaussQuadratureTables::point_x[numPoints][j])
            );
        }
    }
    printf("\n");
    return sum;
}