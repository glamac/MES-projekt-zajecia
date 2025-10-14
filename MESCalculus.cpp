#include "MESCalculus.hpp"

MES::real_t MES::CalculusHelper::integrate1D(integrable1D ff, int numPoints) {
    real_t sum = 0.0;
    for(int i = 1; i <= numPoints; i++) {
        sum += 
            GaussQuadratureTables::point_weight[numPoints][i] * 
            ff(GaussQuadratureTables::point_x[numPoints][i]);
    }
    return sum/2.0;
}