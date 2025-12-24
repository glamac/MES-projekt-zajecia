#ifndef _MES_CALCULUS_H
#define _MES_CALCULUS_H
#include <functional>

#include "MESCore.hpp"
#include "MESGaussQuadratures.h"
#include "MESMatrix.hpp"
#include "MESShapeFunctions.hpp"

#define SIDE_BOTTOM 0
#define SIDE_RIGHT 1
#define SIDE_TOP 2
#define SIDE_LEFT 3
#define DEFAULT_INTEGRATION_POINTS 2

namespace MES {
    typedef std::function<real_t(real_t)> integrable1D;
    typedef std::function<real_t(real_t, real_t)> integrable2D;

	struct CalculusHelper {
        real_t static integrate(integrable1D ff, int numPoints = 3) {
            return _integrate1D(ff, numPoints);
        };
        real_t static integrate(integrable2D ff, int numPoints = 3) {
            return _integrate2D(ff, numPoints);
        };

    private:
    	real_t static _integrate1D(integrable1D ff, int numPoints = 3);
        real_t static _integrate2D(integrable2D ff, int numPoints = 3);
	};
	struct IntegrationPoint {
        node static get(int n, int numPoints = DEFAULT_INTEGRATION_POINTS);
        node static getWeight(int n, int numPoints = DEFAULT_INTEGRATION_POINTS);
    	node static getForSide(int n, int side, int numPoints = DEFAULT_INTEGRATION_POINTS);
        real_t static getWeightForSide(int n, int numPoints = DEFAULT_INTEGRATION_POINTS);
	};

    struct Jakobian {
        Matrix<2, 2> J;
        Matrix<2, 2> J1;
        real_t detJ;
        Jakobian(node integrationP, element E);
    };
}
#endif // _MES_CALCULUS_H
