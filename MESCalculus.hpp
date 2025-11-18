#ifndef _MES_CALCULUS_H
#define _MES_CALCULUS_H
#include <functional>

#include "MESCore.hpp"
#include "MESGaussQuadratures.h"
#include "MESMatrix.hpp"


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

    node getIntegrationPointForSide(int n, int side);
    node getIntegrationPoint(int n);
    node getIntegrationPointWeight(int n);

    struct ShapeFunctions {
        real_t static N1(real_t ksi, real_t eta);
        real_t static N2(real_t ksi, real_t eta);
        real_t static N3(real_t ksi, real_t eta);
        real_t static N4(real_t ksi, real_t eta);
        struct dKsi {
            real_t static N1(real_t eta);
            real_t static N2(real_t eta);
            real_t static N3(real_t eta);
            real_t static N4(real_t eta);
        };
        struct dEta {
            real_t static N1(real_t ksi);
            real_t static N2(real_t ksi);
            real_t static N3(real_t ksi);
            real_t static N4(real_t ksi);
        };
    };

    struct Jakobian {
        Matrix<2, 2> J;
        Matrix<2, 2> J1;
        real_t detJ;
    };
}
#endif // _MES_CALCULUS_H
