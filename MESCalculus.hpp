#include <functional>

#include "MESCore.hpp"
#include "MESGaussQuadratures.h"

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
}