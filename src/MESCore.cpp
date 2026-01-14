#include "MESCore.hpp"
#include "MESCalculus.hpp"
#include "MESShapeFunctions.hpp"
#include "MESTypes.hpp"
#include <cmath>

std::tuple<MES::Matrix, MES::Matrix> calc_dN_dxdy(MES::Jakobian Jak, MES::node pc) {
    using MSF = MES::ShapeFunctions;
    std::array<MES::Matrix, 4> ret = {};
    auto dNdKsi = MSF::dKsi::NMatrix(pc);
    auto dNdEta = MSF::dEta::NMatrix(pc);
    for(int i = 0; i < 4; i++) {
        MES::Matrix a = MES::Matrix(2, 1);
        a[0, 0] =  dNdKsi[i];
        a[1, 0] =  dNdEta[i];
        ret[i] = Jak.J1 * a;
    }
    return {
        MES::Matrix(4, 1, {
            ret[0][0, 0],
            ret[1][0, 0],
            ret[2][0, 0],
            ret[3][0, 0]
        }),
        MES::Matrix(4, 1, {
            ret[0][1, 0],
            ret[1][1, 0],
            ret[2][1, 0],
            ret[3][1, 0]
        }),
    };
}

MES::element& MES::element::calculateMatrices() {
	using MIP = MES::IntegrationPoint;
	using MSF = MES::ShapeFunctions;

	const int numPoints = GLOB.integrationPoints_Surface;
	element& thisElement = (*this);

	for(int i = 0; i < numPoints*numPoints; i++) {
		auto iPoint = MIP::get(i);
		auto pWeight = MIP::getWeight(i);
		auto Jak = MES::Jakobian(iPoint, thisElement);
		auto [dNdx, dNdy] = calc_dN_dxdy(Jak, iPoint);
		// TODO: cache this and do not recalculate each time god damn...
		Matrix NMatrix = MSF::NMatrix(iPoint);

		const real_t weight = (pWeight.x * pWeight.y);

		auto H_i = (dNdx * dNdx.transpose() + dNdy * dNdy.transpose());
		H_i *= Jak.detJ * GLOB.Conductivity * weight;
		H += H_i;

		auto C_i = (NMatrix * NMatrix.transpose());
		C_i *= GLOB.Density * GLOB.SpecificHeat * Jak.detJ * weight;
		C += C_i;

	}
	return thisElement;
}

MES::element& MES::element::calculateBoundaryCondition() {
	using MIP = MES::IntegrationPoint;
	using MSF = MES::ShapeFunctions;

	const int numPoints = GLOB.integrationPoints_Boundary;
	element& thisElement = (*this);
	// Iterate over each side in this element.
	for(int i = 0; i < 4; i++) {
		// The two nodes of this side.
		auto node1 = Node[(i+0)%4];
		auto node2 = Node[(i+1)%4];
		if(!node1->BC || !node2->BC) {
			continue; // No boundary condition on this side.
		}
		Matrix H_BC_i(4, 4, {0});
		Matrix vP_i(4, 1, {0});
		for(int j = 0; j < numPoints; j++) {
			auto iPoint = MIP::getForSide(j, i);
			auto weight = MIP::getWeightForSide(j);
			auto NMatrix = MSF::NMatrix(iPoint);

			H_BC_i += weight * (NMatrix * NMatrix.transpose());
			vP_i += weight * NMatrix;
		}
		real_t detJ = (std::sqrt(
			std::pow(node1->x - node2->x, 2) +
			std::pow(node1->y - node2->y, 2)
		)) / 2.0;
		H_BC_i *= GLOB.Alpha * detJ;
		vP_i *= GLOB.Alpha * GLOB.Tot * detJ;

		H_BC += H_BC_i;
		vP += vP_i;
	}
	H += H_BC;
	return thisElement;
}

MES::element& MES::element::aggregateToGlobal() {
	element& thisElement = (*this);
	for(size_t i = 0; i < ID.size(); i++) {
		GLOB.vP[ID[i], 0] += vP[i, 0];
		for(size_t j = 0; j < ID.size(); j++) {
			GLOB.H[ID[i], ID[j]] += H[i, j];
			GLOB.C[ID[i], ID[j]] += C[i, j];
		}
	}
	return thisElement;
}
