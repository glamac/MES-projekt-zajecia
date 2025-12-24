#ifndef _MES_SHAPEFUNCTIONS_H
#define _MES_SHAPEFUNCTIONS_H
#include "MESCore.hpp"
#include "MESMatrix.hpp"

namespace MES {
    struct ShapeFunctions {
        // Ksi, Eta
        real_t constexpr static N1(real_t ksi, real_t eta) {
            return 0.25 * (1 - ksi) * (1 - eta);
        };
        real_t constexpr static N2(real_t ksi, real_t eta){
            return 0.25 * (1 + ksi) * (1 - eta);
        };
        real_t constexpr static N3(real_t ksi, real_t eta){
            return 0.25 * (1 + ksi) * (1 + eta);
        };
        real_t constexpr static N4(real_t ksi, real_t eta){
            return 0.25 * (1 - ksi) * (1 + eta);
        };
        Matrix<4, 1> constexpr static NMatrix(real_t ksi, real_t eta) {
            return Matrix<4, 1>{
                N1(ksi, eta), N2(ksi, eta), N3(ksi, eta), N4(ksi, eta)
            };
        };
        // Point
        real_t constexpr static N1(node p) {return N1(p.ksi, p.eta);};
        real_t constexpr static N2(node p) {return N2(p.ksi, p.eta);};
        real_t constexpr static N3(node p) {return N3(p.ksi, p.eta);};
        real_t constexpr static N4(node p) {return N4(p.ksi, p.eta);};
        Matrix<4, 1> constexpr static NMatrix(node p) {return NMatrix(p.ksi, p.eta);};
        struct dKsi {
            real_t constexpr static N1(real_t eta){
                return -0.25 * (1 - eta);
            };
            real_t constexpr static N2(real_t eta){
                return 0.25 * (1 - eta);
            };
            real_t constexpr static N3(real_t eta){
                return 0.25 * (1 + eta);
            };
            real_t constexpr static N4(real_t eta){
                return -0.25 * (1 + eta);
            };
            Matrix<4, 1> constexpr static NMatrix(real_t eta) {
                return Matrix<4, 1>{
                    N1(eta), N2(eta), N3(eta), N4(eta)
                };
            };
            real_t constexpr static N1(node p) {return N1(p.eta);};
            real_t constexpr static N2(node p) {return N2(p.eta);};
            real_t constexpr static N3(node p) {return N3(p.eta);};
            real_t constexpr static N4(node p) {return N4(p.eta);};
            Matrix<4, 1> constexpr static NMatrix(node p) {return NMatrix(p.eta);};
        };
        struct dEta {
            real_t constexpr static N1(real_t ksi){
                return -0.25 * (1 - ksi);
            };
            real_t constexpr static N2(real_t ksi){
                return -0.25 * (1 + ksi);
            };
            real_t constexpr static N3(real_t ksi){
                return 0.25 * (1 + ksi);
            };
            real_t constexpr static N4(real_t ksi){
                return 0.25 * (1 - ksi);
            };
            Matrix<4, 1> constexpr static NMatrix(real_t ksi) {
                return Matrix<4, 1>{
                    N1(ksi), N2(ksi), N3(ksi), N4(ksi)
                };
            };
            real_t constexpr static N1(node p) {return N1(p.ksi);};
            real_t constexpr static N2(node p) {return N2(p.ksi);};
            real_t constexpr static N3(node p) {return N3(p.ksi);};
            real_t constexpr static N4(node p) {return N4(p.ksi);};
            Matrix<4, 1> constexpr static NMatrix(node p) {return NMatrix(p.ksi);};
        };
    };
}

#endif // _MES_SHAPEFUNCTIONS_H
