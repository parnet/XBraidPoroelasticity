/*
 * Copyright (c) 2019-2020:  G-CSC, Goethe University Frankfurt
 * Author: Arne Naegel
 *
 * This file is part of UG4.
 *
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "barry_mercer.h"


namespace ug {
    namespace XBraidPoroelasticity {
        const double BarryMercerNondimensional::X0 = 0.25;
        const double BarryMercerNondimensional::Y0 = 0.25;
        const double BarryMercerNondimensional::m_PI = ug::PI;
        double BarryMercerNondimensional::source_strength = 2.0;
        size_t BarryMercerNondimensional::NAPPROX = 64;

//! Computes coefficient from Eq. (24) in Barry & Mercer, ACME, 1999 (for $\omega=1)
        double BarryMercerNondimensional::FourierCoeff_P(int n, int q, double t_hat) const {

            double coeff_factor = source_strength;
            //  double beta = BARRY_MERCER_DATA.BETA
            double x0 = X0;
            double y0 = Y0;

            if ((n % 4 == 0) || (q % 4 == 0)) { return 0.0; }

            double lambda_n = n * m_PI;
            double lambda_q = q * m_PI;
            double _lambda_nq = lambda_n * lambda_n + lambda_q * lambda_q;

            double val1 = - coeff_factor * sin(lambda_n * x0) * sin(lambda_q * y0);
            double val2 = (_lambda_nq * sin(t_hat) - cos(t_hat) + exp(-_lambda_nq * t_hat));
            double val3 = (1 + _lambda_nq * _lambda_nq);

            return (val1 * val2) / val3;
        }


//! double t_hat = m_beta*t;

/*! Pressure is normalized and should be multiplied by (BARRY_MERCER_DATA.LAMBDA+2.0*BARRY_MERCER_DATA.MU) */
        double BarryMercerNondimensional::Pressure2D(double x, double y, double t_hat) const {
            int N = NAPPROX;
            double pp = 0.0;
            for (int m = 2; m <= N; ++m) {
                for (int n = 1; n <= m - 1; ++n) {
                    int q = m - n;

                    double sinnx = sin(m_PI * n * x);
                    double sinqy = sin(m_PI * q * y);
                    double _coeff_nq = FourierCoeff_P(n, q, t_hat);
                    pp += _coeff_nq * sinnx * sinqy;
                }
            }
            return -4.0 * pp;
        }

        double BarryMercerNondimensional::VelX2D(double x, double y, double t_hat) const {
            int N = NAPPROX;
            double u = 0.0;
            for (int m = 2; m <= N; ++m) {
                for (int n = 1; n <= m - 1; ++n) {
                    int q = m - n;

                    double _lambda_n = m_PI * n;
                    double _lambda_q = m_PI * q;
                    double _lambda_nq = (m_PI * m_PI) * (n * n + q * q);

                    double _coskx = cos(_lambda_n * x) ;
                    double _sinqy = sin(_lambda_q * y);

                    double _coeff_kq = FourierCoeff_P(n, q, t_hat);

                    u +=  _coeff_kq * _coskx * _sinqy  * _lambda_n/ _lambda_nq;
                }
            }
            return 4.0 * u;
        }


        double BarryMercerNondimensional::VelY2D(double x, double y, double t_hat) const {
            int N = NAPPROX;
            double u = 0.0;
            for (int m = 2; m <= N; ++m) {
                for (int n = 1; n <= m - 1; ++n) {
                    int q = m - n;
                    double _lambda_q = m_PI * q;
                    double _lambda_n = m_PI * n;
                    double _lambda_nq = (m_PI * m_PI) * (n * n + q * q);

                    double _cosqy = cos(_lambda_q * y) ;
                    double _sinnx = sin(_lambda_n * x);

                    double _coeff_nq = FourierCoeff_P(n, q, t_hat);
                    u += _coeff_nq * _sinnx * _cosqy * _lambda_q / _lambda_nq;
                }
            }
            return 4.0 * u;
        }
    }
}




