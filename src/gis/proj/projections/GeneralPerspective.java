/*
    Snyder Projection Implementation

    Copyright (C) 2012, 2013  Jeremy J. Gibbons

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
package gis.proj.projections;

import static gis.proj.SnyderMath.*;

import gis.proj.Azimuthal;
import gis.proj.Datum;
import gis.proj.Ellipsoid;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;


/**
 * All formulas adopted from USGS Professional Paper 1395.
 *
 * References:
 *
 * Map Projections - A Working Manual, USGS Professional Paper 1395
 * John P. Snyder
 * Pages 169-181
 *
 */

public final class GeneralPerspective implements Azimuthal {

    public String getName() {
        return "General Perspective";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double esq = ellip.getProperty("e^2");

        double h     = datum.getProperty("h"); // Usually 0
        double H     = datum.getProperty("H");
        double h0    = datum.getProperty("h0");
        double lon0  = datum.getProperty("lon0");
        double lat1  = datum.getProperty("lat1");
        double theta = datum.getProperty("theta");
        double x0    = datum.getProperty("x0");
        double y0    = datum.getProperty("y0");
        double gamma = datum.getProperty("gamma");
        double omega = datum.getProperty("omega");

        boolean tlt = datum.getProperty("useTilt")       == 1.0;
        boolean prj = datum.getProperty("useProjective") == 1.0;

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double coslat1 = StrictMath.cos(lat1);
        double sinlat1 = StrictMath.sin(lat1);

        double coslon0 = StrictMath.cos(lon0);
        double sinlon0 = StrictMath.sin(lon0);

        double cosgamma = StrictMath.cos(gamma);
        double singamma = StrictMath.sin(gamma);

        double cosomega = StrictMath.cos(omega);
        double sinomega = StrictMath.sin(omega);

        double costheta = StrictMath.cos(theta);
        double sintheta = StrictMath.sin(theta);

        double M, Q, Xt, Yt;

        if(ellip.isSphere()) {
            double P = H / a + 1.0;

            double P_M_one   = P - 1.0;
            double P_P_one   = P + 1.0;
            double aP_M_one  = a * P_M_one;
            double aaP_M_one = a * aP_M_one;

            double rhoCheck = a * (P - 1.0) / P;

            double rho, c, cosc, sinc;

            for(int i = 0; i < lon.length; i++) {
                Xt = x[i];
                Yt = y[i];

                if(tlt) {
                    M  = H * Xt / (H - Yt * sinomega);
                    Q  = H * Yt * cosomega / (H - Yt * sinomega);
                    Xt = M * cosgamma + Q * singamma;
                    Yt = Q * cosgamma - M * singamma;
                }

                rho = StrictMath.hypot(Xt, Yt);
                c = StrictMath.asin(
                        (P - StrictMath.sqrt(1.0 - rho * rho * P_P_one / aaP_M_one)) / (a * P_M_one / rho + rho / aP_M_one)
                        );

                if(P < 0.0 && rho > rhoCheck)
                    c -= PI_DIV_2;

                cosc = StrictMath.cos(c);
                sinc = StrictMath.sin(c);

                lat[i] = StrictMath.asin(cosc * sinlat1 + (Yt * sinc * coslat1 / rho));
                lon[i] = normalizeLonRad(lon0 + StrictMath.atan2(Xt * sinc, rho * coslat1 * cosc - Yt * sinlat1 * sinc));
            }
        } else {
            double one_M_esq = 1.0 - esq;

            double one_M_esqcoslat1coslat1 = 1.0 - esq * coslat1 * coslat1;

            double E_factor = - (h * h) / (a * a - a * a * esq);

            double N1            = a / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);
            double P_factor      = (H + N1 + h0) / a;
            double phig_dividend = N1 * esq * sinlat1 * coslat1;

            double L = 1.0 - esq * coslat1 * coslat1;
            double G = 1.0 - esq * sinlat1 * sinlat1;
            double J = 2.0 * esq * sinlat1 * coslat1;

            double LH_N_2 = L * H * -2.0;
            double G2     = G * 2.0;
            double HJ     = H * J;
            double LHH    = L * H * H;

            double h_over_a = h / a;

            double B, D, u, v, BIG_E, t, Kp, X, Y, S, P, phig, phig_old, lat1_M_phig, phi, phi_old, sinphi;
            double BIG_E_factor, esqsinphisinphi;
            int j;

            double F_num = 0.0, V_num = 0.0, M_num = 0.0, N_num = 0.0, W_num = 0.0, T_num = 0.0, K7_num = 0.0;

            if(prj) {
                F_num = sinlat1 * sinlon0 * cosgamma - coslon0 * singamma;
                V_num = sinlat1 * sinlon0 * singamma + coslon0 * cosgamma;
                M_num = sinlat1 * coslon0 * singamma - sinlon0 * cosgamma;
                N_num = sinlat1 * coslon0 * cosgamma + sinlon0 * singamma;
                W_num = -singamma * cosomega * costheta - cosgamma * sintheta;
                T_num = -singamma * cosomega * sintheta + cosgamma * costheta;
                K7_num = coslat1 * cosgamma * sinomega - sinlat1 * cosomega;// k[6];
            }

            double U, F, V, W, T, sinlat1_M_phig, N;
            double[] k = new double[11]; // K1..11 in text, or k0..10 in code
            double[] A = new double[16]; // A1..16 in text, or A0..15 in code
            double A15_DIV_A14, den;
            phi = 0;
            for(int i = 0; i < lon.length; ++i) {
                Xt = x[i];
                Yt = y[i];

                if(tlt || prj) {
                    M  = H * Xt / (H - Yt * sinomega);
                    Q  = H * Yt * cosomega / (H - Yt * sinomega);
                    Xt = M * cosgamma + Q * singamma;
                    Yt = Q * cosgamma - M * singamma;
                }

                phig = lat1;
                j = 0;
                do {
                    phig_old = phig;
                    P = (coslat1 / StrictMath.cos(phig)) * P_factor;
                    phig = lat1 - StrictMath.asin(phig_dividend / (P * a));
                } while(++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phig - phig_old));

                lat1_M_phig = lat1 - phig;

                B = P * StrictMath.cos(lat1_M_phig);
                D = P * StrictMath.sin(lat1_M_phig);

                u = LH_N_2 * B - G2 * D * Yt + B * J * Yt + D * HJ;
                v = LHH + G * Yt * Yt - HJ * Yt + one_M_esq * Xt * Xt;

                BIG_E = 1.0;
                t  = (P * P) * one_M_esqcoslat1coslat1 - BIG_E * one_M_esq;
                Kp = (-u + StrictMath.sqrt(u * u - 4.0 * t * v)) / (2.0 * t);
                S  = (Yt / Kp - D) * coslat1 + (B - H / Kp) * sinlat1;

                if(h != 0) {
                    phi = StrictMath.asin(S);
                    sinphi = StrictMath.sin(phi);
                    esqsinphisinphi = esq * sinphi * sinphi;
                    j = 0;
                    do {
                        phi_old = phi;
                        phi = StrictMath.asin(S / (one_M_esq / StrictMath.sqrt(1.0 - esqsinphisinphi) + h_over_a));

                        sinphi = StrictMath.sin(phi);

                        esqsinphisinphi = esq * sinphi * sinphi;
                        BIG_E_factor = (1.0 / StrictMath.sqrt(1.0 - esqsinphisinphi)) + h_over_a;

                        BIG_E = BIG_E_factor * BIG_E_factor - esqsinphisinphi * (1.0 / (1-esqsinphisinphi) + E_factor);//- (h * h) / (a * a - a * a * esq));
                    } while(++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phi - phi_old));
                }

                if(prj) {
                    sinlat1_M_phig = StrictMath.sin(lat1_M_phig);
                    U = P * ((sinlat1_M_phig) * cosgamma * sinomega + StrictMath.cos(lat1_M_phig) * cosomega);
                    F = F_num / U;
                    V = V_num * cosomega / U;
                    M = M_num * cosomega / U;
                    N = N_num / U;
                    W = W_num / U;
                    T = T_num / U;

                    k[4]  = -N * sinomega - coslat1 * coslon0 * cosomega / U;
                    k[5]  = -F * sinomega - coslat1 * sinlon0 * cosomega / U;
                    k[6]  = K7_num / U;
                    k[0]  = H * (M * costheta + N * sintheta) + k[4] * x0;
                    k[1]  = H * (V * costheta + F * sintheta) + k[5] * x0;
                    k[2]  = H * W * coslat1 + k[6] * x0;
                    k[3]  = H * W * P * sinlat1_M_phig + x0;
                    k[7]  = H * (M * sintheta - N * costheta) + k[4] * y0;
                    k[8]  = H * (V * sintheta - F * costheta) + k[5] * y0;
                    k[9]  = H * T * coslat1 + k[6] * y0;
                    k[10] = H * T * P * sinlat1_M_phig + y0;

                    A[0] = x[i] * k[4] - k[0];
                    A[1] = x[i] * k[5] - k[1];
                    A[2] = x[i] * k[6] - k[2];
                    A[3] = k[3] - x[i];
                    A[4] = y[i] * k[4] - k[7];
                    A[5] = y[i] * k[5] - k[8];
                    A[6] = y[i] * k[6] - k[9];
                    A[7] = k[10] - y[i];
                    A[8] = A[0] * A[7] - A[3] * A[4];
                    A[9] = A[0] * A[6] - A[2] * A[4];
                    A[10] = A[1] * A[4] - A[0] * A[5];
                    A[11] = A[1] * A[6] - A[2] * A[5];
                    A[12] = A[1] * A[7] - A[3] * A[5];
                    A[13] = A[9] * A[9] + A[10] * A[10] / one_M_esq + A[11] * A[11];
                    A[14] = A[8] * A[9] + A[11] * A[12];
                    A[15] = A[8] * A[8] - BIG_E * A[10] * A[10] + A[12] * A[12];

                    A15_DIV_A14 = A[14] / A[13];

                    S = A15_DIV_A14 + StrictMath.sqrt(A15_DIV_A14 * A15_DIV_A14 - A[15] / A[13]);
                }

                if(h != 0) {

                    t  = P * P * (1.0 - esq * coslat1 * coslat1) - BIG_E * one_M_esq;
                    Kp = (-u + StrictMath.sqrt(u * u - 4.0 * t * v)) / (2.0 * t);

                    lat[i] = phi;
                }
                else
                    lat[i] = StrictMath.atan2(S, StrictMath.sqrt(one_M_esq * (one_M_esq - S * S)));

                if(prj) {
                    den = A[11]*S-A[12];
                    lon[i] = normalizeLonRad(StrictMath.atan((A[8]-A[9]*S)/den));
                } else {
                    X  = a * ((B - H / Kp) * coslat1 - (Yt / Kp - D) * sinlat1);
                    Y  = a * Xt / Kp;

                    lon[i] = normalizeLonRad(lon0 + StrictMath.atan2(Y, X));
                }
            }
        }
        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double esq = ellip.getProperty("e^2");

        double h     = datum.getProperty("h"); // Usually 0
        double H     = datum.getProperty("H");
        double h0    = datum.getProperty("h0");
        double lon0  = datum.getProperty("lon0");
        double lat1  = datum.getProperty("lat1");
        double theta = datum.getProperty("theta");
        double x0    = datum.getProperty("x0");
        double y0    = datum.getProperty("y0");
        double gamma = datum.getProperty("gamma");
        double omega = datum.getProperty("omega");

        boolean tlt = datum.getProperty("useTilt")         == 1.0;
        boolean prj = datum.getProperty("useProjective")   == 1.0;
        boolean ret = datum.getProperty("returnAllPoints") == 1.0;

        double x[] = new double[lon.length];
        double y[] = new double[lat.length];

        double coslat1 = StrictMath.cos(lat1);
        double sinlat1 = StrictMath.sin(lat1);

        double coslon0 = StrictMath.cos(lon0);
        double sinlon0 = StrictMath.sin(lon0);

        double cosgamma = StrictMath.cos(gamma);
        double singamma = StrictMath.sin(gamma);

        double cosomega = StrictMath.cos(omega);
        double sinomega = StrictMath.sin(omega);

        double costheta = StrictMath.cos(theta);
        double sintheta = StrictMath.sin(theta);

        double xt, yt, A;

        if(ellip.isSphere()) {
            double P = H / a + 1.0;
            double lon_M_lon0, coslat, sinlat, c, Kp, coslon_M_lon0;

            for(int i = 0; i < lon.length; i++) {
                lon_M_lon0 = normalizeLonRad(lon[i] - lon0);

                coslon_M_lon0 = StrictMath.cos(lon_M_lon0);

                coslat = StrictMath.cos(lat[i]);
                sinlat = StrictMath.sin(lat[i]);

                c = sinlat1 * sinlat + coslat1 * coslat * coslon_M_lon0; // (5-3)

                if((ret == false && 1 / P < c) || ret == true) {
                    Kp = (P - 1.0) / (P - c);
                    x[i] = a * Kp * coslat * StrictMath.sin(lon_M_lon0);
                    y[i] = a * Kp * (coslat1 * sinlat - sinlat1 * coslat * coslon_M_lon0);

                    if(tlt) {
                        A = ((y[i] * cosgamma + x[i] * singamma) * sinomega / H) + cosomega;

                        xt = (x[i] * cosgamma - y[i] * singamma) * cosomega / A;
                        yt = (y[i] * cosgamma + x[i] * singamma) / A;

                        x[i] = xt;
                        y[i] = yt;
                    }
                }
            }
        }
        else {
            double one_M_esq = 1.0 - esq;
            double N1 = a / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);
            double P_factor = (H + N1 + h0) / a;
            double phig_dividend = N1 * esq * sinlat1 * coslat1;

            double lat1_M_phig, lon_M_lon0, coslon_M_lon0, c;
            double K, S, P, C, N, sinlat, coslat, phig, phig_old;
            int j;

            double X, Y, Z, prj_den, sinlat1_M_phig, coslat1_M_phig;

            double F_num = 0.0, V_num = 0.0, M_num = 0.0, N_num = 0.0, W_num = 0.0, T_num = 0.0, K7_num = 0.0;

            if(prj) {
                F_num = sinlat1 * sinlon0 * cosgamma - coslon0 * singamma;
                V_num = sinlat1 * sinlon0 * singamma + coslon0 * cosgamma;
                M_num = sinlat1 * coslon0 * singamma - sinlon0 * cosgamma;
                N_num = sinlat1 * coslon0 * cosgamma + sinlon0 * singamma;
                W_num = -singamma * cosomega * costheta - cosgamma * sintheta;
                T_num = -singamma * cosomega * sintheta + cosgamma * costheta;
                K7_num = coslat1 * cosgamma * sinomega - sinlat1 * cosomega;// k[6];
            }

            double U, F, V, M, W, T;
            double[] k = new double[11]; // K1..11 in text, or k0..10 in code

            for(int i = 0; i < lon.length; ++i) {
                coslat = StrictMath.cos(lat[i]);
                sinlat = StrictMath.sin(lat[i]);

                phig = lat1;
                j = 0;
                do {
                    phig_old = phig;
                    P = (coslat1 / StrictMath.cos(phig)) * P_factor;
                    phig = lat1 - StrictMath.asin(phig_dividend / (P * a));
                } while(++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phig - phig_old));

                lat1_M_phig = lat1 - phig;
                lon_M_lon0 = lon[i] - lon0;

                coslon_M_lon0 = StrictMath.cos(lon_M_lon0);
                coslat1_M_phig = StrictMath.cos(lat1_M_phig);

                N = a / StrictMath.sqrt(1.0 - esq * sinlat * sinlat);
                C = ((N + h) / a) * coslat;
                S = ((N * one_M_esq + h) / a) * sinlat;
                K = H / (P * coslat1_M_phig - S * sinlat1 - C * coslat1 * coslon_M_lon0);

                c = sinlat1 * sinlat + coslat1 * coslat * coslon_M_lon0; // (5-3)

                if((ret == false && 1 / P < c) || ret == true) {
                    sinlat1_M_phig = StrictMath.sin(lat1_M_phig);

                    x[i] = K * C * StrictMath.sin(lon_M_lon0);
                    y[i] = K * (P * sinlat1_M_phig + S * coslat1 - C * sinlat1 * coslon_M_lon0);

                    if(tlt) {
                        A = ((y[i] * cosgamma + x[i] * singamma) * sinomega / H) + cosomega;

                        xt = (x[i] * cosgamma - y[i] * singamma) * cosomega / A;
                        yt = (y[i] * cosgamma + x[i] * singamma) / A;

                        x[i] = xt;
                        y[i] = yt;
                    }
                    if(prj) {
                        U = P * ((sinlat1_M_phig) * cosgamma * sinomega + coslat1_M_phig * cosomega);
                        F = F_num / U;
                        V = V_num * cosomega / U;
                        M = M_num * cosomega / U;
                        N = N_num / U;
                        W = W_num / U;
                        T = T_num / U;

                        k[4]  = -N * sinomega - coslat1 * coslon0 * cosomega / U;
                        k[5]  = -F * sinomega - coslat1 * sinlon0 * cosomega / U;
                        k[6]  = K7_num / U;
                        k[0]  = H * (M * costheta + N * sintheta) + k[4] * x0;
                        k[1]  = H * (V * costheta + F * sintheta) + k[5] * x0;
                        k[2]  = H * W * coslat1 + k[6] * x0;
                        k[3]  = H * W * P * sinlat1_M_phig + x0;
                        k[7]  = H * (M * sintheta - N * costheta) + k[4] * y0;
                        k[8]  = H * (V * sintheta - F * costheta) + k[5] * y0;
                        k[9]  = H * T * coslat1 + k[6] * y0;
                        k[10] = H * T * P * sinlat1_M_phig + y0;

                        X = C * StrictMath.cos(lon[i]);
                        Y = C * StrictMath.sin(lon[i]);
                        Z = S;

                        prj_den = k[4] * X + k[5] * Y + k[6] * Z + 1.0;

                        x[i] = (k[0] * X + k[1] * Y + k[2] * Z + k[3]) / prj_den;
                        y[i] = (k[7] * X + k[8] * Y + k[9] * Z + k[10]) / prj_den;
                    }
                }
            }
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "h",
                "H",
                "h0",
                "lon0",
                "lat1",
                "theta",
                "x0",
                "y0",
                "gamma",
                "omega",
                "useTilt",
                "useProjective",
                "returnAllPoints"}));
    }

}