/*
Snyder Projection Implementation

Copyright (c) 2012-2015, APIS Point, LLC

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
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
 * Pages 191-202
 *
 */

public final class AzimuthalEquidistant implements Azimuthal {

    public String getName() {
        return "Azimuthal Equidistant";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double R   = ellip.getProperty("R");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat1 = datum.getProperty("lat1");
        double x0   = datum.getProperty("x0");
        double y0   = datum.getProperty("y0");

        boolean useGuam       = datum.getProperty("useGuam")       == 1.0;
        boolean useMicronesia = datum.getProperty("useMicronesia") == 1.0;

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        // Polar aspect
        if (POLE_RAD - StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
            double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
            double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
            double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
            double M_factor6 = (esq * esq * esq) * 0.011393229167;

            double poleFactor = lat1 < 0.0 ? -1 : 1.0;

            double sqrt_1_M_esq = StrictMath.sqrt(1.0 - esq);
            double E1 = (1.0 - sqrt_1_M_esq) / (1.0 + sqrt_1_M_esq);
            double MU_DIV = (a * (1.0 - esq * 0.25 - esq * esq * (3.0 / 64.0) - esq * esq * esq * (5.0 / 256.0)));

            double Mp  = PI_DIV_2 * M_factor0;
            Mp -= StrictMath.sin(2.0 * PI_DIV_2) * M_factor2;
            Mp += StrictMath.sin(4.0 * PI_DIV_2) * M_factor4;
            Mp -= StrictMath.sin(6.0 * PI_DIV_2) * M_factor6;
            Mp *= a;

            double M, mu, rho;

            for(int i = 0; i < x.length; ++i) {
                rho = StrictMath.hypot(x[i], y[i]);
                M = -rho * poleFactor + Mp * poleFactor;
                mu = M / MU_DIV;

                lat[i] = mu + ((3.0/2.0) * E1 - (27.0 / 32.0)*(E1 * E1 * E1)) * StrictMath.sin(2.0 * mu);
                lat[i] += (((21.0/16.0) * E1 * E1) - ((55.0/32.0) * E1 * E1 * E1 * E1)) * StrictMath.sin(4.0 * mu);
                lat[i] += ((151.0/96.0) * E1 * E1 * E1) * StrictMath.sin(6.0 * mu);
                lat[i] += ((1097.0/512.0) * E1 * E1 * E1 * E1) * StrictMath.sin(8.0 * mu);

                lon[i] = normalizeLonRad(lon0 + StrictMath.atan2(x[i], -y[i] * poleFactor));
            }
        }
        else if(useGuam == true) {
            double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
            double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
            double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
            double M_factor6 = (esq * esq * esq) * 0.011393229167;

            double M1  = lat1 * M_factor0;
            M1 -= StrictMath.sin(2.0 * lat1) * M_factor2;
            M1 += StrictMath.sin(4.0 * lat1) * M_factor4;
            M1 -= StrictMath.sin(6.0 * lat1) * M_factor6;
            M1 *= a;

            double sqrt_1_M_esq = StrictMath.sqrt(1.0 - esq);
            double E1 = (1.0 - sqrt_1_M_esq) / (1.0 + sqrt_1_M_esq);

            double MU_DIV = (a * (1.0 - esq * 0.25 - esq * esq * (3.0 / 64.0) - esq * esq * esq * (5.0 / 256.0)));

            double TWO_a = 2.0 * a;

            double M, xtmp, ytmp, phi, phi_old, xtmp2, sinphi, mu;
            int j;

            for(int i = 0; i < lon.length; ++i) {
                xtmp = x[i] - x0;
                ytmp = y[i] - y0;

                xtmp2 = xtmp * xtmp;

                phi = lat1;
                j = 0;
                do {
                    phi_old = phi;
                    sinphi = StrictMath.sin(phi);

                    M = M1 + ytmp - xtmp2 * StrictMath.tan(phi) * StrictMath.sqrt(1.0 - esq * sinphi * sinphi) / TWO_a;
                    mu = M / MU_DIV;

                    phi = mu + ((3.0/2.0) * E1 - (27.0 / 32.0)*(E1 * E1 * E1)) * StrictMath.sin(2.0 * mu);
                    phi += (((21.0/16.0) * E1 * E1) - ((55.0/32.0) * E1 * E1 * E1 * E1)) * StrictMath.sin(4.0 * mu);
                    phi += ((151.0/96.0) * E1 * E1 * E1) * StrictMath.sin(6.0 * mu);
                    phi += ((1097.0/512.0) * E1 * E1 * E1 * E1) * StrictMath.sin(8.0 * mu);
                } while(++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phi - phi_old));
                lat[i] = phi;
                sinphi = StrictMath.sin(phi);
                lon[i] = normalizeLonRad(lon0 + xtmp * StrictMath.sqrt(1.0 - esq * sinphi * sinphi) / (a * StrictMath.cos(phi)));
            }
        }
        else if(useMicronesia == true) {
            double sinlat1 = StrictMath.sin(lat1);
            double coslat1 = StrictMath.cos(lat1);
            double N1 = a / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);
            double one_M_esq = (1.0 - esq);
            double M_esq_coslat1coslat1 = -esq * coslat1 * coslat1;

            double psi, Az, c, cosAz, A, B, D, E, F, x_M_x0, y_M_y0;

            for(int i = 0; i < x.length; ++i) {
                x_M_x0 = x[i] - x0;
                y_M_y0 = y[i] - y0;

                c = StrictMath.hypot(x_M_x0, y_M_y0);
                Az = StrictMath.atan2(x_M_x0, y_M_y0);
                cosAz = StrictMath.cos(Az);
                A = M_esq_coslat1coslat1 * cosAz * cosAz / one_M_esq;
                B = 3.0 * esq * (1.0 - A) * sinlat1 * coslat1 * cosAz / one_M_esq;
                D = c/N1;
                E = D - A * (1.0 + A) * (D*D*D) / 6.0 - B * (1.0 + 3 *A)*(D*D*D*D)/24;
                F = 1.0 - A*E*E*0.5-(B*E*E*E)/6;
                psi = StrictMath.asin(sinlat1 * StrictMath.cos(E) + coslat1 * StrictMath.sin(E) * cosAz);

                lon[i] = normalizeLonRad(lon0 + StrictMath.asin(StrictMath.sin(Az) * StrictMath.sin(E) / StrictMath.cos(psi)));
                lat[i] = StrictMath.atan((1-esq*F*sinlat1/StrictMath.sin(psi))*StrictMath.tan(psi)/one_M_esq);
            }
        }
        // Oblique et Equatorial aspects
        else {

            double
            coslat1 = StrictMath.cos(lat1),
            sinlat1 = StrictMath.sin(lat1);

            double rho, c, cosc, sinc;

            for(int i = 0; i < x.length; ++i) {
                rho = StrictMath.hypot(x[i], y[i]);
                c = rho / R;

                cosc = StrictMath.cos(c);
                sinc = StrictMath.sin(c);

                if(StrictMath.abs(rho) < NEAR_ZERO_DEG) {
                    lon[i] = lon0;
                    lat[i] = lat1;
                } else {
                    lon[i] = normalizeLonRad(lon0 + StrictMath.atan2(
                            x[i] * sinc,
                            rho * coslat1 *cosc - y[i] * sinlat1 * sinc
                            ));

                    lat[i] = StrictMath.asin(cosc * sinlat1 + (y[i] * sinc * coslat1 / rho));
                }
            }
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double R   = ellip.getProperty("R");
        double e   = ellip.getProperty("e");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat1 = datum.getProperty("lat1");
        double x0   = datum.getProperty("x0");
        double y0   = datum.getProperty("y0");

        boolean useGuam       = datum.getProperty("useGuam")       == 1.0;
        boolean useMicronesia = datum.getProperty("useMicronesia") == 1.0;

        double x[] = new double[lon.length];
        double y[] = new double[lat.length];

        // Equatorial aspect
        if(StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
            double lon_M_lon0, coslat, kp, c, coslon_M_lon0;

            for(int i = 0; i < lon.length; ++i) {
                coslat = StrictMath.cos(lat[i]);

                lon_M_lon0 = normalizeLonRad(lon[i] - lon0);
                coslon_M_lon0 = StrictMath.cos(lon_M_lon0);

                // eq5_3
                c = StrictMath.acos(coslat * coslon_M_lon0);
                kp = c / StrictMath.sin(c);

                if(Double.isInfinite(kp) || Double.isNaN(kp)) {
                    kp = 1.0;
                    x[i] = 0.0;
                    y[i] = 0.0;
                } else {
                    x[i] = R * kp * coslat * StrictMath.sin(lon_M_lon0);
                    y[i] = R * kp * StrictMath.sin(lat[i]);
                }
            }
        }
        // Polar aspect
        else if (POLE_RAD - StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
            double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
            double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
            double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
            double M_factor6 = (esq * esq * esq) * 0.011393229167;

            double poleFactor = lat1 < 0.0 ? -1.0 : 1.0;

            //
            //TODO
            //
            // Activate Maths elliptic curve approximation of latitude since it
            // used in Bonne, Albers, etc. Also, generalize the inverse elliptic curve
            // approximation of latitude.
            double Mp  = PI_DIV_2 * M_factor0;
            Mp -= StrictMath.sin(2.0 * PI_DIV_2) * M_factor2;
            Mp += StrictMath.sin(4.0 * PI_DIV_2) * M_factor4;
            Mp -= StrictMath.sin(6.0 * PI_DIV_2) * M_factor6;
            Mp *= a;

            double M, rho;

            for(int i = 0; i < lon.length; ++i) {
                M  = lat[i] * M_factor0;
                M -= StrictMath.sin(2.0 * lat[i]) * M_factor2;
                M += StrictMath.sin(4.0 * lat[i]) * M_factor4;
                M -= StrictMath.sin(6.0 * lat[i]) * M_factor6;
                M *= a;

                rho = Mp + (-M * poleFactor);

                x[i] = rho * StrictMath.sin(lon[i] - lon0);
                y[i] = (-rho * poleFactor) * StrictMath.cos(lon[i] - lon0);
            }
        }
        else if(useGuam == true) {
            double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
            double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
            double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
            double M_factor6 = (esq * esq * esq) * 0.011393229167;

            double M1  = lat1 * M_factor0;
            M1 -= StrictMath.sin(2.0 * lat1) * M_factor2;
            M1 += StrictMath.sin(4.0 * lat1) * M_factor4;
            M1 -= StrictMath.sin(6.0 * lat1) * M_factor6;
            M1 *= a;

            double TWO_a = 2.0 * a;

            double sinlat, sqrt_1_M_esqsinlatsinlat, M;

            for(int i = 0; i < lon.length; ++i) {
                M  = lat[i] * M_factor0;
                M -= StrictMath.sin(2.0 * lat[i]) * M_factor2;
                M += StrictMath.sin(4.0 * lat[i]) * M_factor4;
                M -= StrictMath.sin(6.0 * lat[i]) * M_factor6;
                M *= a;

                sinlat = StrictMath.sin(lat[i]);
                sqrt_1_M_esqsinlatsinlat = StrictMath.sqrt(1.0 - esq * sinlat * sinlat);

                x[i] = a * normalizeLonRad(lon[i] - lon0) * StrictMath.cos(lat[i]) / sqrt_1_M_esqsinlatsinlat;
                y[i] = M - M1 + x[i] * x[i] * StrictMath.tan(lat[i]) * sqrt_1_M_esqsinlatsinlat / TWO_a;
                x[i] += x0;
                y[i] += y0;
            }
        }
        else if(useMicronesia == true) {
            double sinlat1 = StrictMath.sin(lat1);
            double coslat1 = StrictMath.cos(lat1);
            double N1 = a / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);
            double one_M_esq = (1.0 - esq);
            double sqrt_one_M_esq = StrictMath.sqrt(1.0 - esq);

            double sinlat, sqrt_1_M_esqsinlatsinlat, N, psi, Az, s, G, H, c;

            for(int i = 0; i < lon.length; ++i) {
                sinlat = StrictMath.sin(lat[i]);
                sqrt_1_M_esqsinlatsinlat = StrictMath.sqrt(1.0 - esq * sinlat * sinlat);

                N = a / sqrt_1_M_esqsinlatsinlat;
                psi = StrictMath.atan(
                        one_M_esq * StrictMath.tan(lat[i]) + esq * N1 * sinlat1 / (N * StrictMath.cos(lat[i]))
                        );
                Az = StrictMath.atan2(
                        StrictMath.sin(normalizeLonRad(lon[i] - lon0)),
                        StrictMath.cos(lat1) * StrictMath.tan(psi) - sinlat1 * StrictMath.cos(normalizeLonRad(lon[i] - lon0))
                        );
                s = StrictMath.asin(StrictMath.sin(normalizeLonRad(lon[i] - lon0)) * StrictMath.cos(psi) / StrictMath.sin(Az));
                G = e * sinlat1 / sqrt_one_M_esq;
                H = e * coslat1 * StrictMath.cos(Az) / sqrt_one_M_esq;
                c = N1 * s * (
                        1.0 - s *s * H * H * (1.0 -H * H) /6.0 + ((s*s*s)/8.0)*G*H*(1.0 - 2.0 *H*H) +
                        ((s*s*s*s)/120.0) * (H* H*(4.0-7*H*H)-3.0*G*G*(1.0 -7.0*H*H))-((s*s*s*s*s)/48.0) *G*H
                        );

                x[i] = c * StrictMath.sin(Az) + x0;
                y[i] = c * StrictMath.cos(Az) + y0;
            }
        }
        // Oblique aspect
        else {
            double coslat1 = StrictMath.cos(lat1);
            double sinlat1 = StrictMath.sin(lat1);

            double lon_M_lon0, coslat, sinlat, kp, c, coslon_M_lon0, rhc;

            for(int i = 0; i < lon.length; ++i) {
                coslat = StrictMath.cos(lat[i]);
                sinlat = StrictMath.sin(lat[i]);

                lon_M_lon0 = normalizeLonRad(lon[i] - lon0);
                coslon_M_lon0 = StrictMath.cos(lon_M_lon0);

                // eq5_3
                rhc = sinlat1 * sinlat + coslat1 * coslat * coslon_M_lon0;
                c = StrictMath.acos(rhc);
                if(StrictMath.abs(1.0 - rhc) < NEAR_ZERO_DEG) {
                    y[i] = x[i] = 0;
                }
                else {
                    kp = c / StrictMath.sin(c);

                    x[i] = R * kp * coslat * StrictMath.sin(lon_M_lon0);
                    y[i] = R * kp * (coslat1 * sinlat - sinlat1 * coslat * coslon_M_lon0);
                }
            }
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0",
                "lat1",
                "useGuam",
                "useMicronesia",
                "x0",
                "y0"}));
    }

}