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

import gis.proj.Cylindrical;
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
 * Pages 48-65
 *
 */

public final class TransverseMercator implements Cylindrical {

    private static final double ONE_6TH    = 1.0 / 6.0;
    private static final double ONE_24TH   = 1.0 / 24.0;
    private static final double ONE_120TH  = 1.0 / 120.0;
    private static final double ONE_720TH  = 1.0 / 720.0;
    private static final double THR_SCNDTH = 3.0 / 2.0;

    public String getName() {
        return "Transverse Mercator";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a    = ellip.getProperty("a");
        double esq  = ellip.getProperty("e^2");
        double epsq = ellip.getProperty("e'^2");

        double k0   = datum.getProperty("k0");
        double lon0 = datum.getProperty("lon0");
        double lat0 = datum.getProperty("lat0");

        double epsq8       = 8.0 * epsq;
        double neg_epsq9   = -9.0 * epsq;
        double neg_epsq252 = -252.0 * epsq;

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double E1     = (1.0 - StrictMath.sqrt(1.0 - esq)) / (1.0 + StrictMath.sqrt(1.0 - esq));
        double MU_DIV = (a * (1.0 - esq * 0.25 - esq * esq * (3.0 / 64.0) - esq * esq * esq * (5.0 / 256.0)));

        double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
        double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
        double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
        double M_factor6 = (esq * esq * esq) * 0.011393229167;

        // Series
        double
        M0  = lat0 * M_factor0;
        M0 -= StrictMath.sin(2.0 * lat0) * M_factor2;
        M0 += StrictMath.sin(4.0 * lat0) * M_factor4;
        M0 -= StrictMath.sin(6.0 * lat0) * M_factor6;
        M0 *= a;

        double tanlat1, coslat1, sinlat1;
        double C1, T1, N1, R1, D, M, mu, lat1, TWOT1, T1T1, DD, DDD, C1C1;

        for(int i = 0; i < lon.length; ++i) {
            M = M0 + y[i] / k0;
            mu = M / MU_DIV;

            lat1  = mu + ((3.0/2.0) * E1 - (27.0 / 32.0)*(E1 * E1 * E1)) * StrictMath.sin(2.0 * mu);
            lat1 += (((21.0/16.0) * E1 * E1) - ((55.0/32.0) * E1 * E1 * E1 * E1)) * StrictMath.sin(4.0 * mu);
            lat1 += ((151.0/96.0) * E1 * E1 * E1) * StrictMath.sin(6.0 * mu);
            lat1 += ((1097.0/512.0) * E1 * E1 * E1 * E1) * StrictMath.sin(8.0 * mu);

            tanlat1 = StrictMath.tan(lat1);
            coslat1 = StrictMath.cos(lat1);
            sinlat1 = StrictMath.sin(lat1);

            C1 = epsq * coslat1 * coslat1;
            T1 = tanlat1 * tanlat1;
            N1 = a / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);
            R1 = a * (1.0 - esq) / StrictMath.pow(1.0 - esq * sinlat1 * sinlat1, THR_SCNDTH);
            D = x[i] / (N1 * k0);

            C1C1 = C1 * C1;
            T1T1 = T1 * T1;
            TWOT1 = 2.0 * T1;
            DD = D * D;

            if(NEAR_ZERO_RAD <= POLE_RAD - StrictMath.abs(lat[i])) {
                lon[i] = lon0 + (D - (1.0 + TWOT1 + C1) * ((DDD = DD * D) * ONE_6TH) +
                        (5.0 - 2.0 * C1 + 28.0 * T1 - 3.0 * C1C1 + epsq8 + 24.0 * T1T1) *
                        (DDD * DD * ONE_120TH)) / coslat1;

                lat[i] = lat1 - (N1 * tanlat1 / R1) * (
                        DD * 0.5 - (5.0 + 3.0 * T1 + 10.0 * C1 - 4.0 * C1C1 + neg_epsq9) *
                        (DDD * D * ONE_24TH) + (61.0 + 90.0 * T1 + 298.0 * C1 + 45.0 * T1T1 +
                                neg_epsq252 - 3.0 * C1C1) * (DDD * DDD * ONE_720TH)
                        );
            }
            else {
                x[i] = lon0;
                y[i] = y[i] < 0 ? N_PI_DIV_2 : PI_DIV_2;
            }
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a    = ellip.getProperty("a");
        double esq  = ellip.getProperty("e^2");
        double epsq = ellip.getProperty("e'^2");

        double k0   = datum.getProperty("k0");
        double lon0 = datum.getProperty("lon0");
        double lat0 = datum.getProperty("lat0");

        double neg_esq     = -esq;
        double neg_epsq58  = -58.0 * epsq;
        double neg_epsq330 = -330.0 * epsq;

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double N, T, C, A, M, M0;

        double sinlat, tanlat, coslat;

        double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
        double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
        double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
        double M_factor6 = (esq * esq * esq) * 0.011393229167;

        // Series
        M0  = lat0 * M_factor0;
        M0 -= StrictMath.sin(2.0 * lat0) * M_factor2;
        M0 += StrictMath.sin(4.0 * lat0) * M_factor4;
        M0 -= StrictMath.sin(6.0 * lat0) * M_factor6;
        M0 *= a;

        double TT, AA, AAA;

        for(int i = 0; i < lon.length; ++i) {
            sinlat = StrictMath.sin(lat[i]);
            tanlat = StrictMath.tan(lat[i]);
            coslat = StrictMath.cos(lat[i]);

            N = a / StrictMath.sqrt(1.0 + neg_esq * sinlat * sinlat);
            T = tanlat * tanlat;
            C = epsq * coslat * coslat;
            A = normalizeLonRad(lon[i] - lon0) * coslat;

            // Series
            M  = lat[i] * M_factor0;
            M -= StrictMath.sin(2.0 * lat[i]) * M_factor2;
            M += StrictMath.sin(4.0 * lat[i]) * M_factor4;
            M -= StrictMath.sin(6.0 * lat[i]) * M_factor6;
            M *= a;

            TT = T * T;
            AA = A * A;

            if(NEAR_ZERO_RAD <= POLE_RAD - StrictMath.abs(lat[i])) {
                x[i] = k0 * N * (
                        A + (1.0 - T + C) * ((AAA = AA * A) * ONE_6TH) +
                        (5.0 - 18.0 * T + TT + 72.0 * C + neg_epsq58) * (AAA * AA * ONE_120TH)
                        );
                y[i] = k0 * (
                        M - M0 + N * tanlat * (AA * 0.5 + (5.0 - T + 9.0 * C + 4.0 * C * C) *
                        (AAA * A * ONE_24TH) + (61.0 - 58.0 * T + TT + 600.0 * C + neg_epsq330) *
                        (AAA * AAA * ONE_720TH)) 
                        );
            }
            else {
                x[i] = 0;
                y[i] = k0 * (M - M0);
            }
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "k0",
                "lon0",
                "lat0"}));
    }

}