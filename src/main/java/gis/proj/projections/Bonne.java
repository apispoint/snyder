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

import gis.proj.Conic;
import gis.proj.Datum;
import gis.proj.Ellipsoid;
import gis.proj.SnyderMath;

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
 * Pages 138-140
 *
 */

public final class Bonne implements Conic {

    public String getName() {
        return "Bonne";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a    = ellip.getProperty("a");
        double esq  = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat1 = datum.getProperty("lat1");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
        double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
        double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
        double M_factor6 = (esq * esq * esq) * 0.011393229167;

        double M1, m1;

        // Avoid division by Zero
        // TODO When lat1 == 0; use Sinusoidal projection
        //        if(StrictMath.abs(lat1) < NEAR_ZERO_RAD)
        //            lat1 = NEAR_ZERO_RAD;

        M1  = lat1 * M_factor0;
        M1 -= StrictMath.sin(2.0 * lat1) * M_factor2;
        M1 += StrictMath.sin(4.0 * lat1) * M_factor4;
        M1 -= StrictMath.sin(6.0 * lat1) * M_factor6;
        M1 *= a;

        double sinlat1 = StrictMath.sin(lat1);
        m1 = StrictMath.cos(lat1) / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);

        double rhoFactor = (lat1 < 0 ? -1.0 : 1.0), rho, M, mu, sinlat, m;

        double sqrt_1_M_esq = StrictMath.sqrt(1.0 - esq);
        double E1 = (1.0 - sqrt_1_M_esq) / (1.0 + sqrt_1_M_esq);
        double MU_DIV = (a * (1.0 - esq * 0.25 - 3.0 * esq * esq / 64.0 - 5.0/256.0 * esq * esq * esq));

        for(int i = 0; i < lon.length; ++i) {
            rho = rhoFactor * StrictMath.sqrt(x[i] * x[i] + StrictMath.pow(a * m1 / sinlat1 - y[i], 2.0));
            M = a * m1 / sinlat1 + M1 - rho;
            mu = M / MU_DIV;

            // (3-21)
            lat[i] = mu + ((3.0/2.0) * E1 - (27.0 / 32.0)*(E1 * E1 * E1)) * StrictMath.sin(2.0 * mu);
            lat[i] += (((21.0/16.0) * E1 * E1) - ((55.0/32.0) * E1 * E1 * E1 * E1)) * StrictMath.sin(4.0 * mu);
            lat[i] += ((151.0/96.0) * E1 * E1 * E1) * StrictMath.sin(6.0 * mu);
            lat[i] += ((1097.0/512.0) * E1 * E1 * E1 * E1) * StrictMath.sin(8.0 * mu);

            sinlat = StrictMath.sin(lat[i]);
            m = StrictMath.cos(lat[i]) / StrictMath.sqrt(1.0 - esq * sinlat * sinlat);

            if(SnyderMath.POLE_RAD - StrictMath.abs(lat[i]) < NEAR_ZERO_RAD)
                lon[i] = lon0;
            else
                lon[i] = normalizeLonRad(lon0 + rho * (StrictMath.atan2(rhoFactor * x[i], rhoFactor * (a * m1 / sinlat1 - y[i]))) / (a * m));
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a    = ellip.getProperty("a");
        double esq  = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat1 = datum.getProperty("lat1");

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
        double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
        double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
        double M_factor6 = (esq * esq * esq) * 0.011393229167;

        double M1, m1;

        // Avoid division by Zero
        // TODO When lat1 == 0; use Sinusoidal projection
        //        if(StrictMath.abs(lat1) < NEAR_ZERO_RAD)
        //            lat1 = NEAR_ZERO_RAD;

        M1  = lat1 * M_factor0;
        M1 -= StrictMath.sin(2.0 * lat1) * M_factor2;
        M1 += StrictMath.sin(4.0 * lat1) * M_factor4;
        M1 -= StrictMath.sin(6.0 * lat1) * M_factor6;
        M1 *= a;

        double sinlat1 = StrictMath.sin(lat1);
        m1 = StrictMath.cos(lat1) / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);

        double sinlat, m, M, rho, EE;

        for(int i = 0; i < lon.length; ++i) {
            sinlat = StrictMath.sin(lat[i]);

            m = StrictMath.cos(lat[i]) / StrictMath.sqrt(1.0 - esq * sinlat * sinlat);

            // (3-21)
            M  = lat[i] * M_factor0;
            M -= StrictMath.sin(2.0 * lat[i]) * M_factor2;
            M += StrictMath.sin(4.0 * lat[i]) * M_factor4;
            M -= StrictMath.sin(6.0 * lat[i]) * M_factor6;
            M *= a;

            rho = a * m1 / sinlat1 + M1 - M;
            EE = a * m * normalizeLonRad(lon[i] - lon0) / rho;

            x[i] = rho * StrictMath.sin(EE);
            y[i] = a * m1 / sinlat1 - rho * StrictMath.cos(EE);
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0",
                "lat1"}));
    }

}