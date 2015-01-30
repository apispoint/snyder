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
 * Pages 111-115
 *
 */

public final class EquidistantConic implements Conic {

    public String getName() {
        return "Equidistant Conic";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat0 = datum.getProperty("lat0");
        double lat1 = datum.getProperty("lat1");
        double lat2 = datum.getProperty("lat2");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double
               sinlat1 = StrictMath.sin(lat1),
               sinlat2 = StrictMath.sin(lat2);

        double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
        double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
        double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
        double M_factor6 = (esq * esq * esq) * 0.011393229167;

        //
        //TODO
        //
        // Activate Maths elliptic curve approximation of latitude since it
        // used in Bonne, Albers, etc. Also, generalize the inverse elliptic curve
        // approximation of latitude.
        double M0  = lat0 * M_factor0;
        M0 -= StrictMath.sin(2.0 * lat0) * M_factor2;
        M0 += StrictMath.sin(4.0 * lat0) * M_factor4;
        M0 -= StrictMath.sin(6.0 * lat0) * M_factor6;
        M0 *= a;

        double M1  = lat1 * M_factor0;
        M1 -= StrictMath.sin(2.0 * lat1) * M_factor2;
        M1 += StrictMath.sin(4.0 * lat1) * M_factor4;
        M1 -= StrictMath.sin(6.0 * lat1) * M_factor6;
        M1 *= a;

        double M2  = lat2 * M_factor0;
        M2 -= StrictMath.sin(2.0 * lat2) * M_factor2;
        M2 += StrictMath.sin(4.0 * lat2) * M_factor4;
        M2 -= StrictMath.sin(6.0 * lat2) * M_factor6;
        M2 *= a;

        double m1 = StrictMath.cos(lat1) / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);
        double m2 = StrictMath.cos(lat2) / StrictMath.sqrt(1.0 - esq * sinlat2 * sinlat2);

        double nu = sinlat1;

        if(NEAR_ZERO_RAD < StrictMath.abs(lat1 - lat2))
            nu = a * (m1 - m2) / (M2 - M1);

        double nuFactor = nu < 0.0 ? -1.0 : 1.0;

        double G = m1 / nu + M1 / a;
        double rho0 = a * G - M0;

        double sqrt_1_M_esq = StrictMath.sqrt(1.0 - esq);
        double E1 = (1.0 - sqrt_1_M_esq) / (1.0 + sqrt_1_M_esq);

        double M, mu, rho, theta;

        for(int i = 0; i < lon.length; ++i) {
            theta = StrictMath.atan2(x[i] * nuFactor, rho0 * nuFactor - y[i] * nuFactor);
            lon[i] = normalizeLonRad(lon0 + theta / nu);

            rho = StrictMath.copySign(StrictMath.hypot(x[i], rho0 - y[i]), nu);

            M = a * G - rho;
            mu = M / (a * (1.0 - esq * 0.25 - 3.0 * esq * esq / 64.0 - 5.0/256.0 * esq * esq * esq));

            lat[i] = mu + ((3.0/2.0) * E1 - (27.0 / 32.0)*(E1 * E1 * E1)) * StrictMath.sin(2.0 * mu);
            lat[i] += (((21.0/16.0) * E1 * E1) - ((55.0/32.0) * E1 * E1 * E1 * E1)) * StrictMath.sin(4.0 * mu);
            lat[i] += ((151.0/96.0) * E1 * E1 * E1) * StrictMath.sin(6.0 * mu);
            lat[i] += ((1097.0/512.0) * E1 * E1 * E1 * E1) * StrictMath.sin(8.0 * mu);
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat0 = datum.getProperty("lat0");
        double lat1 = datum.getProperty("lat1");
        double lat2 = datum.getProperty("lat2");

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double rho, theta;

        double
        sinlat1 = StrictMath.sin(lat1),
        sinlat2 = StrictMath.sin(lat2);

        double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
        double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
        double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
        double M_factor6 = (esq * esq * esq) * 0.011393229167;

        //
        //TODO
        //
        // Activate Maths elliptic curve approximation of latitude since it
        // used in Bonne, Albers, etc. Also, generalize the inverse elliptic curve
        // approximation of latitude.
        double M0  = lat0 * M_factor0;
        M0 -= StrictMath.sin(2.0 * lat0) * M_factor2;
        M0 += StrictMath.sin(4.0 * lat0) * M_factor4;
        M0 -= StrictMath.sin(6.0 * lat0) * M_factor6;
        M0 *= a;

        double M1  = lat1 * M_factor0;
        M1 -= StrictMath.sin(2.0 * lat1) * M_factor2;
        M1 += StrictMath.sin(4.0 * lat1) * M_factor4;
        M1 -= StrictMath.sin(6.0 * lat1) * M_factor6;
        M1 *= a;

        double M2  = lat2 * M_factor0;
        M2 -= StrictMath.sin(2.0 * lat2) * M_factor2;
        M2 += StrictMath.sin(4.0 * lat2) * M_factor4;
        M2 -= StrictMath.sin(6.0 * lat2) * M_factor6;
        M2 *= a;

        double m1 = StrictMath.cos(lat1) / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);
        double m2 = StrictMath.cos(lat2) / StrictMath.sqrt(1.0 - esq * sinlat2 * sinlat2);

        double nu = sinlat1;

        if(NEAR_ZERO_RAD < StrictMath.abs(lat1 - lat2))
            nu = a * (m1 - m2) / (M2 - M1);

        double G = m1 / nu + M1 / a;
        double rho0 = a * G - M0;

        double M;

        for(int i = 0; i < lon.length; ++i) {
            M  = lat[i] * M_factor0;
            M -= StrictMath.sin(2.0 * lat[i]) * M_factor2;
            M += StrictMath.sin(4.0 * lat[i]) * M_factor4;
            M -= StrictMath.sin(6.0 * lat[i]) * M_factor6;
            M *= a;

            rho = a * G - M;
            theta = nu * normalizeLonRad(lon[i] - lon0);

            x[i] = rho * StrictMath.sin(theta);
            y[i] = rho0 - rho * StrictMath.cos(theta);
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0",
                "lat0",
                "lat1",
                "lat2"}));
    }

}