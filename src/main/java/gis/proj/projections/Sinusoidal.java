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

import gis.proj.Datum;
import gis.proj.Ellipsoid;
import gis.proj.Pseudocylindrical;

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
 * Pages 243-248
 *
 */

public final class Sinusoidal implements Pseudocylindrical {

    public String getName() {
        return "Sinusoidal";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a    = ellip.getProperty("a");
        double esq  = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double sqrt_1_M_esq = StrictMath.sqrt(1.0 - esq);
        double E1 = (1.0 - sqrt_1_M_esq) / (1.0 + sqrt_1_M_esq);

        double mu, sinlat;

        for(int i = 0; i < lon.length; ++i) {
            // Where M = y[i]
            mu = y[i] / (a * (1.0 - esq * 0.25 - 3.0 * esq * esq / 64.0 - 5.0/256.0 * esq * esq * esq));

            lat[i] = mu + ((3.0/2.0) * E1 - (27.0 / 32.0)*(E1 * E1 * E1)) * StrictMath.sin(2.0 * mu);
            lat[i] += (((21.0/16.0) * E1 * E1) - ((55.0/32.0) * E1 * E1 * E1 * E1)) * StrictMath.sin(4.0 * mu);
            lat[i] += ((151.0/96.0) * E1 * E1 * E1) * StrictMath.sin(6.0 * mu);
            lat[i] += ((1097.0/512.0) * E1 * E1 * E1 * E1) * StrictMath.sin(8.0 * mu);

            if(NEAR_ZERO_RAD < StrictMath.abs(PI_DIV_2 - lat[i])) {
                sinlat = StrictMath.sin(lat[i]);
                lon[i] = lon0 + x[i] * StrictMath.sqrt(1.0 - esq * sinlat * sinlat) / (a * StrictMath.cos(lat[i]));
            }
            else
                lon[i] = lon0;

            lon[i] = normalizeLonRad(lon[i]);
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a    = ellip.getProperty("a");
        double esq  = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
        double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
        double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
        double M_factor6 = (esq * esq * esq) * 0.011393229167;

        double sinlat, M;

        for(int i = 0; i < lon.length; ++i) {
            // Series
            M  = lat[i] * M_factor0;
            M -= StrictMath.sin(2.0 * lat[i]) * M_factor2;
            M += StrictMath.sin(4.0 * lat[i]) * M_factor4;
            M -= StrictMath.sin(6.0 * lat[i]) * M_factor6;
            M *= a;

            sinlat = StrictMath.sin(lat[i]);

            x[i] = a * normalizeLonRad(lon[i] - lon0) * StrictMath.cos(lat[i]) / StrictMath.sqrt(1.0 - esq * sinlat * sinlat);
            y[i] = M;
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0"}));
    }

}