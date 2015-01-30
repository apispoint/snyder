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
 * Pages 38-47
 *
 */

public final class Mercator implements Cylindrical {

    public String getName() {
        return "Mercator";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a    = ellip.getProperty("a");
        double e    = ellip.getProperty("e");

        double lon0 = datum.getProperty("lon0");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double hlf_e = e * 0.5;

        double esinlat, t, phi, phi_old;
        int j;

        for(int i = 0; i < lon.length; ++i) {
            lon[i] = normalizeLonRad(x[i] / a + lon0);

            t = StrictMath.pow(_e_, -y[i] / a);

            phi = PI_DIV_2 - 2.0 * StrictMath.atan(t);
            j = 0;
            do {
                phi_old = phi;
                esinlat = e * StrictMath.sin(phi);
                phi = PI_DIV_2 - 2.0 * StrictMath.atan(t * StrictMath.pow((1.0 - esinlat) / (1.0 + esinlat), hlf_e));
            } while (++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phi - phi_old));

            lat[i] = phi;
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a    = ellip.getProperty("a");
        double e    = ellip.getProperty("e");

        double lon0 = datum.getProperty("lon0");

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double hlf_a = a * 0.5;

        double sinlat, esinlat;

        for(int i = 0; i < lon.length; ++i) {
            x[i] = a * normalizeLonRad(lon[i] - lon0);

            sinlat = StrictMath.sin(lat[i]);
            esinlat = e * sinlat;
            y[i] = hlf_a * StrictMath.log(((1.0 + sinlat) / (1.0 - sinlat)) * StrictMath.pow((1.0 - esinlat) / (1.0 + esinlat), e));
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0"}));
    }

}