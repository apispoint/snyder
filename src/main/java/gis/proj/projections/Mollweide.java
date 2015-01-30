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
import gis.proj.Spherical;

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
 * Pages 249-252
 *
 */

public final class Mollweide implements Pseudocylindrical, Spherical {

    public String getName() {
        return "Mollweide";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double R    = ellip.getProperty("R");

        double lon0 = datum.getProperty("lon0");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double phi, phi2;

        for(int i = 0; i < x.length; ++i) {
            phi  = StrictMath.asin(y[i] / (SQRT_2 * R));
            phi2 = 2.0 * phi;

            lon[i] = normalizeLonRad(lon0 + (PI * x[i]) / (SQRT_8 * R * StrictMath.cos(phi)));
            lat[i] = StrictMath.asin((phi2 + StrictMath.sin(phi2)) * _1_DIV_PI);
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double R    = ellip.getProperty("R");

        double lon0 = datum.getProperty("lon0");

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double theta, theta_old, sinLatPI;
        int j;

        for(int i = 0; i < lon.length; ++i) {
            sinLatPI = StrictMath.sin(lat[i]) * PI;

            j = 0;
            theta = lat[i];
            do {
                theta_old = theta;
                theta -= (theta + StrictMath.sin(theta) - sinLatPI) / (1.0 + StrictMath.cos(theta));
            } while(++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(theta - theta_old));

            theta *= 0.5;

            x[i] = R * SQRT_8_DIV_PI * normalizeLonRad(lon[i] - lon0) * StrictMath.cos(theta);
            y[i] = R * SQRT_2 * StrictMath.sin(theta);
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0"}));
    }

}