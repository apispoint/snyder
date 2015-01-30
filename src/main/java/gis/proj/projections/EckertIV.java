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
 * Pages 253-258
 *
 */

public final class EckertIV implements Pseudocylindrical, Spherical {

    private static final double X_FACTOR = 2.0 / StrictMath.sqrt(PI * (4.0 + PI));
    private static final double Y_FACTOR = 2.0 * StrictMath.sqrt(PI / (4.0 + PI));

    private static final double POS_TWO_P_HLF_PI = 2.0 + PI_DIV_2;
    private static final double NEG_TWO_P_HLF_PI = -POS_TWO_P_HLF_PI;

    public String getName() {
        return "Eckert IV";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double R    = ellip.getProperty("R");

        double lon0 = datum.getProperty("lon0");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double x_factor_R = X_FACTOR * R, y_factor_R = Y_FACTOR * R;

        double theta, cosTheta, sinTheta;

        for(int i = 0; i < lon.length; ++i) {
            theta = StrictMath.asin(y[i] / y_factor_R);
            cosTheta = StrictMath.cos(theta);
            sinTheta = StrictMath.sin(theta);

            lon[i] = normalizeLonRad(lon0 + x[i] / (x_factor_R * (1.0 + cosTheta)));
            lat[i] = StrictMath.asin((theta + sinTheta * cosTheta + 2.0 * sinTheta) / POS_TWO_P_HLF_PI);
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double R    = ellip.getProperty("R");

        double lon0 = datum.getProperty("lon0");

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        int j;
        double x_factor_R = X_FACTOR * R;
        double y_factor_R = Y_FACTOR * R;

        double theta, dTheta, sinTheta, cosTheta, sinlat;

        for(int i = 0; i < lon.length; ++i) {
            sinlat = StrictMath.sin(lat[i]);
            theta = lat[i] * 0.5;
            j = 0;
            do {
                cosTheta = StrictMath.cos(theta);
                sinTheta = StrictMath.sin(theta);
                theta +=
                        dTheta = -(theta + sinTheta * cosTheta + 2.0 * sinTheta + NEG_TWO_P_HLF_PI * sinlat) /
                        (2.0 * cosTheta * (1.0 + cosTheta));
            } while(++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(dTheta));

            x[i] = x_factor_R * normalizeLonRad(lon[i] - lon0) * (1.0 + StrictMath.cos(theta));
            y[i] = y_factor_R * StrictMath.sin(theta);
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0"}));
    }

}