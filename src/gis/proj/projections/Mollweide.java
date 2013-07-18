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