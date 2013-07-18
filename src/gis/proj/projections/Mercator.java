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