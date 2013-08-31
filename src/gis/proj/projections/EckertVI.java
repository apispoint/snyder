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
 * Pages 253-258
 *
 */

public final class EckertVI implements Pseudocylindrical, Spherical {

    private static final double POS_ONE_P_HLF_PI = 1.0 + PI_DIV_2;
    private static final double NEG_ONE_P_HLF_PI = -POS_ONE_P_HLF_PI;

    private static final double SQRT_TWO_P_PI = StrictMath.sqrt(2.0 + PI);

    public String getName() {
        return "Eckert VI";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double R    = ellip.getProperty("R");

        double lon0 = datum.getProperty("lon0");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double TWO_R = 2.0 * R;

        double theta;

        for(int i = 0; i < lon.length; ++i) {
            theta = SQRT_TWO_P_PI * y[i] / TWO_R;

            lon[i] = normalizeLonRad(lon0 + SQRT_TWO_P_PI * x[i] / (R * (1.0 + StrictMath.cos(theta))));
            lat[i] = StrictMath.asin((theta + StrictMath.sin(theta)) / POS_ONE_P_HLF_PI);
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double R    = ellip.getProperty("R");

        double lon0 = datum.getProperty("lon0");

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double TWO_R = 2.0 * R;

        double theta, dTheta, sinTheta, cosTheta, sinlat;
        int j;

        for(int i = 0; i < lon.length; ++i) {
            sinlat = StrictMath.sin(lat[i]);
            theta = lat[i];
            j = 0;
            do {
                cosTheta = StrictMath.cos(theta);
                sinTheta = StrictMath.sin(theta);
                theta += dTheta = -(theta + sinTheta + NEG_ONE_P_HLF_PI * sinlat) / (1.0 + cosTheta);
            } while(++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(dTheta));

            x[i] = R * normalizeLonRad(lon[i] - lon0) * (1.0 + StrictMath.cos(theta)) / SQRT_TWO_P_PI;
            y[i] = TWO_R * theta / SQRT_TWO_P_PI;
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0"}));
    }

}