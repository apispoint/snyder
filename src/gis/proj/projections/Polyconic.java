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
 * Pages 124-137
 *
 */

public final class Polyconic implements Conic {

    private static final double CUTOFF_LAT_DIST = PI_DIV_2;

    public String getName() {
        return "Polyconic";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double asq = ellip.getProperty("a^2");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat0 = datum.getProperty("lat0");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

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

        double A, B, C, phi, phi_old, Mn, Mpn, sinphi, sin2phi;
        int j;

        for(int i = 0; i < lon.length; ++i) {
            if(StrictMath.abs(y[i] - M0) < NEAR_ZERO_DEG) {
                lon[i] = normalizeLonRad(x[i] / a + lon0);
                lat[i] = 0.0;
            }
            else {
                A = (M0 + y[i]) / a;
                B = (x[i] * x[i]) / asq + (A * A);

                phi = A;
                j = 0;
                do {
                    phi_old = phi;
                    sinphi = StrictMath.sin(phi);
                    sin2phi = StrictMath.sin(2.0 * phi);

                    C = StrictMath.sqrt(1.0 - esq * sinphi * sinphi) * StrictMath.tan(phi);

                    // Ma is not required because Mn is not multiplied by a
                    // so there is no need to divide a
                    //
                    // However, once Maths.eq3_21 goes active then Ma will need to be
                    // divided out
                    Mn  = phi * M_factor0;
                    Mn -= StrictMath.sin(2.0 * phi) * M_factor2;
                    Mn += StrictMath.sin(4.0 * phi) * M_factor4;
                    Mn -= StrictMath.sin(6.0 * phi) * M_factor6;

                    // eq 18-17
                    Mpn  = 1.0 - esq * 0.25 - esq * esq * (3.0 / 64.0) - esq * esq * esq * (5.0 / 256.0);
                    Mpn -= 2.0 * StrictMath.cos(2.0 * phi) * ((3.0 / 8.0) * esq + (3.0 / 32.0) * esq * esq + (45.0 / 1024.0 * esq * esq * esq));
                    Mpn += 4.0 * StrictMath.cos(4.0 * phi) * ((15.0 / 256.0) * esq * esq + (45.0 / 1024.0) * esq * esq * esq);
                    Mpn -= 6.0 * StrictMath.cos(6.0 * phi) * ((35.0 / 3072.0) * esq * esq * esq);

                    phi = phi -
                            (A * (C * Mn + 1.0) - Mn - 0.5 * (Mn * Mn + B) * C) /
                            (esq * sin2phi * (Mn * Mn + B - 2.0 * A * Mn) / 4.0 *
                                    C + (A - Mn) * (C * Mpn - 2.0 / sin2phi) - Mpn);

                } while(++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phi - phi_old));

                lat[i] = phi;
                lon[i] = normalizeLonRad(StrictMath.asin(x[i] * C / a) / StrictMath.sin(phi) + lon0);
            }
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat0 = datum.getProperty("lat0");
        boolean ret = datum.getProperty("returnAllPoints") == 1.0;

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

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

        double M, cotlat, E, N, sinlat, lon_M_lon0;
        boolean forwardOk;

        for(int i = 0; i < lon.length; ++i) {
            lon_M_lon0 = normalizeLonRad(lon[i] - lon0);

            forwardOk = StrictMath.abs(lon_M_lon0) < CUTOFF_LAT_DIST | ret;

            // To actually approach 90 degrees difference lon[i] - lon0 requires an almost exponential amount of
            // time to inverse project.  Going +/- 80 degrees allows the software to back project properly
            // within the 50 iterations defined for SERIES_EXPANSION_LIMIT
            if(forwardOk && StrictMath.abs(lat[i]) < NEAR_ZERO_RAD) {
                x[i] = a * lon_M_lon0;
                y[i] = -M0;
            }
            else if(forwardOk) {
                M  = lat[i] * M_factor0;
                M -= StrictMath.sin(2.0 * lat[i]) * M_factor2;
                M += StrictMath.sin(4.0 * lat[i]) * M_factor4;
                M -= StrictMath.sin(6.0 * lat[i]) * M_factor6;
                M *= a;

                sinlat = StrictMath.sin(lat[i]);
                cotlat = cot(lat[i]);

                N = a / StrictMath.sqrt(1.0 - esq * sinlat * sinlat);
                E = lon_M_lon0 * sinlat;

                x[i] = N * cotlat * StrictMath.sin(E);
                y[i] = M - M0 + N * cotlat * (1.0 - StrictMath.cos(E));
            }
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0",
                "lat0",
                "returnAllPoints"}));
    }

}