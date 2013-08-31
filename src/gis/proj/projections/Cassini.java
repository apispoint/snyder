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
 * Pages 92-95
 *
 */

public final class Cassini implements Cylindrical {

    public String getName() {
        return "Cassini";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat0 = datum.getProperty("lat0");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
        double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
        double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
        double M_factor6 = (esq * esq * esq) * 0.011393229167;

        double M0  = lat0 * M_factor0;
        M0 -= StrictMath.sin(2.0 * lat0) * M_factor2;
        M0 += StrictMath.sin(4.0 * lat0) * M_factor4;
        M0 -= StrictMath.sin(6.0 * lat0) * M_factor6;
        M0 *= a;

        double E1     = (1.0 - StrictMath.sqrt(1.0 - esq)) / (1.0 + StrictMath.sqrt(1.0 - esq));
        double MU_DIV = (a * (1.0 - esq * 0.25 - esq * esq * (3.0 / 64.0) - esq * esq * esq * (5.0 / 256.0)));

        double M1, mu1, lat1, sinlat1, T1, N1, R1, D;

        for(int i = 0; i < lon.length; ++i) {
            M1  = M0 + y[i];
            mu1 = M1 / MU_DIV;

            lat1  = mu1 + ((3.0/2.0) * E1 - (27.0 / 32.0)*(E1 * E1 * E1)) * StrictMath.sin(2.0 * mu1);
            lat1 += (((21.0/16.0) * E1 * E1) - ((55.0/32.0) * E1 * E1 * E1 * E1)) * StrictMath.sin(4.0 * mu1);
            lat1 += ((151.0/96.0) * E1 * E1 * E1) * StrictMath.sin(6.0 * mu1);
            lat1 += ((1097.0/512.0) * E1 * E1 * E1 * E1) * StrictMath.sin(8.0 * mu1);

            sinlat1 = StrictMath.sin(lat1);

            T1 = StrictMath.tan(lat1) * StrictMath.tan(lat1);
            N1 = a / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);
            R1 = a * (1.0 - esq) / StrictMath.pow(1.0 - esq * sinlat1 * sinlat1, 3.0 /2.0);
            D = x[i] / N1;

            if(
                    (PI_DIV_2 - StrictMath.abs(lat1) < NEAR_ZERO_RAD && !(lat1 < 0.0 ^ y[i] < 0.0)) ||
                    (PI_DIV_2 - StrictMath.abs(lat[i]) < NEAR_ZERO_RAD && !(lat[i] < 0.0 ^ y[i] < 0.0))
                    )
                lon[i] = lon0;
            else
                lon[i] = lon0 + (D - T1 * (D*D*D)/3.0 + (1.0 + 3.0 * T1) * T1 * (D*D*D*D*D) / 15.0) / StrictMath.cos(lat1);

            lat[i] = lat1 - (N1 * StrictMath.tan(lat1) / R1) * (D*D *0.5 - (1.0 + 3.0 * T1) * (D*D*D*D) / 24.0);
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat0 = datum.getProperty("lat0");

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
        double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
        double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
        double M_factor6 = (esq * esq * esq) * 0.011393229167;

        double M0  = lat0 * M_factor0;
        M0 -= StrictMath.sin(2.0 * lat0) * M_factor2;
        M0 += StrictMath.sin(4.0 * lat0) * M_factor4;
        M0 -= StrictMath.sin(6.0 * lat0) * M_factor6;
        M0 *= a;

        double one_M_esq = 1.0 - esq;

        double coslat, sinlat, tanlat, M, N, T, A, C;

        for(int i = 0; i < lon.length; ++i) {
            M  = lat[i] * M_factor0;
            M -= StrictMath.sin(2.0 * lat[i]) * M_factor2;
            M += StrictMath.sin(4.0 * lat[i]) * M_factor4;
            M -= StrictMath.sin(6.0 * lat[i]) * M_factor6;
            M *= a;

            coslat = StrictMath.cos(lat[i]);
            sinlat = StrictMath.sin(lat[i]);
            tanlat = StrictMath.tan(lat[i]);

            N = a / StrictMath.sqrt(1.0 - esq * sinlat * sinlat);
            T = tanlat * tanlat;
            A = normalizeLonRad(lon[i] - lon0) * coslat;
            C = esq * (coslat * coslat) / one_M_esq;

            x[i] = N * (A - T * (A * A * A) / 6.0 - (8.0 - T + 8.0 * C) * T * (A * A * A * A * A) / 120.0);
            y[i] = M - M0 + N * tanlat * (A * A * 0.5 + (5.0 - T + 6.0 * C) * (A * A * A * A) / 24.0);
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0",
                "lat0"}));
    }

}