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

import gis.proj.Azimuthal;
import gis.proj.Datum;
import gis.proj.Ellipsoid;
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
 * Pages 164-168
 *
 */

public final class Gnomonic implements Azimuthal, Spherical {

    public String getName() {
        return "Gnomonic";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double R    = ellip.getProperty("R");

        double lon0 = datum.getProperty("lon0");
        double lat1 = datum.getProperty("lat1");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double coslat1 = StrictMath.cos(lat1);
        double sinlat1 = StrictMath.sin(lat1);

        double rho, c;

        // Polar aspect
        if (POLE_RAD - StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
            double factorY = lat1 > 0.0 ? -1.0 : 1.0;

            for(int i = 0; i < x.length; ++i) {
                rho = StrictMath.hypot(x[i], y[i]);
                c = StrictMath.atan(rho / R);

                lon[i] = normalizeLonRad(lon0 + StrictMath.atan2(x[i], factorY * y[i]));
                lat[i] = StrictMath.asin(StrictMath.cos(c) * sinlat1 + (y[i] * StrictMath.sin(c) * coslat1 / rho));
            }
        }
        // Oblique et Equatorial aspects
        else {
            double cosc, sinc;

            for(int i = 0; i < x.length; ++i) {
                rho = StrictMath.hypot(x[i], y[i]);
                c = StrictMath.atan(rho / R);

                cosc = StrictMath.cos(c);
                sinc = StrictMath.sin(c);

                lon[i] = normalizeLonRad(lon0 + StrictMath.atan2(x[i] * sinc, (rho * coslat1 * cosc - y[i] * sinlat1 * sinc)));
                lat[i] = StrictMath.asin(cosc * sinlat1 + (y[i] * sinc * coslat1 / rho));
            }
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double R    = ellip.getProperty("R");

        double lon0 = datum.getProperty("lon0");
        double lat1 = datum.getProperty("lat1");
        boolean ret = datum.getProperty("returnAllPoints")                == 1.0;
        boolean mrk = datum.getProperty("returnNonVisiblePointsAsNegInf") == 1.0;

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double sinlat1 = StrictMath.sin(lat1);
        double coslat1 = StrictMath.cos(lat1);

        double lon_M_lon0, coslat, sinlat, c, coslon_M_lon0;

        // Equatorial aspect
        if(StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
            for(int i = 0; i < lon.length; ++i) {
                coslat = StrictMath.cos(lat[i]);
                sinlat = StrictMath.sin(lat[i]);

                lon_M_lon0 = normalizeLonRad(lon[i] - lon0);
                coslon_M_lon0 = StrictMath.cos(lon_M_lon0);

                c = sinlat1 * sinlat + coslat1 * coslat * coslon_M_lon0;
                if((ret == false && c >= 0.0) || ret == true) {
                    x[i] = R * StrictMath.tan(lon_M_lon0);
                    y[i] = R * StrictMath.tan(lat[i]) / coslon_M_lon0;
                }
                else if (mrk == true) {
                    x[i] = Double.NEGATIVE_INFINITY;
                    y[i] = Double.NEGATIVE_INFINITY;
                }
            }
        }
        // Polar aspect
        else if (POLE_RAD - StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
            double Rp = lat1 < 0.0 ? -R : R, cotlat;

            for(int i = 0; i < lon.length; ++i) {
                coslat = StrictMath.cos(lat[i]);
                sinlat = StrictMath.sin(lat[i]);
                cotlat = cot(lat[i]);

                lon_M_lon0 = normalizeLonRad(lon[i] - lon0);
                coslon_M_lon0 = StrictMath.cos(lon_M_lon0);

                c = sinlat1 * sinlat + coslat1 * coslat * coslon_M_lon0;
                if((ret == false && c >= 0.0) || ret == true) {
                    x[i] =  Rp * cotlat * StrictMath.sin(lon_M_lon0);
                    y[i] = -Rp * cotlat * coslon_M_lon0;
                }
                else if (mrk == true) {
                    x[i] = Double.NEGATIVE_INFINITY;
                    y[i] = Double.NEGATIVE_INFINITY;
                }
            }
        }
        // Oblique aspect
        else {
            double kp;

            for(int i = 0; i < lon.length; ++i) {
                coslat = StrictMath.cos(lat[i]);
                sinlat = StrictMath.sin(lat[i]);

                lon_M_lon0 = normalizeLonRad(lon[i] - lon0);
                coslon_M_lon0 = StrictMath.cos(lon_M_lon0);

                lon_M_lon0 = normalizeLonRad(lon[i] - lon0);
                coslat = StrictMath.cos(lat[i]);
                c = sinlat1 * sinlat + coslat1 * coslat * coslon_M_lon0; //eq5_3

                if((ret == false && c >= 0.0) || ret == true) {
                    kp = 1.0 / c;
                    x[i] = R * kp * coslat * StrictMath.sin(lon_M_lon0);
                    y[i] = R * kp * (coslat1 * sinlat - sinlat1 * coslat * coslon_M_lon0);
                }
                else if (mrk == true) {
                    x[i] = Double.NEGATIVE_INFINITY;
                    y[i] = Double.NEGATIVE_INFINITY;
                }
            }
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0",
                "lat1",
                "returnAllPoints",
                "returnNonVisiblePointsAsNegInf"}));
    }

}