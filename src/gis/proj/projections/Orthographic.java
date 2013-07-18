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
 * Pages 145-153
 *
 */

public final class Orthographic implements Azimuthal, Spherical {

    public String getName() {
        return "Orthographic";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double R = ellip.getProperty("R");

        double lon0 = datum.getProperty("lon0");
        double lat1 = datum.getProperty("lat1");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double coslat1 = StrictMath.cos(lat1);
        double sinlat1 = StrictMath.sin(lat1);

        ASPECT aspect = getAspectRad(lat1);
        double Rp = lat1 < 0.0 ? 1.0 : -1.0;

        double rho, c, sinc, cosc;

        for(int i = 0; i < x.length; ++i) {
            rho = StrictMath.hypot(x[i], y[i]);

            if(rho < NEAR_ZERO_DEG == false) {
                c = StrictMath.asin(rho / R);

                sinc = StrictMath.sin(c);
                cosc = StrictMath.cos(c);

                switch(aspect) {
                case POLAR:
                    lon[i] = lon0 + StrictMath.atan2(x[i], y[i] * Rp);
                    break;
                default:
                    lon[i] = lon0 + StrictMath.atan2(x[i] * sinc, rho * coslat1 * cosc - y[i] * sinlat1 * sinc);
                }

                lat[i] = StrictMath.asin(cosc * sinlat1 + y[i] * sinc * coslat1 / rho);
            }
            else {
                lon[i] = lon0;
                lat[i] = lat1;
            }

            lon[i] = normalizeLonRad(lon[i]);
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double R = ellip.getProperty("R");

        double lon0 = datum.getProperty("lon0");
        double lat1 = datum.getProperty("lat1");
        boolean ret = datum.getProperty("returnAllPoints") == 1.0;

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double sinlat1 = StrictMath.sin(lat1);
        double coslat1 = StrictMath.cos(lat1);

        double c, sinlat, coslat;
        double lon_M_lon0, coslon_M_lon0;

        ASPECT aspect = getAspectRad(lat1);
        double Rp = lat1 < 0.0 ? R : -R;

        for(int i = 0; i < lon.length; ++i) {
            coslat = StrictMath.cos(lat[i]);
            sinlat = StrictMath.sin(lat[i]);

            lon_M_lon0 = normalizeLonRad(lon[i] - lon0);
            coslon_M_lon0 = StrictMath.cos(lon_M_lon0);

            // (5-3)
            c = sinlat1 * sinlat + coslat1 * coslat * coslon_M_lon0;

            if((ret == false && c >= 0.0) || ret == true) {
                x[i] = R * coslat * StrictMath.sin(lon_M_lon0);

                switch(aspect) {
                case OBLIQUE:
                    y[i] = R * (coslat1 * sinlat - sinlat1 * coslat * coslon_M_lon0);
                    break;
                case EQUATORIAL:
                    y[i] = R * sinlat;
                    break;
                default:
                    y[i] = Rp * coslat * coslon_M_lon0;
                }
            }
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0",
                "lat1",
                "returnAllPoints"}));
    }

}