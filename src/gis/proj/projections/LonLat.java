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

import gis.proj.Datum;
import gis.proj.Ellipsoid;
import gis.proj.Projection;

import java.util.Collections;
import java.util.Set;


/**
 * Longitude and Latitude pass through.
 *
 */

public final class LonLat implements Projection {

    public String getName() {
        return "Longitude & Latitude";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a = ellip.getProperty("a");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        for(int i = 0; i < lon.length; ++i) {
            lon[i] = x[i] / a;
            lat[i] = y[i] / a;
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a = ellip.getProperty("a");

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        for(int i = 0; i < lon.length; ++i) {
            x[i] = lon[i] * a;
            y[i] = lat[i] * a;
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return Collections.emptySet();
    }

}