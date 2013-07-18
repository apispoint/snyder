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
 * Null projection (e.g. pass through).
 *
 */

public final class Null implements Projection {

    public String getName() {
        return "Null";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        System.arraycopy(x, 0, lon, 0, x.length);
        System.arraycopy(y, 0, lat, 0, y.length);

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        System.arraycopy(lon, 0, x, 0, lon.length);
        System.arraycopy(lat, 0, y, 0, lat.length);

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return Collections.emptySet();
    }

}