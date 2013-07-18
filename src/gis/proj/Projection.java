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
package gis.proj;

import java.util.Set;

public interface Projection {

    /**
     * 
     * @param x[] (meters)
     * @param y[] (meters)
     * @param e Ellipsoid
     * @param d Datum
     * @return
     *      [0][N] of lon or lambda values {lon0, lon1, ..., lonN} (radians)
     *      [1][N] of lat or phi    values {lat0, lat1, ..., latN} (radians)
     */
    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum);

    /**
     * 
     * @param lon[] or lambda (radians)
     * @param lat[] or phi    (radians)
     * @param e Ellipsoid
     * @param d Datum
     * @return
     *      [0][N] of X values {x0, x1, ..., xN} (meters)
     *      [1][N] of Y values {y0, y1, ..., yN} (meters)
     */
    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum);

    public String getName();

    public Set<String> getDatumProperties();

}