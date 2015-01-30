/*
Snyder Projection Implementation

Copyright (c) 2012-2015, APIS Point, LLC

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
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