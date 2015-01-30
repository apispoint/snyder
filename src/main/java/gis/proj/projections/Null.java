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