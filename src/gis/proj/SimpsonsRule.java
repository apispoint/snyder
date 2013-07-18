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

public final class SimpsonsRule {

    public static interface Function {
        public double[] f(double[] parameters, double x);
    }

    private static SimpsonsRule sr = new SimpsonsRule();

    private SimpsonsRule() {}

    public static SimpsonsRule getInstance() {
        return sr;
    }

    public double[] approximate(Function function, double[] parameters, double a, double b, int n)
            throws IllegalArgumentException {
        if(n < 2 || (n & 0x01) == 1)
            throw new IllegalArgumentException("'n' must be greater than zero and an even integer.");

        double[] value = vectorAdd(function.f(parameters, a), function.f(parameters, b));

        double dX = (b - a) / (double) n;
        double iX = a + dX;

        for(int i = 1, scalar = 2; i < n; i++, iX += dX) {
            value = vectorAdd(value, vectorMulScalar(function.f(parameters, iX), (scalar = scalar % 4 + 2)));
        }

        return vectorMulScalar(value, dX / 3.0);
    }

    /*
     * Vector operations
     */
    private double[] vectorAdd(double[] v1, double[] v2) {
        if(v1 == null || v2 == null || v1.length != v2.length)
            throw new IllegalArgumentException("Vectors must not be null and equal length.");

        double[] v = new double[v1.length];

        for(int i = 0; i < v.length; i++)
            v[i] = v1[i] + v2[i];

        return v;
    }

    private double[] vectorMulScalar(double[] v1, double scalar) {
        if(v1 == null)
            throw new IllegalArgumentException("Vector must not be null.");

        double[] v = new double[v1.length];

        for(int i = 0; i < v.length; i++)
            v[i] = v1[i] * scalar;

        return v;
    }

}