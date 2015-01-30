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