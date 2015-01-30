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

public class ComplexNumber implements Cloneable {

    private double re;
    private double im;

    public ComplexNumber() {
        this(0.0, 0.0);
    }

    public ComplexNumber(double re, double im) {
        this.re = re;
        this.im = im;
    }

    public ComplexNumber add(ComplexNumber b) {
        return new ComplexNumber(re + b.re, im + b.im);
    }

    public ComplexNumber sub(ComplexNumber b) {
        return new ComplexNumber(re - b.re, im - b.im);
    }

    public ComplexNumber mul(ComplexNumber b) {
        return new ComplexNumber(
                re * b.re - im * b.im,
                re * b.im + im * b.re);
    }

    public ComplexNumber mul(double d) {
        return new ComplexNumber(d * re, d * im);
    }

    public ComplexNumber div(ComplexNumber b) {
        double dem = b.re * b.re + b.im * b.im;
        return new ComplexNumber(
                (re * b.re + im * b.im) / dem,
                (im * b.re - re * b.im) / dem);
    }

    public ComplexNumber div(double d) {
        return new ComplexNumber(re / d, im / d);
    }

    public ComplexNumber conjugate() {
        return new ComplexNumber(re, -im);
    }

    public double abs() {
        return StrictMath.sqrt(re * re + im * im);
    }

    public double absSq() {
        return re * re + im * im;
    }

    public double phase() {
        return StrictMath.atan2(im, re);
    }

    public ComplexNumber clone() {
        return new ComplexNumber(re, im);
    }

    public double re() {
        return re;
    }

    public double im() {
        return im;
    }

    public void setReal(double r) {
        re = r;
    }

    public void setImaginary(double i) {
        im = i;
    }

}