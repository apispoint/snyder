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