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

public final class SnyderMath {

    private SnyderMath() {}

    public static enum ASPECT { POLAR, EQUATORIAL, OBLIQUE }

    // ~0.1 millimeters in degrees
    public static final double NEAR_ZERO_DEG = 1.0E-9;

    // ~0.1 millimeters in radians
    public static final double NEAR_ZERO_RAD = StrictMath.toRadians(NEAR_ZERO_DEG);

    // 50 expansions ought to be enough for most projection series expansions
    // as defined in Snyder's text
    public static final int SERIES_EXPANSION_LIMIT = 50;

    // Generalized numerical constants
    public static final double _1_DIV_2 = 1.0 / 2.0;
    public static final double _1_DIV_4 = 1.0 / 4.0;

    // Generalized PI constants
    public static final double PI         =  StrictMath.PI;
    public static final double N_PI       = -PI;
    public static final double PI_DIV_2   =  PI * _1_DIV_2;
    public static final double N_PI_DIV_2 = -PI_DIV_2;
    public static final double PI_DIV_4   =  PI * _1_DIV_4;
    public static final double PI_MUL_2   =  2.0 * PI;
    public static final double _1_DIV_PI  =  1.0 / PI;

    // Generalized E (base of natural logarithm)
    public static final double _e_ = StrictMath.E;

    // Pole Constants (RADIANS)
    public static final double POLE_RAD   =  PI_DIV_2;
    public static final double N_POLE_RAD = -POLE_RAD;

    // Generalized SQRT constants
    public static final double SQRT_2        = StrictMath.sqrt(2);
    public static final double SQRT_8        = 2.0 * SQRT_2;
    public static final double SQRT_8_DIV_PI = SQRT_8 * _1_DIV_PI;

    public static final double DEG_TO_RAD = PI / 180.0;
    public static final double RAD_TO_DEG = 180.0 / PI;

    public static ASPECT getAspectRad(double a) {
        ASPECT aspect = ASPECT.OBLIQUE;
        double pos_a = StrictMath.abs(a);

             if(POLE_RAD - pos_a < NEAR_ZERO_RAD) aspect = ASPECT.POLAR;
        else if(           pos_a < NEAR_ZERO_RAD) aspect = ASPECT.EQUATORIAL;

        return aspect;
    }

    // In proj4 this is known as mlfn
    // Closed form equation
//    public static double eq_3_21(Ellipsoid ellip, double lat) {
//        double a    = ellip.getProperty("a");
//        double esq  = ellip.getProperty("e^2");
//
//        return
//              a * (
//                      lat * (1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125)
//                      -StrictMath.sin(2.0 * lat) * (esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125)
//                      +StrictMath.sin(4.0 * lat) * ((esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125)
//                      -StrictMath.sin(6.0 * lat) * ((esq * esq * esq) * 0.011393229167)
//              );
//    }
//
//    public static double[] eq_3_21(Ellipsoid ellip, double[] lat) {
//      if(lat == null || lat.length == 0)
//          return new double[]{};
//
//      double a   = ellip.getProperty("a");
//      double esq = ellip.getProperty("esq");
//
//      double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
//      double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
//      double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
//      double M_factor6 = (esq * esq * esq) * 0.011393229167;
//
//      double M[] = new double[lat.length];
//      for(int i = 0; i < lat.length; i++) {
//          M[i]  = lat[i] * M_factor0;
//          M[i] -= StrictMath.sin(2.0 * lat[i]) * M_factor2;
//          M[i] += StrictMath.sin(4.0 * lat[i]) * M_factor4;
//          M[i] -= StrictMath.sin(6.0 * lat[i]) * M_factor6;
//          M[i] *= a;
//      }
//
//      return M;
//    }

    public static double normalizeLonDeg(double n) {
        return ((n + 180.0) % 360.0) + (n < -180.0 ? 180.0 : -180.0);
    }

    public static double normalizeLonRad(double n) {
        return ((n + PI) % PI_MUL_2) + (n < N_PI ? PI : N_PI);
    }

//    public static double exp(double a) {
//      return StrictMath.exp(a);
//    }
//
//    public static double copySign(double magnitude, double sign) {
//      return StrictMath.copySign(magnitude, sign);
//    }

    public static double cot(double a) {
        return 1.0 / StrictMath.tan(a);
    }

//    public static double abs(double a) {
//        return StrictMath.abs(a);
//    }
//
//    public static double acos(double a) {
//        return StrictMath.acos(a);
//    }
//
//    public static double asin(double a) {
//        return StrictMath.asin(a);
//    }
//
//    public static double atan(double a) {
//        return StrictMath.atan(a);
//    }
//
//    public static double atan2(double y, double x) {
//        return StrictMath.atan2(y, x);
//    }
//
//    public static double cos(double a) {
//        return StrictMath.cos(a);
//    }
//
//    public static double sin(double a) {
//        return StrictMath.sin(a);
//    }
//
//    public static double sinh(double a) {
//        return StrictMath.sinh(a);
//    }
//
//    public static double tan(double a) {
//        return StrictMath.tan(a);
//    }
//
//    public static double ln(double a) {
//        return StrictMath.log(a);
//    }
//
//    public static double sqrt(double a) {
//        return StrictMath.sqrt(a);
//    }
//
//    public static double hypot(double x, double y) {
//        // Using StrictMath.hypot is very slow compared to below
//        return StrictMath.sqrt(x * x + y * y);
//    }
//
//    public static double pow(double a, double b) {
//        return StrictMath.pow(a, b);
//    }
//
//    public static double normalizeLatDeg(double n) {
//        if(n < -90.0)
//            return -90.0;
//        if(n > 90.0)
//            return 90.0;
//        return n;
//    }
//
//    public static double normalizeLatRad(double n) {
//        if(n < NEG_HLF_PI)
//            return NEG_HLF_PI;
//        if(n > POS_HLF_PI)
//            return POS_HLF_PI;
//        return n;
//    }

}