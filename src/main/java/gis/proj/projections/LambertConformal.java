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

import static gis.proj.SnyderMath.*;

import gis.proj.Conic;
import gis.proj.Datum;
import gis.proj.Ellipsoid;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;


/**
 * All formulas adopted from USGS Professional Paper 1395.
 *
 * References:
 *
 * Map Projections - A Working Manual, USGS Professional Paper 1395
 * John P. Snyder
 * Pages 104-110
 *
 */

public final class LambertConformal implements Conic {

    public String getName() {
        return "Lambert Conformal";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double e   = ellip.getProperty("e");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat0 = datum.getProperty("lat0");
        double lat1 = datum.getProperty("lat1");
        double lat2 = datum.getProperty("lat2");

        //        boolean useGeocentricLat = datum.getProperty("useGeocentricLat") == 1.0;

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double hlf_e = e * 0.5;

        double rho, theta, phi, phi_old, esinphi;
        int j;

        double sinlat0 = StrictMath.sin(lat0),
               sinlat1 = StrictMath.sin(lat1),
               sinlat2 = StrictMath.sin(lat2);

        double m1 = StrictMath.cos(lat1) / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);
        double m2 = StrictMath.cos(lat2) / StrictMath.sqrt(1.0 - esq * sinlat2 * sinlat2);

        double t, esinlat0, esinlat1, esinlat2;
        double t0 = StrictMath.tan(PI_DIV_4 - lat0 * 0.5) / StrictMath.pow((1.0 - (esinlat0 = e * sinlat0)) / (1.0 + esinlat0), hlf_e);
        double t1 = StrictMath.tan(PI_DIV_4 - lat1 * 0.5) / StrictMath.pow((1.0 - (esinlat1 = e * sinlat1)) / (1.0 + esinlat1), hlf_e);
        double t2 = StrictMath.tan(PI_DIV_4 - lat2 * 0.5) / StrictMath.pow((1.0 - (esinlat2 = e * sinlat2)) / (1.0 + esinlat2), hlf_e);

        double nu = sinlat1;

        if(NEAR_ZERO_RAD < StrictMath.abs(lat1 - lat2))
            nu = (StrictMath.log(m1) - StrictMath.log(m2)) / (StrictMath.log(t1) - StrictMath.log(t2));

        double inv_nu = 1.0 / nu;
        double nuFactor = nu < 0.0 ? -1.0 : 1.0;

        double F = m1 / (nu * StrictMath.pow(t1, nu));
        double aF = a * F;
        double rho0 = aF * StrictMath.pow(t0, nu);

        for(int i = 0; i < lon.length; ++i) {
            rho = StrictMath.copySign(StrictMath.hypot(x[i], rho0 - y[i]), nu);

            t = StrictMath.pow(rho / aF, inv_nu);
            theta = StrictMath.atan2(x[i] * nuFactor, rho0 * nuFactor - y[i] * nuFactor);

            lon[i] = normalizeLonRad(theta / nu + lon0);

            phi = PI_DIV_2 - 2.0 * StrictMath.atan(t);
            j = 0;
            do {
                esinphi = e * StrictMath.sin(phi);
                phi_old = phi;
                phi = PI_DIV_2 - 2.0 * StrictMath.atan(t * StrictMath.pow((1.0 - esinphi) / (1.0 + esinphi), hlf_e));
            } while (++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phi - phi_old));

            lat[i] = phi;
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double e   = ellip.getProperty("e");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat0 = datum.getProperty("lat0");
        double lat1 = datum.getProperty("lat1");
        double lat2 = datum.getProperty("lat2");

        //        boolean useGeocentricLat = datum.getProperty("useGeocentricLat") == 1.0;

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double hlf_e = e * 0.5;

        double rho, theta, esinlat;

        //TODO undefined when it's a cylinder
        //        if(StrictMath.abs(POS_POLE_RADIANS - lat1) < NEAR_ZERO_RAD)
        //            lat1 -= NEAR_ZERO_RAD;

        double sinlat,
        sinlat0 = StrictMath.sin(lat0),
        sinlat1 = StrictMath.sin(lat1),
        sinlat2 = StrictMath.sin(lat2);

        double m1 = StrictMath.cos(lat1) / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);
        double m2 = StrictMath.cos(lat2) / StrictMath.sqrt(1.0 - esq * sinlat2 * sinlat2);

        double t, esinlat0, esinlat1, esinlat2;
        double t0 = StrictMath.tan(PI_DIV_4 - lat0 * 0.5) / StrictMath.pow((1.0 - (esinlat0 = e * sinlat0)) / (1.0 + esinlat0), hlf_e);
        double t1 = StrictMath.tan(PI_DIV_4 - lat1 * 0.5) / StrictMath.pow((1.0 - (esinlat1 = e * sinlat1)) / (1.0 + esinlat1), hlf_e);
        double t2 = StrictMath.tan(PI_DIV_4 - lat2 * 0.5) / StrictMath.pow((1.0 - (esinlat2 = e * sinlat2)) / (1.0 + esinlat2), hlf_e);

        double nu = sinlat1;

        if(NEAR_ZERO_RAD < StrictMath.abs(lat1 - lat2))
            nu = (StrictMath.log(m1) - StrictMath.log(m2)) / (StrictMath.log(t1) - StrictMath.log(t2));

        double F = m1 / (nu * StrictMath.pow(t1, nu));
        double aF = a * F;
        double rho0 = aF * StrictMath.pow(t0, nu);

        for(int i = 0; i < lon.length; ++i) {
            sinlat = StrictMath.sin(lat[i]);
            esinlat = e * sinlat;
            t = StrictMath.tan(PI_DIV_4 - lat[i] * 0.5) / StrictMath.pow((1.0 - esinlat) / (1.0 + esinlat), hlf_e);

            theta = nu * normalizeLonRad(lon[i] - lon0);
            rho = (NEAR_ZERO_RAD <= POLE_RAD - StrictMath.abs(lat[i])) ? aF * StrictMath.pow(t, nu) : 0;

            x[i] = rho * StrictMath.sin(theta);
            y[i] = rho0 - rho * StrictMath.cos(theta);
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0",
                "lat0",
                "lat1",
                "lat2"}));
    }

}