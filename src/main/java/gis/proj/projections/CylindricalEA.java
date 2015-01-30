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

import gis.proj.Cylindrical;
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
 * Pages 76-85
 *
 */

public final class CylindricalEA implements Cylindrical {

    public String getName() {
        return "Cylindrical Equal-Area";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double e   = ellip.getProperty("e");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lats = datum.getProperty("lats");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        // Handle the spherical case when 1/(2e) or (1/0 -> Inf) given
        // eccentricity (e) = 0 (sphere or circle)
        if(ellip.isSphere()) {
            double coslat1 = StrictMath.cos(lats);
            double acoslat1 = a * coslat1;

            // Not using R because this is still considered an ellipsoidal
            // solution regardless of eccentricity equaling 0
            for(int i = 0; i < lon.length; ++i) {
                lon[i] = normalizeLonRad(x[i] / acoslat1 + lon0);
                lat[i] = StrictMath.asin((y[i] / a) * coslat1);
            }
        } else {
            double sinlat1 = StrictMath.sin(lats);

            double k0 = StrictMath.cos(lats) / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);
            double TWO_k0 = 2.0 * k0;
            double ak0 = a * k0;

            double one_M_esq = 1.0 - esq;
            double one_over_2e = 1.0 / (2.0 * e);

            double qp = one_M_esq * (
                    StrictMath.sin(PI_DIV_2) / (1.0 - esq * StrictMath.sin(PI_DIV_2) * StrictMath.sin(PI_DIV_2)) - one_over_2e * StrictMath.log((1.0 - e * StrictMath.sin(PI_DIV_2)) / (1.0 + e * StrictMath.sin(PI_DIV_2)))
                    );
            double aqp = a * qp;

            double q, beta, phi, phi_old, factor, sinphi;
            int j;

            for(int i = 0; i < lon.length; ++i) {
                beta = StrictMath.asin(TWO_k0 * y[i] / aqp);

                q = qp * StrictMath.sin(beta);
                phi = q * 0.5;
                j = 0;
                do {
                    phi_old = phi;
                    sinphi = StrictMath.sin(phi);
                    // (3-16)
                    factor = (1.0 - esq * sinphi * sinphi);
                    phi = phi + ((factor * factor) / (2.0 * StrictMath.cos(phi))) * (
                            q / one_M_esq - sinphi / factor + one_over_2e * StrictMath.log((1.0 - e * sinphi) / (1.0 + e * sinphi))
                            );
                } while (++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phi - phi_old));

                lon[i] = normalizeLonRad(lon0 + x[i] / ak0);
                lat[i] = phi;
            }
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double e   = ellip.getProperty("e");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lats = datum.getProperty("lats");

        double x[] = new double[lon.length];
        double y[] = new double[lat.length];

        // Handle the spherical case when 1/(2e) or (1/0 -> Inf) given
        // eccentricity (e) = 0 (sphere or circle)
        if(ellip.isSphere()) {
            double coslat1 = StrictMath.cos(lats);

            // Not using R because this is still considered an ellipsoidal
            // solution regardless of eccentricity equaling 0
            for(int i = 0; i < lon.length; ++i) {
                x[i] = a * normalizeLonRad(lon[i] - lon0) * coslat1;
                y[i] = a * StrictMath.sin(lat[i]) / coslat1;
            }
        } else {
            double sinlat1 = StrictMath.sin(lats);

            double k0 = StrictMath.cos(lats) / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);
            double TWO_k0 = 2.0 * k0;
            double ak0 = a * k0;

            double one_M_esq = 1.0 - esq;
            double one_over_2e = 1.0 / (2.0 * e);

            double sinlat, q, esinlat;

            for(int i = 0; i < lon.length; ++i) {
                sinlat = StrictMath.sin(lat[i]);
                esinlat = e * sinlat;

                //eq3_12
                q = one_M_esq * (
                        sinlat / (1.0 - esq * sinlat * sinlat) - one_over_2e * StrictMath.log((1.0 - esinlat) / (1.0 + esinlat))
                        );

                x[i] = ak0 * normalizeLonRad(lon[i] - lon0);
                y[i] = a * q / TWO_k0;
            }
        }
        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0",
                "lats"}));
    }

}