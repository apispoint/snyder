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
 * Pages 98-103
 *
 */

public final class AlbersEA implements Conic {

    public String getName() {
        return "Albers Equal-Area";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double asq = ellip.getProperty("a^2");
        double e   = ellip.getProperty("e");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat0 = datum.getProperty("lat0");
        double lat1 = datum.getProperty("lat1");
        double lat2 = datum.getProperty("lat2");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double twoe = 2.0 * e;
        double rho, theta, phi, phi_old, sinphi;
        int j;

        double sinlat0 = StrictMath.sin(lat0),
               sinlat1 = StrictMath.sin(lat1),
               sinlat2 = StrictMath.sin(lat2);

        double one_M_esq = 1.0 - esq;
        double one_over_2e = 1.0 / twoe;

        double q, q0, q1, q2, qTest, factor;

        boolean eNearZeroOrSphere = ellip.isSphere();

        if(eNearZeroOrSphere == false) {
            q0 = qFn(one_M_esq, one_over_2e, esq, sinlat0, sinlat0 * e);
            q1 = qFn(one_M_esq, one_over_2e, esq, sinlat1, sinlat1 * e);
            q2 = qFn(one_M_esq, one_over_2e, esq, sinlat2, sinlat2 * e);
        } else {
            q0 = 2.0 * sinlat0;
            q1 = 2.0 * sinlat1;
            q2 = 2.0 * sinlat2;
        }

        double m1 = StrictMath.cos(lat1) / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);
        double m2 = StrictMath.cos(lat2) / StrictMath.sqrt(1.0 - esq * sinlat2 * sinlat2);

        double nu = sinlat1;

        if(NEAR_ZERO_RAD < StrictMath.abs(lat1 - lat2))
            nu = (m1 * m1 - m2 * m2) / (q2 - q1);

        double twonu = 2.0 * nu;

        double C = m1 * m1 + nu * q1;
        double rho0 = a * StrictMath.sqrt(C - nu * q0) / nu;

        double nuFactor = nu < 0.0 ? -1.0 : 1.0;

        for(int i = 0; i < lon.length; ++i) {
            rho = StrictMath.hypot(x[i], rho0 - y[i]);

            q = (C - rho * rho * nu * nu / asq) / nu;

            theta = StrictMath.atan2(x[i] * nuFactor, rho0 * nuFactor - y[i] * nuFactor);

            lon[i] = normalizeLonRad(lon0 + theta / nu);

            if(eNearZeroOrSphere == false) {
                qTest = (1.0 - ((1.0 - esq) / twoe) * StrictMath.log((1.0 - e) / (1.0 + e)));

                if(StrictMath.abs(q - qTest) < NEAR_ZERO_DEG)
                {
                    phi = PI_DIV_2 * q < 0 ? -1 : 1;
                } else {
                	phi = StrictMath.sin(q * 0.5);
                    j = 0;
                    do {
                        sinphi = StrictMath.sin(phi);
                        // eq3_16
                        factor = (1.0 - esq * sinphi * sinphi);

                        phi_old = phi;
                        phi = phi + ((factor * factor) / (2.0 * StrictMath.cos(phi))) * (
                                q / one_M_esq - sinphi / factor + one_over_2e * StrictMath.log((1.0 - e * sinphi) / (1.0 + e * sinphi))
                                );
                    } while (++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phi - phi_old));
                }
            }
            else 
                phi = StrictMath.asin((C - StrictMath.pow(rho * nu / a, 2)) / twonu);

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

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double rho, theta, esinlat;

        double sinlat,
        sinlat0 = StrictMath.sin(lat0),
        sinlat1 = StrictMath.sin(lat1),
        sinlat2 = StrictMath.sin(lat2);

        double one_M_esq = 1.0 - esq;
        double one_over_2e = 1.0 / (2.0 * e);

        double q, q0, q1, q2;

        boolean eNearZeroOrSphere = ellip.isSphere();

        if(eNearZeroOrSphere == false) {
            q0 = qFn(one_M_esq, one_over_2e, esq, sinlat0, sinlat0 * e);
            q1 = qFn(one_M_esq, one_over_2e, esq, sinlat1, sinlat1 * e);
            q2 = qFn(one_M_esq, one_over_2e, esq, sinlat2, sinlat2 * e);
        } else {
            q0 = 2.0 * sinlat0;
            q1 = 2.0 * sinlat1;
            q2 = 2.0 * sinlat2;
        }

        double m1 = StrictMath.cos(lat1) / StrictMath.sqrt(1.0 - esq * sinlat1 * sinlat1);
        double m2 = StrictMath.cos(lat2) / StrictMath.sqrt(1.0 - esq * sinlat2 * sinlat2);

        double nu = sinlat1;

        if(NEAR_ZERO_RAD < StrictMath.abs(lat1 - lat2))
            nu = (m1 * m1 - m2 * m2) / (q2 - q1);

        double C = m1 * m1 + nu * q1;
        double rho0 = a * StrictMath.sqrt(C - nu * q0) / nu;

        for(int i = 0; i < lon.length; ++i) {
            sinlat = StrictMath.sin(lat[i]);
            esinlat = e * sinlat;

            q = eNearZeroOrSphere == false ? one_M_esq * (
            		sinlat / (1.0 - esq * sinlat * sinlat) - one_over_2e *
            		StrictMath.log((1.0 - esinlat) / (1.0 + esinlat))
            		) : 2.0 * sinlat;

            theta = nu * normalizeLonRad(lon[i] - lon0);
            rho = a * StrictMath.sqrt(C - nu * q) / nu;

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

    private double qFn(double one_M_esq, double one_over_2e, double esq, double sinlat, double esinlat) {
        // (3-12)
        return one_M_esq * (
                        sinlat / (1.0 - esq * sinlat * sinlat) - one_over_2e *
                        StrictMath.log((1.0 - esinlat) / (1.0 + esinlat))
                        );
    }

}