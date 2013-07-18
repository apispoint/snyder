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
 * Pages 66-75
 *
 */

public final class ObliqueMercator implements Cylindrical {

    public String getName() {
        return "Oblique Mercator";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double e   = ellip.getProperty("e");
        double esq = ellip.getProperty("e^2");

        double k0   = datum.getProperty("k0");
        double lat0 = datum.getProperty("lat0");
        double x0   = datum.getProperty("x0");
        double y0   = datum.getProperty("y0");

        // Used for Central Line
        double lon1 = datum.getProperty("lon1");
        double lat1 = datum.getProperty("lat1");
        double lon2 = datum.getProperty("lon2");
        double lat2 = datum.getProperty("lat2");

        // Used for Central Point
        double lonc   = datum.getProperty("lonc");
        double alphac = datum.getProperty("alphac");

        boolean useCentralPoint = datum.getProperty("useCentralPoint") == 1.0;

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double hlf_e = e * 0.5;

        double one_M_esq = 1.0 - esq;
        double sqrt_one_M_esq = StrictMath.sqrt(one_M_esq);

        double coslat0 = StrictMath.cos(lat0);
        double sinlat0 = StrictMath.sin(lat0);

        double esinlat0                = e * sinlat0;
        double one_M_esqsinlat0sinlat0 = 1.0 - esq * sinlat0 * sinlat0;
        double coslat0_pow_4           = coslat0 * coslat0 * coslat0 * coslat0;

        // Common variable for both projection alternatives
        double B  = StrictMath.sqrt(1.0 + esq * coslat0_pow_4 / one_M_esq);
        double A  = a * B * k0 * sqrt_one_M_esq / one_M_esqsinlat0sinlat0;
        double t0 = StrictMath.tan(PI_DIV_4 - (lat0 * 0.5)) / StrictMath.pow((1.0 - esinlat0) / (1.0 + esinlat0), hlf_e);
        double D  = B * sqrt_one_M_esq / (coslat0 * StrictMath.sqrt(one_M_esqsinlat0sinlat0));
        double DD = ((DD = D * D) < 1.0 ? 1.0 : DD);

        double D_P_sqrt_DD_M_one = D + StrictMath.copySign(StrictMath.sqrt(DD - 1.0), lat0);
        double E                 = D_P_sqrt_DD_M_one * StrictMath.pow(t0, B);

        double one_over_B = 1.0 / B;

        double gamma0, lon0;

        if(useCentralPoint == false) {
            double esinlat1 = e * StrictMath.sin(lat1);
            double esinlat2 = e * StrictMath.sin(lat2);

            double t1 = StrictMath.tan(PI_DIV_4 - (lat1 * 0.5)) / StrictMath.pow((1.0 - esinlat1) / (1.0 + esinlat1), hlf_e);
            double t2 = StrictMath.tan(PI_DIV_4 - (lat2 * 0.5)) / StrictMath.pow((1.0 - esinlat2) / (1.0 + esinlat2), hlf_e);

            double H = StrictMath.pow(t1, B);
            double L = StrictMath.pow(t2, B);
            double F = E / H;
            double G = (F - 1.0 / F) * 0.5;
            double J = (E * E - L * H) / (E * E + L * H);
            double P = (L - H) / (L + H);

            lon0   = (lon1 + lon2) * 0.5 - StrictMath.atan(J * StrictMath.tan(B * normalizeLonRad(lon1 - lon2) * 0.5) / P) / B;
            gamma0 = StrictMath.atan(StrictMath.sin(B * (lon1 - lon0)) / G);
            alphac = StrictMath.asin(D * StrictMath.sin(gamma0));
        } else {
            double F = D_P_sqrt_DD_M_one;
            double G = (F - 1.0 / F) * 0.5;

            gamma0 = StrictMath.asin(StrictMath.sin(alphac) / D);
            lon0   = lonc - StrictMath.asin(G * StrictMath.tan(gamma0)) / B;
        }

        double cosgamma0 = StrictMath.cos(gamma0);
        double singamma0 = StrictMath.sin(gamma0);

        double cosalphac = StrictMath.cos(alphac);
        double sinalphac = StrictMath.sin(alphac);

        double v, u, Qp, Sp, Tp, Vp, Up, t;

        double phi, phi_old, esinlat;
        int j;

        //
        // TODO add special case checks as defined on pg 71 paragraph following Alternate A givens.
        //
        for(int i = 0; i < lon.length; ++i) {
            v = (x[i] - x0) * cosalphac - (y[i] - y0) * sinalphac;
            u = (y[i] - y0) * cosalphac + (x[i] - x0) * sinalphac;

            Qp = StrictMath.exp(-(B * v / A));
            Sp = (Qp - 1.0 / Qp) * 0.5;
            Tp = (Qp + 1.0 / Qp) * 0.5;
            Vp = StrictMath.sin(B * u / A);
            Up = (Vp * cosgamma0 + Sp * singamma0) / Tp;
            t  = StrictMath.pow(E / StrictMath.sqrt((1.0 + Up) / (1.0 - Up)), one_over_B);

            if(1.0 - StrictMath.abs(Up) < NEAR_ZERO_DEG) {
                lon[i] = lon0;
                lat[i] = StrictMath.copySign(PI_DIV_2, Up);
            }
            else {
                lon[i] = normalizeLonRad(
                        lon0 - 
                                normalizeLonRad(StrictMath.atan2(Sp * cosgamma0 - Vp * singamma0, StrictMath.cos(B * u / A))
                                 ) / B
                        );

                // (7-9)
                phi = PI_DIV_2 - 2.0 * StrictMath.atan(t);
                j = 0;
                do {
                    phi_old = phi;
                    esinlat = e * StrictMath.sin(phi);
                    phi = PI_DIV_2 - 2.0 * StrictMath.atan(t * StrictMath.pow((1.0 - esinlat) / (1.0 + esinlat), hlf_e));
                } while (++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phi - phi_old));

                lat[i] = phi;
            }
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double e   = ellip.getProperty("e");
        double esq = ellip.getProperty("e^2");

        double k0   = datum.getProperty("k0");
        double lat0 = datum.getProperty("lat0");
        double x0   = datum.getProperty("x0");
        double y0   = datum.getProperty("y0");

        // Used for Central Line
        double lon1 = datum.getProperty("lon1");
        double lat1 = datum.getProperty("lat1");
        double lon2 = datum.getProperty("lon2");
        double lat2 = datum.getProperty("lat2");

        // Used for Central Point
        double lonc   = datum.getProperty("lonc");
        double alphac = datum.getProperty("alphac");

        boolean useCentralPoint = datum.getProperty("useCentralPoint") == 1.0;

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double hlf_e = e * 0.5;

        double one_M_esq = 1.0 - esq;
        double sqrt_one_M_esq = StrictMath.sqrt(one_M_esq);

        double coslat0 = StrictMath.cos(lat0);
        double sinlat0 = StrictMath.sin(lat0);

        double esinlat0                = e * sinlat0;
        double one_M_esqsinlat0sinlat0 = 1.0 - esq * sinlat0 * sinlat0;
        double coslat0_pow_4           = coslat0 * coslat0 * coslat0 * coslat0;

        // Common variable for both projection alternatives
        double B  = StrictMath.sqrt(1.0 + esq * coslat0_pow_4 / one_M_esq);
        double A  = a * B * k0 * sqrt_one_M_esq / one_M_esqsinlat0sinlat0;
        double t0 = StrictMath.tan(PI_DIV_4 - (lat0 * 0.5)) / StrictMath.pow((1.0 - esinlat0) / (1.0 + esinlat0), hlf_e);
        double D  = B * sqrt_one_M_esq / (coslat0 * StrictMath.sqrt(one_M_esqsinlat0sinlat0));
        double DD = ((DD = D * D) < 1.0 ? 1.0 : DD);

        double D_P_sqrt_DD_M_one = D + StrictMath.copySign(StrictMath.sqrt(DD - 1.0), lat0);
        double E                 = D_P_sqrt_DD_M_one * StrictMath.pow(t0, B);

        double
        B2       = 2.0 * B,
        AB       = A * B,
        A_over_B = A / B;

        double gamma0, lon0, v = 0.0, u = 0.0;

        if(useCentralPoint == false) {
            double esinlat1 = e * StrictMath.sin(lat1);
            double esinlat2 = e * StrictMath.sin(lat2);

            double t1 = StrictMath.tan(PI_DIV_4 - (lat1 * 0.5)) / StrictMath.pow((1.0 - esinlat1) / (1.0 + esinlat1), hlf_e);
            double t2 = StrictMath.tan(PI_DIV_4 - (lat2 * 0.5)) / StrictMath.pow((1.0 - esinlat2) / (1.0 + esinlat2), hlf_e);

            double H = StrictMath.pow(t1, B);
            double L = StrictMath.pow(t2, B);
            double F = E / H;
            double G = (F - 1.0 / F) * 0.5;
            double J = (E * E - L * H) / (E * E + L * H);
            double P = (L - H) / (L + H);

            lon0   = (lon1 + lon2) * 0.5 - StrictMath.atan(J * StrictMath.tan(B * normalizeLonRad(lon1 - lon2) * 0.5) / P) / B;
            gamma0 = StrictMath.atan(StrictMath.sin(B * (lon1 - lon0)) / G);
            alphac = StrictMath.asin(D * StrictMath.sin(gamma0));
        } else {
            double F = D_P_sqrt_DD_M_one;
            double G = (F - 1.0 / F) * 0.5;

            gamma0 = StrictMath.asin(StrictMath.sin(alphac) / D);
            lon0   = lonc - StrictMath.asin(G * StrictMath.tan(gamma0)) / B;

            u = StrictMath.copySign(
            		A_over_B * StrictMath.atan2(StrictMath.sqrt(DD - 1.0),
            		StrictMath.cos(alphac)), lat0);
        }

        double cosgamma0 = StrictMath.cos(gamma0);
        double singamma0 = StrictMath.sin(gamma0);

        double cosalphac = StrictMath.cos(alphac);
        double sinalphac = StrictMath.sin(alphac);

        double t, esinlat, Q, S, T, V, U, lon_M_lon0, cosBlon_M_lon0;

        //
        // TODO add special case checks as defined on pg 71 paragraph following Alternate A givens.
        //
        for(int i = 0; i < lon.length; ++i) {
        	if(PI_DIV_2 - StrictMath.abs(lat[i]) < NEAR_ZERO_RAD) {
        		v = A_over_B * StrictMath.log(StrictMath.tan(PI_DIV_4 + gamma0 * StrictMath.copySign(0.5, lat[i])));
        		u = A * lat[i] / B;
        	} else {
            	lon_M_lon0 = lon[i] - lon0;
            	esinlat = e * StrictMath.sin(lat[i]);
            	cosBlon_M_lon0 = StrictMath.cos(B * lon_M_lon0);

            	t = StrictMath.tan(PI_DIV_4 - lat[i] * 0.5) / StrictMath.pow((1.0 - esinlat) / (1.0 + esinlat), hlf_e);
        		Q = E / StrictMath.pow(t, B);
        		S = (Q - 1.0 / Q) * 0.5;
        		T = (Q + 1.0 / Q) * 0.5;
        		V = StrictMath.sin(B * lon_M_lon0);
        		U = (-V * cosgamma0 + S * singamma0) / T;
        		v = A * StrictMath.log((1.0 - U) / (1.0 + U)) / B2;

        		if(StrictMath.abs(cosBlon_M_lon0) < NEAR_ZERO_RAD)
        			u = AB * lon_M_lon0;
        		else
        			u = A * StrictMath.atan2(S * cosgamma0 + V * singamma0, cosBlon_M_lon0) / B;
        	}

        	x[i] = v * cosalphac + u * sinalphac + x0;
        	y[i] = u * cosalphac - v * sinalphac + y0;
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "k0",
                "lat0",
                "x0",
                "y0",
                "lon1",
                "lat1",
                "lon2",
                "lat2",
                "lonc",
                "alphac",
                "useCentralPoint"}));
    }

}