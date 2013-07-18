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

import gis.proj.Azimuthal;
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
 * Pages 154-163
 *
 */

public final class Stereographic implements Azimuthal {

    public String getName() {
        return "Stereographic";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a     = ellip.getProperty("a");
        double e     = ellip.getProperty("e");

        double k0   = datum.getProperty("k0");
        double lon0 = datum.getProperty("lon0");
        double lat1 = datum.getProperty("lat1");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double hlf_e = e * 0.5;

        double qA = 2.0 * a * k0;

        double phi, phi_old, rho, esinphi;
        int j;

        // Polar aspect
        if (POLE_RAD - StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
            double t, xt, yt;

            double o_P_e = 1.0 + e;
            double o_M_e = 1.0 - e;

            double sqrt_ope_M_ome = StrictMath.sqrt(
                    StrictMath.pow(o_P_e, o_P_e) *
                    StrictMath.pow(o_M_e, o_M_e));

            double poleFactor = lat1 < 0.0 ? -1 : 1.0;
            lon0 *= poleFactor;
            lat1 *= poleFactor;

            for(int i = 0; i < x.length; ++i) {
                xt = x[i] * poleFactor;
                yt = y[i] * poleFactor;

                rho = StrictMath.hypot(xt, yt);

                t = (rho * sqrt_ope_M_ome) / qA;

                j = 0;
                phi = PI_DIV_2 - 2.0 * StrictMath.atan(t);
                do {
                    phi_old = phi;
                    esinphi = e * StrictMath.sin(phi);
                    phi = PI_DIV_2 - 2.0 * StrictMath.atan(
                            t * StrictMath.pow(
                                    (1.0 - esinphi) /
                                    (1.0 + esinphi),
                                    hlf_e)
                            );
                } while (++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phi - phi_old));

                lon[i] = normalizeLonRad(lon0 + StrictMath.atan2(xt, -yt)) * poleFactor;
                lat[i] = phi * poleFactor;
            }
        // Oblique et Equatorial aspects
        } else {
            double Ce, chi, sinCe, cosCe;

            double sinlat1 = StrictMath.sin(lat1);
            double esinlat1 = e * sinlat1;

            double m1 = StrictMath.cos(lat1) /
                    StrictMath.sqrt(1.0 - ellip.getProperty("e^2") * (sinlat1 * sinlat1));

            double qAm1 = qA * m1;

            double chi1 = 2.0 *
                    (
                            StrictMath.atan(
                                    StrictMath.tan(PI_DIV_4 + lat1 * 0.5) *
                                    StrictMath.pow(
                                            (1.0 - esinlat1) /
                                            (1.0 + esinlat1),
                                            hlf_e)
                                    )
                            ) + N_PI_DIV_2;

            double cosX1 = StrictMath.cos(chi1);
            double sinX1 = StrictMath.sin(chi1);

            for(int i = 0; i < x.length; ++i) {
                rho = StrictMath.hypot(x[i], y[i]);

                Ce = 2.0 * StrictMath.atan((rho * cosX1) / qAm1);

                sinCe = StrictMath.sin(Ce);
                cosCe = StrictMath.cos(Ce);

                chi = StrictMath.asin(cosCe * sinX1 + (y[i] * sinCe * cosX1) / rho);

                if(NEAR_ZERO_DEG <= rho)
                    lon[i] = lon0 + StrictMath.atan2((x[i] * sinCe), (rho * cosX1 * cosCe - y[i] * sinX1 * sinCe));
                else {
                    chi = chi1;
                    lon[i] = lon0;
                }

                j = 0;
                phi = chi;
                do {
                    phi_old = phi;
                    esinphi = e * StrictMath.sin(phi);
                    phi = (2.0 * StrictMath.atan(
                            StrictMath.tan(PI_DIV_4 + chi * 0.5) *
                            StrictMath.pow(
                                    (1.0 + esinphi) /
                                    (1.0 - esinphi),
                                    hlf_e)
                            )) + N_PI_DIV_2;
                } while (++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phi - phi_old));

                lon[i] = normalizeLonRad(lon[i]);
                lat[i] = phi;
            }
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a = ellip.getProperty("a");
        double e = ellip.getProperty("e");

        double k0   = datum.getProperty("k0");
        double lon0 = datum.getProperty("lon0");
        double lat1 = datum.getProperty("lat1");

        double x[] = new double[lon.length];
        double y[] = new double[lat.length];

        double hlf_e = e * 0.5;

        double A = 0.0;
        double qA = 2.0 * a * k0;

        double chi, lon_M_lon0, esinlat, coschi;

        // Equatorial aspect
        if(StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
            for(int i = 0; i < lon.length; ++i) {

                esinlat = e * StrictMath.sin(lat[i]);

                chi = 2.0 *
                        (
                                StrictMath.atan(
                                        StrictMath.tan(PI_DIV_4 + lat[i] * 0.5) *
                                        StrictMath.pow(
                                                (1.0 - esinlat) /
                                                (1.0 + esinlat),
                                                hlf_e)
                                        )
                                ) + N_PI_DIV_2;

                coschi = StrictMath.cos(chi);

                lon_M_lon0 = lon[i] - lon0;

                A = qA / (1.0 + coschi * StrictMath.cos(lon_M_lon0));

                x[i] = A * coschi * StrictMath.sin(lon_M_lon0);
                y[i] = A * StrictMath.sin(chi);
            }
        }
        // Polar aspect
        else if (POLE_RAD - StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
            double t, rho, latt, lont;

            double o_P_e = 1.0 + e;
            double o_M_e = 1.0 - e;

            double sqrt_ope_M_ome = StrictMath.sqrt(
                    StrictMath.pow(o_P_e, o_P_e) *
                    StrictMath.pow(o_M_e, o_M_e));

            double poleFactor = lat1 < 0.0 ? -1.0 : 1.0;
            lon0 *= poleFactor;
            lat1 *= poleFactor;

            for(int i = 0; i < lon.length; ++i) {

                lont = lon[i] * poleFactor;
                latt = lat[i] * poleFactor;

                esinlat = e * StrictMath.sin(latt);

                t = StrictMath.tan(PI_DIV_4 - latt * 0.5) /
                        StrictMath.pow(
                                (1.0 - esinlat) /
                                (1.0 + esinlat)
                                , hlf_e);

                rho = (qA * t) / sqrt_ope_M_ome;

                lon_M_lon0 = lont - lon0;

                x[i] =  rho * StrictMath.sin(lon_M_lon0) * poleFactor;
                y[i] = -rho * StrictMath.cos(lon_M_lon0) * poleFactor;
            }
        }
        // Oblique aspect
        else {
            double sinlat1 = StrictMath.sin(lat1);
            double esinlat1 = e * sinlat1;

            double m1 = StrictMath.cos(lat1) /
                    StrictMath.sqrt(1.0 - ellip.getProperty("e^2") * (sinlat1 * sinlat1));

            double qAm1 = qA * m1;

            double chi1 = 2.0 *
                    (
                            StrictMath.atan(
                                    StrictMath.tan(PI_DIV_4 + lat1 * 0.5) *
                                    StrictMath.pow(
                                            (1.0 - esinlat1) /
                                            (1.0 + esinlat1),
                                            hlf_e)
                                    )
                            ) + N_PI_DIV_2;

            double cosX1 = StrictMath.cos(chi1);
            double sinX1 = StrictMath.sin(chi1);

            double sinchi, coslon_M_lon0;

            for(int i = 0; i < lon.length; ++i) {

                esinlat = e * StrictMath.sin(lat[i]);

                chi = 2.0 *
                        (
                                StrictMath.atan(
                                        StrictMath.tan(PI_DIV_4 + lat[i] * 0.5) *
                                        StrictMath.pow(
                                                (1.0 - esinlat) /
                                                (1.0 + esinlat),
                                                hlf_e)
                                        )
                                ) + N_PI_DIV_2;

                coschi = StrictMath.cos(chi);
                sinchi = StrictMath.sin(chi);

                lon_M_lon0 = lon[i] - lon0;

                coslon_M_lon0 = StrictMath.cos(lon_M_lon0);

                A = qAm1 /
                        (
                                cosX1 *
                                (
                                        1.0 +
                                        sinX1 * sinchi +
                                        cosX1 * coschi * coslon_M_lon0
                                        )
                                );

                x[i] = A * coschi * StrictMath.sin(lon_M_lon0);
                y[i] = A * (
                        cosX1 * sinchi -
                        sinX1 * coschi * coslon_M_lon0
                        );
            }
        }

        return new double[][] {x, y};
    }
 
    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "k0",
                "lon0",
                "lat1"}));
    }

}