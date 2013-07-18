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
 * Pages 182-190
 *
 */

public final class LambertAzimuthalEA implements Azimuthal {

    public String getName() {
        return "Lambert Azimuthal Equal-Area";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double e   = ellip.getProperty("e");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat1 = datum.getProperty("lat1");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        if(ellip.isSphere()) {
            double coslat1 = StrictMath.cos(lat1);
            double sinlat1 = StrictMath.sin(lat1);

            double a2 = a * 2.0;

            double rho, c, cosc, sinc;

            // Polar aspect
            if (POLE_RAD - StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
                double poleFactor = lat1 < 0.0 ? 1.0 : -1.0;

                for(int i = 0; i < x.length; ++i) {
                    rho = StrictMath.hypot(x[i], y[i]);
                    c = 2.0 * StrictMath.asin(rho / a2);

                    if(rho < NEAR_ZERO_DEG) {
                        lon[i] = lon0;
                        lat[i] = lat1;
                    }
                    else {
                        cosc = StrictMath.cos(c);
                        sinc = StrictMath.sin(c);

                        lon[i] = normalizeLonRad(lon0 + StrictMath.atan2(x[i], y[i] * poleFactor));
                        lat[i] = StrictMath.asin(cosc * sinlat1 + (y[i] * sinc * coslat1 / rho));
                    }
                }
            }
            // Oblique et Equatorial aspects
            else {
                for(int i = 0; i < x.length; ++i) {
                    rho = StrictMath.hypot(x[i], y[i]);
                    c = 2.0 * StrictMath.asin(rho / a2);

                    if(rho < NEAR_ZERO_DEG) {
                        lon[i] = lon0;
                        lat[i] = lat1;
                    }
                    else {
                        cosc = StrictMath.cos(c);
                        sinc = StrictMath.sin(c);

                        lon[i] = normalizeLonRad(
                                lon0 + StrictMath.atan2(
                                        x[i] * sinc,
                                        rho * coslat1 * cosc - y[i] * sinlat1 * sinc
                                        )
                                );

                        lat[i] = StrictMath.asin(cosc * sinlat1 + (y[i] * sinc * coslat1 / rho));
                    }
                }
            }
        } else {
            double one_M_esq = 1.0 - esq;
            double one_over_2e = 1.0 / (2.0 * e);

            double sinlat1 = StrictMath.sin(lat1);
            double esinlat1 = e * sinlat1;
            double esqsinlat1sinlat1 = esq * sinlat1 * sinlat1;

            double m1 = StrictMath.cos(lat1) /
                    StrictMath.sqrt(1.0 - esqsinlat1sinlat1);

            double q =
                    one_M_esq * (sinlat1 / (1.0 - esqsinlat1sinlat1) - one_over_2e * StrictMath.log((1.0 - esinlat1) / (1.0 + esinlat1))
                            );

            double sin_hlf_pi = StrictMath.sin(PI_DIV_2);
            double esin_hlf_pi = e * sin_hlf_pi;

            double qp = one_M_esq * (
                    sin_hlf_pi / (1.0 - esq * sin_hlf_pi * sin_hlf_pi) - one_over_2e * StrictMath.log((1.0 - esin_hlf_pi) / (1.0 + esin_hlf_pi))
                    );

            double phi_old, phi, factor, sinphi, esinphi;
            int j;

            double rho_over_a, rho;

            // Polar aspect
            if (POLE_RAD - StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
                double poleFactor = lat1 < 0.0 ? 1.0 : -1.0;

                for(int i = 0; i < x.length; ++i) {
                    lon[i] = normalizeLonRad(lon0 + StrictMath.atan2(x[i], (y[i] * poleFactor)));

                    rho = StrictMath.hypot(x[i], y[i]);
                    rho_over_a = rho / a;

                    q = qp - rho_over_a * rho_over_a;

                    phi = q * 0.5;
                    j = 0;
                    do {
                        phi_old = phi;
                        sinphi = StrictMath.sin(phi);
                        esinphi = e * sinphi;
                        // (3-16)
                        factor = (1.0 - esq * sinphi * sinphi);
                        phi = phi + ((factor * factor) / (2.0 * StrictMath.cos(phi))) * (
                                q / one_M_esq - sinphi / factor + one_over_2e * StrictMath.log((1.0 - esinphi) / (1.0 + esinphi))
                                );
                    } while (++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phi - phi_old));

                    lat[i] = phi * -poleFactor;
                }
            // Oblique et Equatorial aspects
            } else {
                // (3-11)
                double beta1 = StrictMath.asin(q / qp);

                double cosbeta1 = StrictMath.cos(beta1);
                double sinbeta1 = StrictMath.sin(beta1);

                // (3-13)
                double Rq = a * StrictMath.sqrt(qp * 0.5);
                double D  = a * m1 / (Rq * cosbeta1);
                double DD = D * D;

                double Ce, cosCe, sinCe, Dy;

                for(int i = 0; i < x.length; ++i) {
                    Dy = D * y[i];
                    rho = StrictMath.hypot(x[i] / D, Dy);

                    if(NEAR_ZERO_DEG < rho) {
                        //
                        // Snyder's text implies that (24-29) be interpreted as
                        // Ce = 2 * asin((rho * Rq) / 2); however, that is incorrect
                        // and will result in NaN for Ce.  The correct equation is
                        // below as Ce = 2. * asin(rho / (2 * Rq))
                        //
                        Ce = 2.0 * StrictMath.asin(rho * 0.5 / Rq);
                        cosCe = StrictMath.cos(Ce);
                        sinCe = StrictMath.sin(Ce);

                        q = qp * (cosCe * sinbeta1 + (Dy * sinCe * cosbeta1 / rho));

                        lon[i] = normalizeLonRad(lon0 +
                                StrictMath.atan2(
                                        x[i] * sinCe,
                                        D * rho * cosbeta1 * cosCe - DD * y[i] * sinbeta1 * sinCe)
                                );
                    }
                    else {
                        q = qp * sinbeta1;
                        lon[i] = lon0;
                    }

                    phi = q * 0.5;
                    j = 0;
                    do {
                        phi_old = phi;
                        sinphi = StrictMath.sin(phi);
                        esinphi = e * sinphi;
                        // (3-16)
                        factor = (1.0 - esq * sinphi * sinphi);
                        phi = phi + ((factor * factor) / (2.0 * StrictMath.cos(phi))) * (
                                q / one_M_esq - sinphi / factor + one_over_2e * StrictMath.log((1.0 - esinphi) / (1.0 + esinphi))
                                );
                    } while (++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phi - phi_old));

                    lat[i] = phi;
                }
            }
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double e   = ellip.getProperty("e");
        double esq = ellip.getProperty("e^2");

        double lon0 = datum.getProperty("lon0");
        double lat1 = datum.getProperty("lat1");

        double x[] = new double[lon.length];
        double y[] = new double[lat.length];

        if(ellip.isSphere()) {
            // Equatorial aspect
            if(StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
                double coslat, lon_M_lon0, kp;

                for(int i = 0; i < lon.length; ++i) {
                    coslat = StrictMath.cos(lat[i]);

                    lon_M_lon0 = normalizeLonRad(lon[i] - lon0);

                    kp = StrictMath.sqrt(2.0 / (1.0 + coslat * StrictMath.cos(lon_M_lon0)));

                    x[i] = a * kp * coslat * StrictMath.sin(lon_M_lon0);
                    y[i] = a * kp * StrictMath.sin(lat[i]);
                }
            }
            // Polar aspect
            else if (POLE_RAD - StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
                double a2 = a * 2.0;

                double lon_M_lon0, QP_M_Hlat;

                if(lat1 < 0.0) {
                    for(int i = 0; i < lon.length; ++i) {
                        QP_M_Hlat = StrictMath.cos(PI_DIV_4 - lat[i] * 0.5);
                        lon_M_lon0 = normalizeLonRad(lon[i] - lon0);

                        x[i] = a2 * QP_M_Hlat * StrictMath.sin(lon_M_lon0);
                        y[i] = a2 * QP_M_Hlat * StrictMath.cos(lon_M_lon0);
                    }
                } else {
                    for(int i = 0; i < lon.length; ++i) {
                        QP_M_Hlat = StrictMath.sin(PI_DIV_4 - lat[i] * 0.5);
                        lon_M_lon0 = normalizeLonRad(lon[i] - lon0);

                        x[i] =  a2 * QP_M_Hlat * StrictMath.sin(lon_M_lon0);
                        y[i] = -a2 * QP_M_Hlat * StrictMath.cos(lon_M_lon0);
                    }
                }
            }
            // Oblique aspect
            else {
                double coslat1 = StrictMath.cos(lat1);
                double sinlat1 = StrictMath.sin(lat1);

                double coslat, sinlat, coslon_M_lon0, lon_M_lon0, kp;

                for(int i = 0; i < lon.length; ++i) {
                    coslat = StrictMath.cos(lat[i]);
                    sinlat = StrictMath.sin(lat[i]);

                    lon_M_lon0 = normalizeLonRad(lon[i] - lon0);
                    coslon_M_lon0 = StrictMath.cos(lon_M_lon0);

                    kp = StrictMath.sqrt(2.0 / (1.0 + sinlat1 * sinlat + coslat1 * coslat * coslon_M_lon0));

                    x[i] = a * kp * coslat * StrictMath.sin(lon_M_lon0);
                    y[i] = a * kp * (coslat1 * sinlat - sinlat1 * coslat * coslon_M_lon0);
                }
            }
        } else {
            double one_M_esq = 1.0 - esq;
            double one_over_2e = 1.0 / (2.0 * e);

            double sin_hlf_pi = StrictMath.sin(PI_DIV_2);
            double esin_hlf_pi = e * sin_hlf_pi;

            double qp = one_M_esq * (
                    sin_hlf_pi / (1.0 - esq * sin_hlf_pi * sin_hlf_pi) - one_over_2e * StrictMath.log((1.0 - esin_hlf_pi) / (1.0 + esin_hlf_pi))
                    );

            // (3-13)
            double Rq = a * StrictMath.sqrt(qp * 0.5);

            double lon_M_lon0, cosbeta, sinbeta, beta, sinlat, esinlat, esqsinlatsinlat;

            // Equatorial aspect
            if(StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
                double xy_factor, q;

                for(int i = 0; i < lon.length; ++i) {
                    lon_M_lon0 = normalizeLonRad(lon[i] - lon0);
                    // (3-12)
                    sinlat = StrictMath.sin(lat[i]);
                    esinlat = e * sinlat;
                    esqsinlatsinlat = esq * sinlat * sinlat;

                    q = one_M_esq * (
                            sinlat / (1.0 - esqsinlatsinlat) - one_over_2e * StrictMath.log((1.0 - esinlat) / (1.0 + esinlat))
                            );

                    // (3-11)
                    beta = StrictMath.asin(q / qp);

                    cosbeta  = StrictMath.cos(beta);
                    sinbeta  = StrictMath.sin(beta);

                    xy_factor = StrictMath.sqrt(
                            2.0 / (1.0 + cosbeta * StrictMath.cos(lon_M_lon0))
                            );

                    x[i] = a * cosbeta * StrictMath.sin(lon_M_lon0) * xy_factor;
                    y[i] = ((Rq * Rq) / a) * sinbeta * xy_factor;
                }
            }
            // Polar aspect
            else if (POLE_RAD - StrictMath.abs(lat1) < NEAR_ZERO_RAD) {
                double poleFactor = lat1 < 0.0 ? 1.0 : -1.0;

                double rho, q;

                for(int i = 0; i < lon.length; ++i) {
                    lon_M_lon0 = normalizeLonRad(lon[i] - lon0);

                    // (3-12)
                    sinlat = StrictMath.sin(lat[i]);
                    esinlat = e * sinlat;
                    esqsinlatsinlat = esq * sinlat * sinlat;

                    q = one_M_esq * (
                            sinlat / (1.0 - esqsinlatsinlat) - one_over_2e * StrictMath.log((1.0 - esinlat) / (1.0 + esinlat))
                            );

                    rho = a * StrictMath.sqrt(qp + q * poleFactor);

                    x[i] = rho * StrictMath.sin(lon_M_lon0);
                    y[i] = rho * StrictMath.cos(lon_M_lon0) * poleFactor;
                }
            }
            // Oblique aspect
            else {
                double sinlat1 = StrictMath.sin(lat1);
                double esinlat1 = e * sinlat1;
                double esqsinlat1sinlat1 = esq * sinlat1 * sinlat1;

                double m1 = StrictMath.cos(lat1) /
                        StrictMath.sqrt(1.0 - esqsinlat1sinlat1);

                double q =
                        one_M_esq * (sinlat1 / (1.0 - esqsinlat1sinlat1) - one_over_2e * StrictMath.log((1.0 - esinlat1) / (1.0 + esinlat1))
                                );

                // (3-11)
                double beta1 = StrictMath.asin(q / qp);

                double cosbeta1 = StrictMath.cos(beta1);
                double sinbeta1 = StrictMath.sin(beta1);

                double D  = a * m1 / (Rq * cosbeta1);

                double B, coslon_M_lon0;

                for(int i = 0; i < lon.length; ++i) {
                    lon_M_lon0 = normalizeLonRad(lon[i] - lon0);
                    coslon_M_lon0 = StrictMath.cos(lon_M_lon0);

                    // (3-12)
                    sinlat = StrictMath.sin(lat[i]);
                    esinlat = e * sinlat;
                    esqsinlatsinlat = esq * sinlat * sinlat;

                    q = one_M_esq * (
                            sinlat / (1.0 - esqsinlatsinlat) - one_over_2e * StrictMath.log((1.0 - esinlat) / (1.0 + esinlat))
                            );

                    // (3-11)
                    beta = StrictMath.asin(q / qp);

                    cosbeta  = StrictMath.cos(beta);
                    sinbeta  = StrictMath.sin(beta);

                    B  = Rq * StrictMath.sqrt(2.0 / (1.0 + sinbeta1 * sinbeta + cosbeta1 * cosbeta * coslon_M_lon0));

                    x[i] = B * D * cosbeta * StrictMath.sin(lon_M_lon0);
                    y[i] = (B / D) * (cosbeta1 * sinbeta - sinbeta1 * cosbeta * coslon_M_lon0);
                }
            }
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "lon0",
                "lat1"}));

    }

}