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

import gis.proj.Datum;
import gis.proj.Ellipsoid;
import gis.proj.Space;
import gis.proj.Spherical;

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
 * Pages 230-238
 *
 */

public final class SatelliteTracking implements Space, Spherical {

    public String getName() {
        return "Satellite Tracking";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double R = ellip.getProperty("R");

        double _i   = datum.getProperty("i");
        double P2   = datum.getProperty("P2");
        double P1   = datum.getProperty("P1");
        double lon0 = datum.getProperty("lon0");
        double lat0 = datum.getProperty("lat0");
        double lat1 = datum.getProperty("lat1");
        double lat2 = datum.getProperty("lat2");

        boolean useConic = datum.getProperty("useConic") == 1.0;

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double coslat1        = StrictMath.cos(lat1);
        double cosi           = StrictMath.cos(_i);
        double sini           = StrictMath.sin(_i);
        double coslat1coslat1 = coslat1 * coslat1;
        double P2_div_P1      = P2 / P1;

        if(useConic == false) {
            double F1p = Fx(P2_div_P1, coslat1coslat1, cosi);

            double Rcoslat1 = R * coslat1;
            double F1p_div_coslat1 = F1p / Rcoslat1;

            double L, lonp;
            double A, delta_lonp, AA;
            double cosicosi = cosi * cosi;
            int j;

            for(int i = 0; i < x.length; ++i) {
                L = y[i] * F1p_div_coslat1;

                lonp = N_POLE_RAD;
                j = 0;
                do {
                    A = StrictMath.tan(L + P2_div_P1 * lonp) / cosi;
                    AA = A * A;
                    delta_lonp = -(lonp - StrictMath.atan(A)) / (1.0 - (AA + 1.0 / cosicosi) * P2_div_P1 * cosi / (AA + 1.0));
                    lonp += delta_lonp;
                } while(++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(delta_lonp));

                lon[i] = normalizeLonRad(lon0 + x[i] / Rcoslat1);
                lat[i] = -StrictMath.asin(StrictMath.sin(lonp) * sini);
            }
        }
        else {
            double coslat2 = StrictMath.cos(lat2);
            double
            F1 = StrictMath.atan(Fx(P2_div_P1, coslat1coslat1, cosi)),
            F2 = StrictMath.atan(Fx(P2_div_P1, coslat2 * coslat2, cosi));

            double
            lonp0 = lonpx(lat0, sini),
            lonp1 = lonpx(lat1, sini),
            lonp2 = lonpx(lat2, sini);

            double
            L0 = lontx(lonp0, cosi) - P2_div_P1 * lonp0,
            L1 = lontx(lonp1, cosi) - P2_div_P1 * lonp1,
            L2 = lontx(lonp2, cosi) - P2_div_P1 * lonp2;

            double n = (F2 - F1) / (L2 - L1);

            // Used for one standard parallel
            if(StrictMath.abs(lat1 - lat2) < NEAR_ZERO_RAD) {
                n = StrictMath.sin(lat1) * (P2_div_P1 * (2.0 * cosi * cosi - coslat1coslat1) - cosi) / (
                        (P2_div_P1 * coslat1coslat1 - cosi) * (P2_div_P1 * (P2_div_P1 * coslat1coslat1 - 2.0 * cosi) + 1.0)
                        );
            }
            // TODO add tracking limits and implement the follow:
            //
            // If one std parallel (lat1) and lat1 == the upper tracking limit, then
            //      n = sini / StrictMath.pow(P2_div_P1 * cosi - 1.0, 2.0); // (28-18)

            double rho_num = R * coslat1 * StrictMath.sin(F1);

            double s0   = F1 - n * L1;
            double rho0 = rho_num / (n * StrictMath.sin(n * L0 + s0));

            double lonp, L, rho, theta;

            double rho0_M_y;

            double A, delta_lonp, AA;
            double cosicosi = cosi * cosi;
            int j;

            for(int i = 0; i < lon.length; ++i) {
                rho0_M_y = rho0 - y[i];

                rho = StrictMath.copySign(StrictMath.hypot(x[i], rho0_M_y), n);
                theta = StrictMath.atan2(x[i], rho0_M_y);
                L = (StrictMath.asin(rho_num / (rho * n)) - s0) / n;

                lonp = N_POLE_RAD;
                j = 0;
                do {
                    A = StrictMath.tan(L + P2_div_P1 * lonp) / cosi;
                    AA = A * A;
                    delta_lonp = -(lonp - StrictMath.atan(A)) / (1.0 - (AA + 1.0 / cosicosi) * P2_div_P1 * cosi / (AA + 1.0));
                    lonp += delta_lonp;
                } while(++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(delta_lonp));

                lon[i] = normalizeLonRad(lon0 + theta / n);
                lat[i] = -StrictMath.asin(StrictMath.sin(lonp) * sini);
            }
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double R = ellip.getProperty("R");

        double _i   = datum.getProperty("i");
        double P2   = datum.getProperty("P2");
        double P1   = datum.getProperty("P1");
        double lon0 = datum.getProperty("lon0");
        double lat0 = datum.getProperty("lat0");
        double lat1 = datum.getProperty("lat1");
        double lat2 = datum.getProperty("lat2");

        boolean useConic = datum.getProperty("useConic") == 1.0;

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double coslat1        = StrictMath.cos(lat1);
        double cosi           = StrictMath.cos(_i);
        double sini           = StrictMath.sin(_i);
        double coslat1coslat1 = coslat1 * coslat1;
        double P2_div_P1      = P2 / P1;

        if(useConic == false) {
            double F1p = Fx(P2_div_P1, coslat1coslat1, cosi);

            double Rcoslat1 = R * coslat1;
            double Rcoslat1_div_F1p = Rcoslat1 / F1p;

            double lon_M_lon0, lonp, lont, L;

            for(int i = 0; i < lon.length; ++i) {
                lon_M_lon0 = normalizeLonRad(lon[i] - lon0);

                lonp = -StrictMath.asin(StrictMath.sin(lat[i]) / sini);
                lont =  StrictMath.atan(StrictMath.tan(lonp) * cosi);
                L = lont - P2_div_P1 * lonp;

                //        		if(NEAR_ZERO_RAD < (81 * DEG_TO_RAD) - StrictMath.abs(lat[i])) {
                x[i] = lon_M_lon0 * Rcoslat1;
                y[i] = L * Rcoslat1_div_F1p;
                //        		}
                //        		else {
                //        		    x[i] = y[i] = Double.NEGATIVE_INFINITY;
                //        		}
            }
        }
        else {
            double coslat2 = StrictMath.cos(lat2);

            double
            F1 = StrictMath.atan(Fx(P2_div_P1, coslat1coslat1, cosi)),
            F2 = StrictMath.atan(Fx(P2_div_P1, coslat2 * coslat2, cosi));

            double
            lonp0 = lonpx(lat0, sini),
            lonp1 = lonpx(lat1, sini),
            lonp2 = lonpx(lat2, sini);

            double
            L0 = lontx(lonp0, cosi) - P2_div_P1 * lonp0,
            L1 = lontx(lonp1, cosi) - P2_div_P1 * lonp1,
            L2 = lontx(lonp2, cosi) - P2_div_P1 * lonp2;

            double n = (F2 - F1) / (L2 - L1);

            // Used for one standard parallel
            if(StrictMath.abs(lat1 - lat2) < NEAR_ZERO_RAD) {
                n = StrictMath.sin(lat1) * (P2_div_P1 * (2.0 * cosi * cosi - coslat1coslat1) - cosi) / (
                        (P2_div_P1 * coslat1coslat1 - cosi) * (P2_div_P1 * (P2_div_P1 * coslat1coslat1 - 2.0 * cosi) + 1.0)
                        );
            }
            // TODO add tracking limits and implement the follow:
            //
            // If one std parallel (lat1) and lat1 == the upper tracking limit, then
            //      n = sini / StrictMath.pow(P2_div_P1 * cosi - 1.0, 2.0); // (28-18)

            double rho_num = R * coslat1 * StrictMath.sin(F1);

            double s0   = F1 - n * L1;
            double rho0 = rho_num / (n * StrictMath.sin(n * L0 + s0));

            double N_s0_div_n = -s0 / n;

            double lonp, lont, L, rho, theta;

            for(int i = 0; i < lon.length; ++i) {
                lonp = -StrictMath.asin(StrictMath.sin(lat[i]) / sini);
                lont = StrictMath.atan(StrictMath.tan(lonp) * cosi);
                L = lont - P2_div_P1 * lonp;

                if(n < 0.0 && L < N_s0_div_n || n < 0.0 == false && L > N_s0_div_n) {
                    rho = rho_num / (n * StrictMath.sin(n * L + s0));
                    theta = n * normalizeLonRad(lon[i] - lon0);

                    x[i] = rho * StrictMath.sin(theta);
                    y[i] = rho0 - rho * StrictMath.cos(theta);
                }
                else {
                    x[i] = y[i] = Double.NEGATIVE_INFINITY;
                }
            }
        }

        return new double[][] {x, y};
    }

    private double Fx (double P2_div_P1, double coslat1coslat1, double cosi) {
        return (P2_div_P1 * coslat1coslat1 - cosi) / StrictMath.sqrt(coslat1coslat1 - cosi * cosi);
    }

    private double lonpx (double lat, double sini) {
        return -StrictMath.asin(StrictMath.sin(lat) / sini);
    }

    private double lontx (double lonp, double cosi) {
        return StrictMath.atan(StrictMath.tan(lonp) * cosi);
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "i",
                "P2",
                "P1",
                "lon0",
                "lat0",
                "lat1",
                "lat2",
                "useConic"}));
    }

}