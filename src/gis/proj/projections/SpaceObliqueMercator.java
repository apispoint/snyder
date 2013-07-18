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
import gis.proj.SimpsonsRule;
import gis.proj.Space;

import java.lang.ref.SoftReference;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;


/**
 * All formulas adopted from USGS Professional Paper 1395.
 *
 * References:
 *
 * Map Projections - A Working Manual, USGS Professional Paper 1395
 * John P. Snyder
 * Pages 214-229
 *
 */

/*
 * Level III implementation definition
 *
 * References:
 *
 * Space Oblique Mercator Projection - Mathematical Development,
 * USGS Bulletin 1518
 * John P. Snyder
 * Pages 2-3
 *
 */
public final class SpaceObliqueMercator implements Space {

    // Default number of computed Fourier Coefficient sub-n's.
    // For sub-n's should be used for satellites with an orbital
    // eccentricity (e') of less than 0.05
    private static final int _N = 4;

    private SoftReference<Map<String, Datum>> fourierCacheMap = null;

    public String getName() {
        return "Space Oblique Mercator (Lvl III)";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double esq = ellip.getProperty("e^2");

        double _i     = datum.getProperty("i");
        double P2     = datum.getProperty("P2");
        double P1     = datum.getProperty("P1");
        double lon0   = datum.getProperty("lon0");
        double sat_ep = datum.getProperty("e'");
        double gamma  = datum.getProperty("gamma");

        // Undocumented on purpose
        double n = datum.getProperty("nsub");

        if(Double.isNaN(n) == true || (n = StrictMath.rint(n)) < _N)
            n = _N;

        Datum fourierCoeDatum = getFourierCoeffiencts(ellip, datum);

        double  J = fourierCoeDatum.getProperty("J");
        double  W = fourierCoeDatum.getProperty("W");
        double  Q = fourierCoeDatum.getProperty("Q");
        double  T = fourierCoeDatum.getProperty("T");
        double B1 = fourierCoeDatum.getProperty("B1");
        double B2 = fourierCoeDatum.getProperty("B2");
        double H1 = fourierCoeDatum.getProperty("H1");
        double S1 = fourierCoeDatum.getProperty("S1");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double P2_DIV_P1   = P2 / P1;
        double _1_S_esq    = 1.0 - esq;
        double _1_S_esq_sq = _1_S_esq * _1_S_esq;

        double cosi     = StrictMath.cos(_i);
        double sini     = StrictMath.sin(_i);
        double U        = esq * (cosi * cosi) / _1_S_esq;
        double _1_A_U   = 1.0 + U;

        double _1_S_esq_sini = _1_S_esq * sini;

        boolean _iNotNearZero = !(StrictMath.abs(_i) < NEAR_ZERO_RAD);

        double sat_ep_sq = sat_ep * sat_ep;
        double sqrt_sat_ep_sq = StrictMath.sqrt(1.0 - sat_ep_sq);

        double sqrt_1msatep_1psatep = StrictMath.sqrt((1 - sat_ep) / (1 + sat_ep));

        double V, lont, lonpp, phipp, coslonpp, sinlonpp, sinphipp;
        double L, Ep, xp, yp, sinphipp_sq;

        double sinlonpp_sq, oldlonpp, Lp, S;

        int k, j, l;

        String nStr;
        double[] AN  = new double[(int) n];
        double[] ApN = new double[(int) n];
        double[] CN  = new double[(int) n];
        double[] CpN = new double[(int) n];

        double[] tmpSums = new double[3];

        for(k = 1; k <= n; k++) {
            nStr = Integer.toString(k);
            AN[k-1]  = fourierCoeDatum.getProperty("A"  + nStr);
            ApN[k-1] = fourierCoeDatum.getProperty("A'" + nStr);
            CN[k-1]  = fourierCoeDatum.getProperty("C"  + nStr);
            CpN[k-1] = fourierCoeDatum.getProperty("C'" + nStr);
        }

        for(int i = 0; i < lon.length; ++i) {
            xp = (x[i] * H1 - y[i] * S1) / a;
            yp = (y[i] * H1 + x[i] * S1) / a;

            j = 0;
            lonpp = xp / B1;
            do {
                oldlonpp = lonpp;
            Ep = 2.0 * StrictMath.atan(
                    StrictMath.tan((lonpp - gamma) * 0.5) *
                    sqrt_1msatep_1psatep
                    );

            Lp = StrictMath.pow(1.0 - sat_ep * StrictMath.cos(Ep), 2) / sqrt_sat_ep_sq;

            sinlonpp    = StrictMath.sin(lonpp);
            sinlonpp_sq = sinlonpp * sinlonpp;

            S = P2_DIV_P1 * Lp * sini * StrictMath.cos(lonpp) * StrictMath.sqrt(
                    (1.0 + T * sinlonpp_sq) /
                    (
                    (1.0 + W * sinlonpp_sq) * (1.0 + Q * sinlonpp_sq)
                    )
                    );

            tmpSums[0] = tmpSums[1] = tmpSums[2] = 0;
            for(k = 1, l = 0; k <= n; k++, l++) {
                tmpSums[0] += (AN[l] + (S/J)*CN[l])*StrictMath.sin(k*lonpp);
                tmpSums[1] += (ApN[l] + (S/J)*CpN[l])*StrictMath.cos(k*lonpp);
                tmpSums[2] += (ApN[l] + (S/J)*CpN[l]);
            }

            lonpp = (xp + (S/J)*yp - tmpSums[0]+tmpSums[1]-tmpSums[2])
                    / (B1 + (S / J) * B2);

            } while (++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(lonpp - oldlonpp));

            tmpSums[0] = tmpSums[1] = tmpSums[2] = 0;
            for(k = 1, j = 0; k <= n; k++, j++) {
                tmpSums[0] += CN[j]  * StrictMath.sin(k * lonpp);
                tmpSums[1] += CpN[j] * StrictMath.cos(k * lonpp);
                tmpSums[2] += CpN[j];
            }

            phipp = (StrictMath.atan(StrictMath.pow(_e_,
                    StrictMath.sqrt(1+(S*S)/(J*J)) * 
                    (yp-B2*lonpp-tmpSums[0]+tmpSums[1]-tmpSums[2])
                    )) - PI_DIV_4) * 2.0;

            sinphipp = StrictMath.sin(phipp);
            coslonpp = StrictMath.cos(lonpp);
            sinlonpp = StrictMath.sin(lonpp);

            sinphipp_sq = sinphipp * sinphipp;

            Ep = 2.0 * StrictMath.atan(
                    StrictMath.tan((lonpp - gamma) * 0.5) *
                    sqrt_1msatep_1psatep
                    );

            L = Ep - sat_ep * StrictMath.sin(Ep);

            V = ((1.0 - sinphipp_sq / _1_S_esq) * cosi * sinlonpp - sini * sinphipp *
                    StrictMath.sqrt((1.0 + Q * sinlonpp * sinlonpp) * (1.0 - sinphipp_sq) - U * sinphipp_sq)) /
                    (1.0 - sinphipp_sq * _1_A_U);

            lont = StrictMath.atan2(V, coslonpp);

            lon[i] = normalizeLonRad(lont - P2_DIV_P1 * (L + gamma) + lon0);

            if(_iNotNearZero == true)
                lat[i] = StrictMath.atan((StrictMath.tan(lonpp) * StrictMath.cos(lont) -
                        cosi * StrictMath.sin(lont)) / _1_S_esq_sini);
            else
                lat[i] = StrictMath.asin(sinphipp / StrictMath.sqrt(_1_S_esq_sq + esq * sinphipp_sq));
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double esq = ellip.getProperty("e^2");

        double _i     = datum.getProperty("i");
        double P2     = datum.getProperty("P2");
        double P1     = datum.getProperty("P1");
        double lon0   = datum.getProperty("lon0");
        double sat_ep = datum.getProperty("e'");
        double gamma  = datum.getProperty("gamma");
        double N      = datum.getProperty("N");

        // Undocumented on purpose
        double n = datum.getProperty("nsub");

        if(Double.isNaN(n) == true || (n = StrictMath.rint(n)) < _N)
            n = _N;

        Datum fourierCoeDatum = getFourierCoeffiencts(ellip, datum);

        double  J = fourierCoeDatum.getProperty("J");
        double  W = fourierCoeDatum.getProperty("W");
        double  Q = fourierCoeDatum.getProperty("Q");
        double  T = fourierCoeDatum.getProperty("T");
        double B1 = fourierCoeDatum.getProperty("B1");
        double B2 = fourierCoeDatum.getProperty("B2");
        double H1 = fourierCoeDatum.getProperty("H1");
        double S1 = fourierCoeDatum.getProperty("S1");

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        int j, k;
        double xp, yp, lonpp, oldlonpp, Ep, L, lont, lonsubpp, sinlonsubpp;
        double phipp, sinlat, tanlat, S, Lp, sinlonpp, sinlonpp_sq;

        double P2_DIV_P1 = P2 / P1;
        double _1_M_esq = 1.0 - esq;

        double cosi = StrictMath.cos(_i);
        double sini = StrictMath.sin(_i);

        double sat_ep_sq = sat_ep * sat_ep;

        double sqrt_1msatep_1psatep = StrictMath.sqrt((1 - sat_ep) / (1 + sat_ep));

        double adjFactor, lontp, hypojs, logtan, cosklonpp, sinklonpp;
        String nStr;

        double[] AN  = new double[(int) n];
        double[] ApN = new double[(int) n];
        double[] CN  = new double[(int) n];
        double[] CpN = new double[(int) n];

        for(k = 1; k <= n; k++) {
            nStr = Integer.toString(k);
            AN[k-1]  = fourierCoeDatum.getProperty("A"  + nStr);
            ApN[k-1] = fourierCoeDatum.getProperty("A'" + nStr);
            CN[k-1]  = fourierCoeDatum.getProperty("C"  + nStr);
            CpN[k-1] = fourierCoeDatum.getProperty("C'" + nStr);
        }

        for(int i = 0; i < lon.length; ++i) {
            sinlat = StrictMath.sin(lat[i]);
            tanlat = StrictMath.tan(lat[i]);

            lonsubpp = POLE_RAD * (4 * N + 2.0 + (lat[i] <= 0 ? 1 : -1));
            sinlonsubpp = StrictMath.sin(lonsubpp);

            lontp = StrictMath.cos(lon[i] - lon0 + P2_DIV_P1 * lonsubpp);
            adjFactor = StrictMath.cos(lon[i] - lon0 + P2_DIV_P1 * lonsubpp) < 0 ? -1 : 1;
            if(StrictMath.abs(lontp) < NEAR_ZERO_RAD)
                adjFactor = 0;

            lonpp = lonsubpp;
            j = 0;
            do {
                oldlonpp = lonpp;

                Ep = 2.0 * StrictMath.atan(
                        StrictMath.tan((lonpp - gamma) * 0.5) *
                        sqrt_1msatep_1psatep
                        );

                L = Ep - sat_ep * StrictMath.sin(Ep);

                lont = lon[i] - lon0 + P2_DIV_P1 * (L + gamma);

                lonpp = StrictMath.atan(cosi * StrictMath.tan(lont) + _1_M_esq * sini * tanlat / StrictMath.cos(lont));

                lonpp += lonsubpp - POLE_RAD * sinlonsubpp * adjFactor;

            } while (++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(lonpp - oldlonpp));
//            lonpp = oldlonpp;

            sinlonpp    = StrictMath.sin(lonpp);
            sinlonpp_sq = sinlonpp * sinlonpp;

            phipp = StrictMath.asin(
                    (_1_M_esq * cosi * sinlat - sini * StrictMath.cos(lat[i]) * StrictMath.sin(lont))
                    /
                    StrictMath.sqrt(1.0 - esq * sinlat * sinlat));

            Lp = StrictMath.pow(1.0 - sat_ep * StrictMath.cos(Ep), 2) / StrictMath.sqrt(1.0 - sat_ep_sq);

            S = P2_DIV_P1 * Lp * sini * StrictMath.cos(lonpp) * StrictMath.sqrt(
                    (1.0 + T * sinlonpp_sq) /
                    (
                    (1.0 + W * sinlonpp_sq) * (1.0 + Q * sinlonpp_sq)
                    )
                    );

            hypojs = StrictMath.hypot(J, S);
            logtan = StrictMath.log(StrictMath.tan(PI_DIV_4 + phipp * 0.5));

            xp = B1 * lonpp - (S / hypojs) * logtan;
            yp = B2 * lonpp + (J / hypojs) * logtan;

            for(k = 1, j = 0; k <= n; k++, j++) {
                cosklonpp = StrictMath.cos(k * lonpp);
                sinklonpp = StrictMath.sin(k * lonpp);

                xp += AN[j]  * sinklonpp;
                xp -= ApN[j] * cosklonpp;
                xp += ApN[j];

                yp += CN[j]  * sinklonpp;
                yp -= CpN[j] * cosklonpp;
                yp += CpN[j];
            }

            x[i] = a * (xp * H1 + yp * S1);
            y[i] = a * (yp * H1 - xp * S1);
        }

        return new double[][] {x, y};
    }

    private synchronized Datum getFourierCoeffiencts(Ellipsoid ellip, Datum datum) {
        if(fourierCacheMap == null || fourierCacheMap.get() == null)
            fourierCacheMap = new SoftReference<Map<String,Datum>>(new HashMap<String, Datum>());

        Map<String, Datum> fourierMap = fourierCacheMap.get();
        Datum fourierCoeDatum = fourierMap.get(ellip.getId());

        if(fourierCoeDatum == null) {
            fourierCoeDatum = new Datum();
            fourierMap.put(ellip.getId(), fourierCoeDatum);

            double esq = ellip.getProperty("e^2");

            double sat_ep = datum.getProperty("e'");

            // Undocumented on purpose
            double n = datum.getProperty("nsub");

            if(Double.isNaN(n) == true || (n = StrictMath.rint(n)) < _N)
                n = _N;

            double gamma = datum.getProperty("gamma");
            double i     = datum.getProperty("i");
            double P2    = datum.getProperty("P2");
            double P1    = datum.getProperty("P1");

            double _1_MIN_ESQ  = 1.0 - esq;
            double _1MINESQ_SQ = _1_MIN_ESQ * _1_MIN_ESQ;

            double cosi     = StrictMath.cos(i);
            double sini     = StrictMath.sin(i);
            double sinisini = sini * sini;

            double J = StrictMath.pow(_1_MIN_ESQ, 3);
            double W = (StrictMath.pow(1.0 - esq * cosi * cosi, 2) / _1MINESQ_SQ) - 1.0;
            double Q = esq * sinisini / _1_MIN_ESQ;
            double T = esq * sinisini * (2.0 - esq) / _1MINESQ_SQ;

            fourierCoeDatum.setUserOverrideProperty("J", J);
            fourierCoeDatum.setUserOverrideProperty("W", W);
            fourierCoeDatum.setUserOverrideProperty("Q", Q);
            fourierCoeDatum.setUserOverrideProperty("T", T);

            final SimpsonsRule sr = SimpsonsRule.getInstance();

            final SimpsonsRule.Function somCoes = new SimpsonsRule.Function() {
                public double[] f(double[] parameters, double x) {
                    double ep        = parameters[0];
                    double gamma     = parameters[1];
                    double P2_DIV_P1 = parameters[2];
                    double J         = parameters[3];
                    double W         = parameters[4];
                    double Q         = parameters[5];
                    double T         = parameters[6];
                    double cosi      = parameters[7];
                    double sini      = parameters[8];
                    int n            = (int) parameters[9];

                    double lonpp = x;
                    double epep = ep * ep;

                    double Ep = 2.0 * StrictMath.atan(
                            StrictMath.tan((lonpp - gamma) * 0.5) *
                            parameters[10]
                            );

                    double Lp = StrictMath.pow(1.0 - ep * StrictMath.cos(Ep), 2) / StrictMath.sqrt(1.0 - epep);

                    double sinlonpp    = StrictMath.sin(lonpp);
                    double sinlonpp_sq = sinlonpp * sinlonpp;

                    double S = P2_DIV_P1 * Lp * sini * StrictMath.cos(lonpp) * StrictMath.sqrt(
                            (1.0 + T * sinlonpp_sq) /
                            (
                            (1.0 + W * sinlonpp_sq) * (1.0 + Q * sinlonpp_sq)
                            )
                            );

                    double H =
                            StrictMath.sqrt((1 + Q * sinlonpp_sq) / (1 + W * sinlonpp_sq)) *
                            (((1 + W * sinlonpp_sq) / StrictMath.pow(1+Q*sinlonpp_sq, 2)) - P2_DIV_P1 * Lp * cosi);

                    double den     = StrictMath.hypot(J, S);
                    double hj_M_ss = ((H * J) - (S * S)) / den;
                    double s_T_hPj = S * (H + J)         / den;

                    // Array size is realized as follows:
                    //
                    //        2 - "static" return values
                    //
                    //    4 * n - "dynamic" return value
                    //
                    // Where:
                    //    4 - represents the distinct values
                    //        of (An, A'n, Cn, C'n)
                    //
                    //    n - represents the number of sets (sub-n's)
                    //        to compute with n being no less than 4
                    //        resulting in a minimum of 16 terms.
                    double[] vector = new double[2 + (4 * n)];

                    vector[0] = hj_M_ss;
                    vector[1] = s_T_hPj;

                    double cosilonpp, sinilonpp;
                    for(int i = 1, j = 2; i <= n; i++, j += 4) {
                        cosilonpp = StrictMath.cos(i * lonpp);
                        sinilonpp = StrictMath.sin(i * lonpp);

                        vector[j]   = hj_M_ss * cosilonpp;
                        vector[j+1] = hj_M_ss * sinilonpp;
                        vector[j+2] = s_T_hPj * cosilonpp;
                        vector[j+3] = s_T_hPj * sinilonpp;
                    }

                    return vector;
                }
            };

            double[] vector = sr.approximate(somCoes,
                    new double[] {
                    sat_ep,
                    gamma,
                    P2 / P1,
                    J, W, Q, T, cosi, sini, n,
                    StrictMath.sqrt((1 - sat_ep) / (1 + sat_ep))
            },
            0, 2.0 * PI, 40);

            double B1 = _1_DIV_PI * 0.5 * vector[0];
            double B2 = _1_DIV_PI * 0.5 * vector[1];
            double H1 = B1 / StrictMath.hypot(B1, B2);
            double S1 = B2 / StrictMath.hypot(B1, B2);

            fourierCoeDatum.setUserOverrideProperty("B1", B1);
            fourierCoeDatum.setUserOverrideProperty("B2", B2);
            fourierCoeDatum.setUserOverrideProperty("H1", H1);
            fourierCoeDatum.setUserOverrideProperty("S1", S1);

            double factor;
            String kStr;

            for(int k = 1, j = 2; k <= n; k++, j += 4) {
                kStr   = Integer.toString(k);
                factor = 1.0 / (PI * k);

                fourierCoeDatum.setUserOverrideProperty("A"  + kStr, factor * vector[j]);
                fourierCoeDatum.setUserOverrideProperty("A'" + kStr, factor * vector[j+1]);
                fourierCoeDatum.setUserOverrideProperty("C"  + kStr, factor * vector[j+2]);
                fourierCoeDatum.setUserOverrideProperty("C'" + kStr, factor * vector[j+3]);
            }
        }

        return fourierCoeDatum;
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "i",
                "P2",
                "P1",
                "lon0",
//                "a'",
                "e'",
                "gamma",
                "N"}));
    }

}