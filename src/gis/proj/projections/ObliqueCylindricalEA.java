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
import gis.proj.SimpsonsRule;
import gis.proj.SnyderMath;

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
 * Pages 76-85
 *
 */

public final class ObliqueCylindricalEA implements Cylindrical {

    private SoftReference<Map<String, Datum>> fourierCacheMap = null;

    public String getName() {
        return "Oblique Cylindrical Equal-Area";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double e   = ellip.getProperty("e");
        double esq = ellip.getProperty("e^2");

        double h0   = datum.getProperty("h0");

        // Used for Central Line
        double lon1 = datum.getProperty("lon1");
        double lat1 = datum.getProperty("lat1");
        double lon2 = datum.getProperty("lon2");
        double lat2 = datum.getProperty("lat2");

        // Used for Central Point
        double lonz  = datum.getProperty("lonz");
        double latz  = datum.getProperty("latz");
        double gamma = datum.getProperty("gamma");

        // Given or calculated from Central (Point | Line)
        double lonp = datum.getProperty("lonp");
        double latp = datum.getProperty("latp");

        boolean useCentralLine  = datum.getProperty("useCentralLine")  == 1.0;
        boolean useCentralPoint = datum.getProperty("useCentralPoint") == 1.0;

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        if(useCentralLine == true && ellip.isSphere() == false) {
            // (9-1) (9-2)
            // Two points
            double sinlat1 = StrictMath.sin(lat1);
            double sinlat2 = StrictMath.sin(lat2);

            double one_M_esq = 1.0 - esq;
            double one_over_2e = 1.0 / (2.0 * e);

            double sin_hlf_pi = StrictMath.sin(POLE_RAD);

            double qp = one_M_esq * (
                    sin_hlf_pi / (1.0 - esq*sin_hlf_pi*sin_hlf_pi) - one_over_2e * StrictMath.log((1.0 - e*sin_hlf_pi) / (1.0 + e*sin_hlf_pi))
                    );

            double q = one_M_esq * (
                    sinlat1 / (1.0 - esq*sinlat1*sinlat1) - one_over_2e * StrictMath.log((1.0 - e*sinlat1) / (1.0 + e*sinlat1))
                    );
            double beta1 = StrictMath.asin(q / qp);

            q = one_M_esq * (
                    sinlat2 / (1.0 - esq*sinlat2*sinlat2) - one_over_2e * StrictMath.log((1.0 - e*sinlat2) / (1.0 + e*sinlat2))
                    );
            double beta2 = StrictMath.asin(q / qp);

            double cosbeta1 = StrictMath.cos(beta1);
            double sinbeta1 = StrictMath.sin(beta1);

            double cosbeta2 = StrictMath.cos(beta2);
            double sinbeta2 = StrictMath.sin(beta2);

            double dom = sinbeta1 * cosbeta2 * StrictMath.sin(lon2) - cosbeta1 * sinbeta2 * StrictMath.sin(lon1);

            lonp = StrictMath.atan2(
                    cosbeta1 * sinbeta2 * StrictMath.cos(lon1) - sinbeta1 * cosbeta2 * StrictMath.cos(lon2),
                    dom
                    );

            if(dom < 0)
                lonp += SnyderMath.PI;

            latp = StrictMath.atan(-StrictMath.cos(lonp - lon1) / StrictMath.tan(beta1));

            // (3-18)
            latp = latp + (esq / 3.0 + (esq *esq) * (31.0/180.0) + (esq*esq*esq) * (517.0/5040.0)) * StrictMath.sin(2.0 * latp) +
                    ((esq *esq) * (23.0/360.0) + (esq*esq*esq) * (251.0/3780.0)) * StrictMath.sin(4.0 * latp) +
                    ((esq*esq*esq) * (761.0/45360.0)) * StrictMath.sin(6.0 * latp);
        }
        else if(useCentralLine == true && ellip.isSphere() == true) {
            // (9-1) (9-2)
            // Two points
            double coslat1 = StrictMath.cos(lat1);
            double sinlat1 = StrictMath.sin(lat1);
            double coslat2 = StrictMath.cos(lat2);
            double sinlat2 = StrictMath.sin(lat2);
            double dom = sinlat1 * coslat2 * StrictMath.sin(lon2) - coslat1 * sinlat2 * StrictMath.sin(lon1);

            lonp = StrictMath.atan2(
                    coslat1 * sinlat2 * StrictMath.cos(lon1) - sinlat1 * coslat2 * StrictMath.cos(lon2),
                    dom
                    );

            if(dom < 0)
                lonp += SnyderMath.PI;

            latp = StrictMath.atan(-StrictMath.cos(lonp - lon1) / StrictMath.tan(lat1));
        }
        else if(useCentralPoint == true) {
            // (9-7) (9-8)
            // Center point and gamma (rotation)
            double singamma = StrictMath.sin(gamma);

            lonp = normalizeLonRad(StrictMath.atan2(-StrictMath.cos(gamma), -StrictMath.sin(latz) * singamma) + lonz);
            latp = StrictMath.asin(StrictMath.cos(latz) * singamma);
        }

        if(ellip.isSphere()) {
            double lon0 = lonp + PI_DIV_2;
            double h0_over_a = h0 / a;
            double ah0 = a * h0;

            double coslatp = StrictMath.cos(latp);
            double sinlatp = StrictMath.sin(latp);

            double xfactor, yfactor, yfactorsq, sqrt_1_M_yfactorsq;

            for(int i = 0; i < lon.length; ++i) {
                xfactor = x[i] / ah0;

                yfactor = (y[i] * h0_over_a);
                yfactorsq = yfactor * yfactor;
                sqrt_1_M_yfactorsq = StrictMath.sqrt(1.0 - yfactorsq);

                lon[i] = normalizeLonRad(lon0 + StrictMath.atan2(
                        sqrt_1_M_yfactorsq * sinlatp * StrictMath.sin(xfactor) - yfactor * coslatp,
                        sqrt_1_M_yfactorsq * StrictMath.cos(xfactor)
                        ));

                lat[i] = StrictMath.asin(yfactor * sinlatp + sqrt_1_M_yfactorsq * coslatp * StrictMath.sin(xfactor));
            }
        }
        else {
            double one_M_esq = 1.0 - esq;
            double one_over_2e = 1.0 / (2.0 * e);

            double sin_hlf_pi = StrictMath.sin(POLE_RAD);

            // (3-12)
            double qp = one_M_esq * (
                    sin_hlf_pi / (1.0 - esq*sin_hlf_pi*sin_hlf_pi) - one_over_2e * StrictMath.log((1.0 - e*sin_hlf_pi) / (1.0 + e*sin_hlf_pi))
                    );

            double sinlatp = StrictMath.sin(latp);
            double esinlatp = e * sinlatp;
            double esqsinlatpsinlatp = esq * sinlatp * sinlatp;

            double q = one_M_esq * (
                    sinlatp / (1.0 - esqsinlatpsinlatp) - one_over_2e * StrictMath.log((1.0 - esinlatp) / (1.0 + esinlatp))
                    );
            double betap = StrictMath.asin(q / qp);

            Datum fourierCoeDatum = getFourierCoeffiencts(ellip);

            double b    = fourierCoeDatum.getProperty("fc_b");
            double a2   = fourierCoeDatum.getProperty("fc_a2");
            double a4   = fourierCoeDatum.getProperty("fc_a4");
            double a6   = fourierCoeDatum.getProperty("fc_a6");
            double b2   = fourierCoeDatum.getProperty("fc_b2");
            double b4   = fourierCoeDatum.getProperty("fc_b4");
            double ap22 = fourierCoeDatum.getProperty("fc_a'22");
            double ap24 = fourierCoeDatum.getProperty("fc_a'24");
            double ap26 = fourierCoeDatum.getProperty("fc_a'26");
            double ap42 = fourierCoeDatum.getProperty("fc_a'42");
            double ap44 = fourierCoeDatum.getProperty("fc_a'44");
            double ap46 = fourierCoeDatum.getProperty("fc_a'46");

            double B = b + a2 * StrictMath.cos(2.0 * latp) + a4 * StrictMath.cos(4.0 * latp) + a6 * StrictMath.cos(6.0 * latp);
            double A2 = b2 + ap22 * StrictMath.cos(2.0 * latp) + ap24 * StrictMath.cos(4.0 * latp) + ap26 * StrictMath.cos(6.0 * latp);
            double A4 = b4 + ap42 * StrictMath.cos(2.0 * latp) + ap44 * StrictMath.cos(4.0 * latp) + ap46 * StrictMath.cos(6.0 * latp);

            double A6 =  0.0;

            double A2_2 = 2.0 * A2;
            double A4_4 = 4.0 * A4;
            double A6_6 = 6.0 * A6;

            double ah0  = a * h0;
            double ah0B = ah0 * B;

            int j;
            double lonPrime, lonPrime_old, F, betaPrime, beta;

            for(int i = 0; i < lon.length; ++i) {
                lonPrime = x[i] / ah0B;
                j = 0;
                do {
                    lonPrime_old = lonPrime;
                    lonPrime = (x[i] / ah0 - A2 * StrictMath.sin(2.0 * lonPrime) - A4 * StrictMath.sin(4.0 * lonPrime) - A6 * StrictMath.sin(6.0 * lonPrime)) / B;
                } while(++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(lonPrime - lonPrime_old));

                F = B + A2_2 * StrictMath.cos(2.0 * lonPrime) + A4_4 * StrictMath.cos(4.0 * lonPrime) + A6_6 * StrictMath.cos(6.0 * lonPrime);

                betaPrime = StrictMath.asin(2.0 * F * h0 * y[i] / (a * qp));
                beta = StrictMath.asin(StrictMath.sin(betap) * StrictMath.sin(betaPrime) + StrictMath.cos(betap) * StrictMath.cos(betaPrime) * StrictMath.sin(lonPrime));

                lon[i] = normalizeLonRad(
                        lonp + StrictMath.atan2(StrictMath.cos(betaPrime) * StrictMath.cos(lonPrime), StrictMath.cos(betap) * StrictMath.sin(betaPrime) - StrictMath.sin(betap) * StrictMath.cos(betaPrime) * StrictMath.sin(lonPrime))
                        );

                // (3-18)
                lat[i] = beta + (esq / 3.0 + (esq *esq) * (31.0/180.0) + (esq*esq*esq) * (517.0/5040.0)) * StrictMath.sin(2.0 * beta) +
                        ((esq *esq) * (23.0/360.0) + (esq*esq*esq) * (251.0/3780.0)) * StrictMath.sin(4.0 * beta) +
                        ((esq*esq*esq) * (761.0/45360.0)) * StrictMath.sin(6.0 * beta);
            }
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a   = ellip.getProperty("a");
        double e   = ellip.getProperty("e");
        double esq = ellip.getProperty("e^2");

        double h0   = datum.getProperty("h0");

        // Used for Central Line
        double lon1 = datum.getProperty("lon1");
        double lat1 = datum.getProperty("lat1");
        double lon2 = datum.getProperty("lon2");
        double lat2 = datum.getProperty("lat2");

        // Used for Central Point
        double lonz  = datum.getProperty("lonz");
        double latz  = datum.getProperty("latz");
        double gamma = datum.getProperty("gamma");

        // Given or calculated from Central (Point | Line)
        double lonp = datum.getProperty("lonp");
        double latp = datum.getProperty("latp");

        boolean useCentralLine  = datum.getProperty("useCentralLine")  == 1.0;
        boolean useCentralPoint = datum.getProperty("useCentralPoint") == 1.0;

        double x[] = new double[lon.length];
        double y[] = new double[lat.length];

        if(useCentralLine == true && ellip.isSphere() == false) {
            // (9-1) (9-2)
            // Two points
            double sinlat1 = StrictMath.sin(lat1);
            double sinlat2 = StrictMath.sin(lat2);

            double one_M_esq = 1.0 - esq;
            double one_over_2e = 1.0 / (2.0 * e);

            double sin_hlf_pi = StrictMath.sin(POLE_RAD);

            double qp = one_M_esq * (
                    sin_hlf_pi / (1.0 - esq*sin_hlf_pi*sin_hlf_pi) - one_over_2e * StrictMath.log((1.0 - e*sin_hlf_pi) / (1.0 + e*sin_hlf_pi))
                    );

            double q = one_M_esq * (
                    sinlat1 / (1.0 - esq*sinlat1*sinlat1) - one_over_2e * StrictMath.log((1.0 - e*sinlat1) / (1.0 + e*sinlat1))
                    );
            double beta1 = StrictMath.asin(q / qp);

            q = one_M_esq * (
                    sinlat2 / (1.0 - esq*sinlat2*sinlat2) - one_over_2e * StrictMath.log((1.0 - e*sinlat2) / (1.0 + e*sinlat2))
                    );
            double beta2 = StrictMath.asin(q / qp);

            double cosbeta1 = StrictMath.cos(beta1);
            double sinbeta1 = StrictMath.sin(beta1);

            double cosbeta2 = StrictMath.cos(beta2);
            double sinbeta2 = StrictMath.sin(beta2);

            double dom = sinbeta1 * cosbeta2 * StrictMath.sin(lon2) - cosbeta1 * sinbeta2 * StrictMath.sin(lon1);

            lonp = StrictMath.atan2(
                    cosbeta1 * sinbeta2 * StrictMath.cos(lon1) - sinbeta1 * cosbeta2 * StrictMath.cos(lon2),
                    dom
                    );

            if(dom < 0)
                lonp += SnyderMath.PI;

            latp = StrictMath.atan(-StrictMath.cos(lonp - lon1) / StrictMath.tan(beta1));

            // (3-18)
            latp = latp + (esq / 3.0 + (esq *esq) * (31.0/180.0) + (esq*esq*esq) * (517.0/5040.0)) * StrictMath.sin(2.0 * latp) +
                    ((esq *esq) * (23.0/360.0) + (esq*esq*esq) * (251.0/3780.0)) * StrictMath.sin(4.0 * latp) +
                    ((esq*esq*esq) * (761.0/45360.0)) * StrictMath.sin(6.0 * latp);
        }
        else if(useCentralLine == true && ellip.isSphere() == true) {
            // (9-1) (9-2)
            // Two points
            double coslat1 = StrictMath.cos(lat1);
            double sinlat1 = StrictMath.sin(lat1);
            double coslat2 = StrictMath.cos(lat2);
            double sinlat2 = StrictMath.sin(lat2);
            double dom = sinlat1 * coslat2 * StrictMath.sin(lon2) - coslat1 * sinlat2 * StrictMath.sin(lon1);

            lonp = StrictMath.atan2(
                    coslat1 * sinlat2 * StrictMath.cos(lon1) - sinlat1 * coslat2 * StrictMath.cos(lon2),
                    dom
                    );

            if(dom < 0)
                lonp += SnyderMath.PI;

            latp = StrictMath.atan(-StrictMath.cos(lonp - lon1) / StrictMath.tan(lat1));
        }
        else if(useCentralPoint == true) {
            // (9-7) (9-8)
            // Center point and gamma (rotation)
            double singamma = StrictMath.sin(gamma);

            lonp = normalizeLonRad(StrictMath.atan2(-StrictMath.cos(gamma), -StrictMath.sin(latz) * singamma) + lonz);
            latp = StrictMath.asin(StrictMath.cos(latz) * singamma);
        }

        if(ellip.isSphere()) {
            double lon0 = lonp + PI_DIV_2;
            double a_over_h0 = a / h0;
            double ah0 = a * h0;

            double coslatp = StrictMath.cos(latp);
            double sinlatp = StrictMath.sin(latp);

            double lon_M_lon0, sinlon_M_lon0, dom;

            for(int i = 0; i < lon.length; ++i) {
                lon_M_lon0 = normalizeLonRad(lon[i] - lon0);
                sinlon_M_lon0 = StrictMath.sin(lon_M_lon0);

                dom = StrictMath.cos(lon_M_lon0);

                x[i] = ah0 * (
                        StrictMath.atan(
                                (StrictMath.tan(lat[i]) * coslatp + sinlatp * sinlon_M_lon0) /
                                dom)
                                + (dom < 0 ? SnyderMath.PI : 0.0));

                y[i] = a_over_h0 * (
                        sinlatp * StrictMath.sin(lat[i]) -
                        coslatp * StrictMath.cos(lat[i]) * sinlon_M_lon0);
            }
        }
        else {
            double one_M_esq = 1.0 - esq;
            double one_over_2e = 1.0 / (2.0 * e);

            double sin_hlf_pi = StrictMath.sin(POLE_RAD);

            double qp = one_M_esq * (
                    sin_hlf_pi / (1.0 - esq*sin_hlf_pi*sin_hlf_pi) - one_over_2e * StrictMath.log((1.0 - e*sin_hlf_pi) / (1.0 + e*sin_hlf_pi))
                    );

            double sinlatp = StrictMath.sin(latp);
            double esinlatp = e * sinlatp;
            double esqsinlatpsinlatp = esq * sinlatp * sinlatp;

            double q = one_M_esq * (
                    sinlatp / (1.0 - esqsinlatpsinlatp) - one_over_2e * StrictMath.log((1.0 - esinlatp) / (1.0 + esinlatp))
                    );
            double betap = StrictMath.asin(q / qp);
            double cosbetap = StrictMath.cos(betap);
            double sinbetap = StrictMath.sin(betap);

            double ah0 = a * h0;
            double hlf_aqp = a * qp * 0.5;

            double F, sinlat, esinlat, esqsinlatsinlat, beta, lonPrime, lon_M_lonp;

            double lonPrime2, lonPrime4, lonPrime6;

            Datum fourierCoeDatum = getFourierCoeffiencts(ellip);

            double b    = fourierCoeDatum.getProperty("fc_b");
            double a2   = fourierCoeDatum.getProperty("fc_a2");
            double a4   = fourierCoeDatum.getProperty("fc_a4");
            double a6   = fourierCoeDatum.getProperty("fc_a6");
            double b2   = fourierCoeDatum.getProperty("fc_b2");
            double b4   = fourierCoeDatum.getProperty("fc_b4");
            double ap22 = fourierCoeDatum.getProperty("fc_a'22");
            double ap24 = fourierCoeDatum.getProperty("fc_a'24");
            double ap26 = fourierCoeDatum.getProperty("fc_a'26");
            double ap42 = fourierCoeDatum.getProperty("fc_a'42");
            double ap44 = fourierCoeDatum.getProperty("fc_a'44");
            double ap46 = fourierCoeDatum.getProperty("fc_a'46");

            double B = b + a2 * StrictMath.cos(2.0 * latp) + a4 * StrictMath.cos(4.0 * latp) + a6 * StrictMath.cos(6.0 * latp);
            double A2 = b2 + ap22 * StrictMath.cos(2.0 * latp) + ap24 * StrictMath.cos(4.0 * latp) + ap26 * StrictMath.cos(6.0 * latp);
            double A4 = b4 + ap42 * StrictMath.cos(2.0 * latp) + ap44 * StrictMath.cos(4.0 * latp) + ap46 * StrictMath.cos(6.0 * latp);

            double A6 =  0.0;

            double A2_2 = 2.0 * A2;
            double A4_4 = 4.0 * A4;
            double A6_6 = 6.0 * A6;

            double coslon_M_lonp, cosbeta, sinbeta, dom;

            for(int i = 0; i < lon.length; ++i) {
                sinlat = StrictMath.sin(lat[i]);
                esinlat = e * sinlat;
                esqsinlatsinlat = esq * sinlat * sinlat;

                //eq3_12
                q = one_M_esq * (
                        sinlat / (1.0 - esqsinlatsinlat) - one_over_2e * StrictMath.log((1.0 - esinlat) / (1.0 + esinlat))
                        );
                beta = StrictMath.asin(q / qp);
                cosbeta = StrictMath.cos(beta);
                sinbeta = StrictMath.sin(beta);

                lon_M_lonp = normalizeLonRad(lon[i] - lonp);
                coslon_M_lonp = StrictMath.cos(lon_M_lonp);

                dom = cosbeta * StrictMath.sin(lon_M_lonp);

                lonPrime = StrictMath.atan(
                        (cosbetap * sinbeta - sinbetap * cosbeta * coslon_M_lonp)/
                        dom
                        );

                if(dom < 0)
                    lonPrime += SnyderMath.PI;

                lonPrime2 = 2.0 * lonPrime;
                lonPrime4 = 4.0 * lonPrime;
                lonPrime6 = 6.0 * lonPrime;

                x[i] = ah0 * (B * lonPrime + A2 * StrictMath.sin(lonPrime2) + A4 * StrictMath.sin(lonPrime4) + A6 * StrictMath.sin(lonPrime6));

                F = B + A2_2 * StrictMath.cos(lonPrime2) + A4_4 * StrictMath.cos(lonPrime4) + A6_6 * StrictMath.cos(lonPrime6);

                y[i] = hlf_aqp * (sinbetap * sinbeta + cosbetap * cosbeta * coslon_M_lonp) / (h0 * F);
            }
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "h0",
                "lon1",
                "lat1",
                "lon2",
                "lat2",
                "lonz",
                "latz",
                "gamma",
                "lonp",
                "latp",
                "useCentralLine",
                "useCentralPoint"}));
    }

    private synchronized Datum getFourierCoeffiencts(Ellipsoid ellip) {
        if(fourierCacheMap == null || fourierCacheMap.get() == null)
            fourierCacheMap = new SoftReference<Map<String,Datum>>(new HashMap<String, Datum>());

        Map<String, Datum> fourierMap = fourierCacheMap.get();
        Datum fourierCoeDatum = fourierMap.get(ellip.getId());

        if(fourierCoeDatum == null) {
            fourierCoeDatum = new Datum();
            fourierMap.put(ellip.getId(), fourierCoeDatum);
            double e   = ellip.getProperty("e");
            double esq = ellip.getProperty("e^2");

            double one_M_esq   = 1.0 - esq;
            double one_over_2e = 1.0 / (2.0 * e);

            double sin_hlf_pi = StrictMath.sin(POLE_RAD);

            double qp = one_M_esq * (
                    sin_hlf_pi / (1.0 - esq * sin_hlf_pi * sin_hlf_pi) -
                    one_over_2e * StrictMath.log(
                            (1.0 - e * sin_hlf_pi) / (1.0 + e * sin_hlf_pi))
                    );

            final double aLim      = 0;
            final double bLim      = PI_DIV_2;
            final int    intervals = 10;

            final SimpsonsRule sr = SimpsonsRule.getInstance();

            final SimpsonsRule.Function F = new SimpsonsRule.Function() {
                public double[] f(double[] parameters, double x) {
                    double cosbetap         = parameters[0];
                    double cosbetapcosbetap = parameters[1];
                    double sinbetapsinbetap = parameters[2];
                    double esq              = parameters[3];
                    double qp               = parameters[4];
                    double phip             = parameters[5];
                    double qpqp             = qp * qp;
                    double lonp             = x;

                    double v;

                    // Short circuit if-return
                    if(
                            phip < NEAR_ZERO_RAD &&
                            StrictMath.abs(POLE_RAD - lonp) < NEAR_ZERO_RAD) {
                        v = StrictMath.sqrt(qp * 0.5);
                        return new double[] {
                                v,
                                v * StrictMath.cos(2.0 * lonp),
                                v * StrictMath.cos(4.0 * lonp)
                        };
                    }

                    double coslonp = StrictMath.cos(lonp);
                    double coslonpcoslonp = coslonp * coslonp;

                    double betac = StrictMath.asin(cosbetap * StrictMath.sin(lonp));

                    // (3-18)
                    double phic = betac + (esq / 3.0 + (esq *esq) * (31.0/180.0) + (esq*esq*esq) * (517.0/5040.0)) * StrictMath.sin(2.0 * betac) +
                            ((esq *esq) * (23.0/360.0) + (esq*esq*esq) * (251.0/3780.0)) * StrictMath.sin(4.0 * betac) +
                            ((esq*esq*esq) * (761.0/45360.0)) * StrictMath.sin(6.0 * betac);

                    double cosphic = StrictMath.cos(phic);
                    double cosphiccosphic = cosphic * cosphic;

                    double sinphic = StrictMath.sin(phic);
                    double one_M_esqsinphicsinphic = 1.0 - esq * sinphic * sinphic;

                    v = StrictMath.sqrt(
                            sinbetapsinbetap * cosphiccosphic / (one_M_esqsinphicsinphic * StrictMath.pow(StrictMath.cos(betac), 4)) +
                            one_M_esqsinphicsinphic * qpqp * cosbetapcosbetap * coslonpcoslonp / (4.0 * cosphiccosphic)
                            );

                    return new double[] {
                            v,
                            v * StrictMath.cos(2.0 * lonp),
                            v * StrictMath.cos(4.0 * lonp)
                    };
                }
            };

            final SimpsonsRule.Function B = new SimpsonsRule.Function() {
                public double[] f(double[] parameters, double x) {
                    double e           = parameters[0];
                    double esq         = parameters[1];
                    double one_M_esq   = parameters[2];
                    double one_over_2e = parameters[3];
                    double qp          = parameters[4];
                    double phip        = x;

                    double sinphip, esinphip, q, betap, cosbetap, sinbetap;
                    double cosbetapcosbetap, sinbetapsinbetap;

                    // (3-12)
                    sinphip  = StrictMath.sin(phip);
                    esinphip = e * sinphip;

                    q = one_M_esq * (
                            sinphip / (1.0 - esq * sinphip * sinphip) - one_over_2e * StrictMath.log((1.0 - esinphip) / (1.0 + esinphip))
                            );

                    // (3-11)
                    betap = StrictMath.asin(q / qp);

                    cosbetap = StrictMath.cos(betap);
                    sinbetap = StrictMath.sin(betap);

                    cosbetapcosbetap = cosbetap * cosbetap;
                    sinbetapsinbetap = sinbetap * sinbetap;

                    double[] F_vect = sr.approximate(F, new double[] {
                            cosbetap,
                            cosbetapcosbetap,
                            sinbetapsinbetap,
                            esq,
                            qp,
                            phip
                    }, aLim, bLim, intervals);

                    double B = (2.0 / PI) * F_vect[0];

                    double _4_DIV_PI2_FV1 = (4.0 / (PI * 2.0)) * F_vect[1];
                    double _4_DIV_PI4_FV2 = (4.0 / (PI * 4.0)) * F_vect[2];

                    double cos2phip = StrictMath.cos(2.0 * phip);
                    double cos4phip = StrictMath.cos(4.0 * phip);
                    double cos6phip = StrictMath.cos(6.0 * phip);

                    return new double[] {
                            B,
                            B * cos2phip,
                            B * cos4phip,
                            B * cos6phip,

                            _4_DIV_PI2_FV1,
                            _4_DIV_PI4_FV2,

                            _4_DIV_PI2_FV1 * cos2phip,
                            _4_DIV_PI2_FV1 * cos4phip,
                            _4_DIV_PI2_FV1 * cos6phip,

                            _4_DIV_PI4_FV2 * cos2phip,
                            _4_DIV_PI4_FV2 * cos4phip,
                            _4_DIV_PI4_FV2 * cos6phip
                    };
                }
            };

            double[] B_vect = sr.approximate(B, new double[] {
                    e,
                    esq,
                    one_M_esq,
                    one_over_2e,
                    qp
            }, aLim, bLim, intervals);

            double _2_DIV_PI = 2.0 / PI;
            double _4_DIV_PI = 4.0 / PI;

            fourierCoeDatum.setUserOverrideProperty("fc_b",    _2_DIV_PI * B_vect[ 0]);
            fourierCoeDatum.setUserOverrideProperty("fc_a2",   _4_DIV_PI * B_vect[ 1]);
            fourierCoeDatum.setUserOverrideProperty("fc_a4",   _4_DIV_PI * B_vect[ 2]);
            fourierCoeDatum.setUserOverrideProperty("fc_a6",   _4_DIV_PI * B_vect[ 3]);
            fourierCoeDatum.setUserOverrideProperty("fc_b2",   _2_DIV_PI * B_vect[ 4]);
            fourierCoeDatum.setUserOverrideProperty("fc_b4",   _2_DIV_PI * B_vect[ 5]);
            fourierCoeDatum.setUserOverrideProperty("fc_a'22", _4_DIV_PI * B_vect[ 6]);
            fourierCoeDatum.setUserOverrideProperty("fc_a'24", _4_DIV_PI * B_vect[ 7]);
            fourierCoeDatum.setUserOverrideProperty("fc_a'26", _4_DIV_PI * B_vect[ 8]);
            fourierCoeDatum.setUserOverrideProperty("fc_a'42", _4_DIV_PI * B_vect[ 9]);
            fourierCoeDatum.setUserOverrideProperty("fc_a'44", _4_DIV_PI * B_vect[10]);
            fourierCoeDatum.setUserOverrideProperty("fc_a'46", _4_DIV_PI * B_vect[11]);
        }

        return fourierCoeDatum;
    }

}