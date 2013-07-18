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
 * Pages 203-212
 *
 */

public final class ModifiedStereographic implements Azimuthal {

    //
    // [0][] >> {lat1, g0, A1, .., Am}
    // [1][] >> {lon0, g0, B1, .., Bm}
    //

	// Miller lat1, lon0 taken from 1994 reprint & corrections
	private final static double[][] MILLER_S = new double [][] {
		{18.0 * DEG_TO_RAD, 0.0, 0.924500, 0.0, 0.019430},
		{18.0 * DEG_TO_RAD, 0.0, 0.0     , 0.0, 0.0     }
	};

	private final static double[][] LEE_S = new double [][] {
		{-010.0 * DEG_TO_RAD, 0.0, 0.721316, 0.0, -0.00881625},
		{-165.0 * DEG_TO_RAD, 0.0, 0.0     , 0.0, -0.00617325}
	};

	private final static double[][] GS50_S = new double [][] {
		{  45.0 * DEG_TO_RAD, 0.0, 0.9842990, 0.0211642, -0.1036018, -0.0329095, 0.0499471, 0.0260460,  0.0007388,  0.0075848, -0.0216473, -0.0225161},
		{-120.0 * DEG_TO_RAD, 0.0, 0.0      , 0.0037608, -0.0575102, -0.0320119, 0.1223335, 0.0899805, -0.1435792, -0.1334108,  0.0776645,  0.0853673}
	};

	private final static double[][] ALASKA_S = new double [][] {
		{  64.0 * DEG_TO_RAD, 0.0, 0.9972523,  0.0052513, 0.0074606, -0.0153783,  0.0636871,  0.3660976},
		{-152.0 * DEG_TO_RAD, 0.0, 0.0      , -0.0041175, 0.0048125, -0.1968253, -0.1408027, -0.2937382}
	};

	private final static double[][] GS48_S = new double [][] {
		{ 39.0 * DEG_TO_RAD, 0.0, 0.98879, 0.0, -0.050909, 0.0, 0.075528},
		{-96.0 * DEG_TO_RAD, 0.0, 0.0    , 0.0,  0.0     , 0.0, 0.0     }
	};

	private final static double[][] GS50_E = new double [][] {
		{  45.0 * DEG_TO_RAD, 0.0, 0.9827497, 0.0210669, -0.1031415, -0.0323337, 0.0502303, 0.0251805, -0.0012315,  0.0072202, -0.0194029, -0.0210072},
		{-120.0 * DEG_TO_RAD, 0.0, 0.0      , 0.0053804, -0.0571664, -0.0322847, 0.1211983, 0.0895678, -0.1416121, -0.1317091,  0.0759677,  0.0834037}
	};

	private final static double[][] ALASKA_E = new double [][] {
		{  64.0 * DEG_TO_RAD, 0.0, 0.9945303,  0.0052083, 0.0072721, -0.0151089,  0.0642675,  0.3582802},
		{-152.0 * DEG_TO_RAD, 0.0, 0.0      , -0.0027404, 0.0048181, -0.1932526, -0.1381226, -0.2884586}
	};

    public String getName() {
        return "Modified Stereographic";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a = ellip.getProperty("a");
        double e = ellip.getProperty("e");

        double[][] consts = ALASKA_E;

        if     (datum.getProperty("useMiller")   == 1.0) consts = MILLER_S;
        else if(datum.getProperty("useLee")      == 1.0) consts = LEE_S;
        else if(datum.getProperty("useGS50_S")   == 1.0) consts = GS50_S;
        else if(datum.getProperty("useAlaska_S") == 1.0) consts = ALASKA_S;
        else if(datum.getProperty("useGS48")     == 1.0) consts = GS48_S;
        else if(datum.getProperty("useGS50_E")   == 1.0) consts = GS50_E;

        int m       = consts[0].length - 1; // used by aM, bM
        int n       = m - 1;                // used by cN, dN

        double lat1 = consts[0][0];
        double lon0 = consts[1][0];

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double hlf_e = e * 0.5;

        double sinlat1 = StrictMath.sin(lat1);
        double esinlat1 = e * sinlat1;

        // (3-1)
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

        double coschi1 = StrictMath.cos(chi1);
        double sinchi1 = StrictMath.sin(chi1);

        double rho, c, chi, cosc, sinc, xp, yp, phi, phi_old, esinphi;
        double oxp, oyp, xp_old, yp_old;

        int j, k;

        // To convert from the rectangular coordinates (xp, yp, with Radius 1 unit)
        // to the base map (x, y) use of complex numbers are required
        //
        // [0][] real
        // [1][] imaginary
        //
        double[][] aM = new double[2][m];
        double[][] bM = new double[2][m];
        double[][] cM = new double[2][m];
        double[][] dM = new double[2][m];

        double r, sp, dxp, dyp, d_denom;
        double[] f = new double [2], F = new double[2];

        for(int i = 0; i < x.length; ++i) {
            oxp = xp = x[i] / a;
            oyp = yp = y[i] / a;

            j = 0;
            do {
            	xp_old = xp;
            	yp_old = yp;

            	r = 2.0 * xp;
                sp = xp * xp + yp * yp;

                aM[0][0] = consts[0][m];
                aM[1][0] = consts[1][m];

                bM[0][0] = consts[0][n];
                bM[1][0] = consts[1][n];

                cM[0][0] = n * aM[0][0];
                cM[1][0] = n * aM[1][0];

                dM[0][0] = (n - 1) * bM[0][0];
                dM[1][0] = (n - 1) * bM[1][0];

                for(k = 1; k < m; k++) {
                    aM[0][k] = bM[0][k - 1] + r * aM[0][k - 1];
                    aM[1][k] = bM[1][k - 1] + r * aM[1][k - 1];

                    // In the last iteration these values are good because
                    // g0 is in that position (usually 0.0)
                    bM[0][k] = consts[0][n - k] - sp * aM[0][k - 1];
                    bM[1][k] = consts[1][n - k] - sp * aM[1][k - 1];

                    cM[0][k] = dM[0][k - 1] + r * cM[0][k - 1];
                    cM[1][k] = dM[1][k - 1] + r * cM[1][k - 1];

                    dM[0][k] = (n - k - 1) * consts[0][n - k] - sp * cM[0][k - 1];
                    dM[1][k] = (n - k - 1) * consts[1][n - k] - sp * cM[1][k - 1];
                }

                // (26-11)
            	f[0] = xp * aM[0][n - 1] - yp * aM[1][n - 1] + bM[0][n - 1] - oxp;
            	f[1] = xp * aM[1][n - 1] + yp * aM[0][n - 1] + bM[1][n - 1] - oyp;

            	// (26-8)
            	F[0] = xp * cM[0][n - 2] - yp * cM[1][n - 2] + dM[0][n - 2];
            	F[1] = xp * cM[1][n - 2] + yp * cM[0][n - 2] + dM[1][n - 2];

            	d_denom = (F[0] * F[0] + F[1] * F[1]);
            	dxp = (f[0] * F[0] + f[1] * F[1]) / d_denom;
            	dyp = (f[1] * F[0] - f[0] * F[1]) / d_denom;

            	xp -= dxp;
            	yp -= dyp;
            } while (++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_DEG < StrictMath.abs(xp - xp_old) && NEAR_ZERO_DEG < StrictMath.abs(yp - yp_old));

            rho = StrictMath.hypot(xp, yp);

            if(NEAR_ZERO_DEG <= rho) {
                c = 2.0 * StrictMath.atan(rho * 0.5);
                cosc = StrictMath.cos(c);
                sinc = StrictMath.sin(c);

                chi = StrictMath.asin(cosc * sinchi1 + (yp * sinc * coschi1 / rho));

                lon[i] = lon0 + StrictMath.atan2(xp * sinc, rho * coschi1 * cosc - yp * sinchi1 * sinc);
            }
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
            } while (++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD <= StrictMath.abs(phi - phi_old));

            lon[i] = normalizeLonRad(lon[i]);
            lat[i] = phi;
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a = ellip.getProperty("a");
        double e = ellip.getProperty("e");

        double[][] consts = ALASKA_E;

        if     (datum.getProperty("useMiller")   == 1.0) consts = MILLER_S;
        else if(datum.getProperty("useLee")      == 1.0) consts = LEE_S;
        else if(datum.getProperty("useGS50_S")   == 1.0) consts = GS50_S;
        else if(datum.getProperty("useAlaska_S") == 1.0) consts = ALASKA_S;
        else if(datum.getProperty("useGS48")     == 1.0) consts = GS48_S;
        else if(datum.getProperty("useGS50_E")   == 1.0) consts = GS50_E;

    	int m       = consts[0].length - 1; // used by aM, bM
    	int n       = m - 1;                // used by cN, dN

    	double lat1 = consts[0][0];
    	double lon0 = consts[1][0];

    	double x[] = new double[lon.length];
        double y[] = new double[lat.length];

        double hlf_e = e * 0.5;

        double esinlat1 = e * StrictMath.sin(lat1);

        // (3-1)
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

        double coschi1 = StrictMath.cos(chi1);
        double sinchi1 = StrictMath.sin(chi1);

        double esinlat, chi, s, xp, yp, lon_M_lon0, coslon_M_lon0, coschi, sinchi;
        double r, sp;

        // To convert from the rectangular coordinates (xp, yp, with Radius 1 unit)
        // to the base map (x, y) use of complex numbers are required
        //
        // [0][] real
        // [1][] imaginary
        //
        double[][] aM = new double[2][m];
        double[][] bM = new double[2][m];

        for(int i = 0; i < lon.length; ++i) {
        	esinlat = e * StrictMath.sin(lat[i]);

            // (3-1)
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

        	lon_M_lon0 = normalizeLonRad(lon[i] - lon0);
        	coslon_M_lon0 = StrictMath.cos(lon_M_lon0);

        	s = 2.0 / (1.0 + sinchi1 * sinchi + coschi1 * coschi * coslon_M_lon0);
        	xp = s * coschi * StrictMath.sin(lon_M_lon0);
        	yp = s * (coschi1 * sinchi - sinchi1 * coschi * coslon_M_lon0);

        	r = 2.0 * xp;
    		sp = xp * xp + yp * yp;

    		aM[0][0] = consts[0][m];
    		aM[1][0] = consts[1][m];

    		bM[0][0] = consts[0][n];
    		bM[1][0] = consts[1][n];

        	for(int j = 1; j < m; j++) {
        		aM[0][j] = bM[0][j - 1] + r * aM[0][j - 1];
        		aM[1][j] = bM[1][j - 1] + r * aM[1][j - 1];

        		// In the last iteration these values are good because
        		// g0 is in that position (usually 0.0)
        		bM[0][j] = consts[0][n - j] - sp * aM[0][j - 1];
        		bM[1][j] = consts[1][n - j] - sp * aM[1][j - 1];
        	}

    		x[i] = a * (xp * aM[0][n - 1] - yp * aM[1][n - 1] + bM[0][n - 1]);
    		y[i] = a * (xp * aM[1][n - 1] + yp * aM[0][n - 1] + bM[1][n - 1]);
        }

        return new double[][] {x, y};
    }
 
    public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "useMiller",
                "useLee",
                "useGS50_S",
                "useAlaska_S",
                "useGS48",
                "useGS50_E"}));
    }

}