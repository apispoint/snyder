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
import gis.proj.ComplexNumber;
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

//
// ModifiedStereographic implementation uses arrays, and this (ALT)
// version uses the ComplexNumber class
// 
//
public final class ModifiedStereographicAlt implements Azimuthal {

    private final static ComplexNumber g0 = new ComplexNumber();

    private final static ComplexNumber[] MILLER_S = new ComplexNumber[] {
        new ComplexNumber(18.0 * DEG_TO_RAD, 18.0 * DEG_TO_RAD),
        g0,
        new ComplexNumber(0.924500, 0.0),
        new ComplexNumber(0.0,      0.0),
        new ComplexNumber(0.019430, 0.0)
    };

    private final static ComplexNumber[] LEE_S = new ComplexNumber[] {
        new ComplexNumber(-10.0 * DEG_TO_RAD, -165.0 * DEG_TO_RAD),
        g0,
        new ComplexNumber( 0.721316,    0.0),
        new ComplexNumber( 0.0,         0.0),
        new ComplexNumber(-0.00881625, -0.00617325)
    };

    private final static ComplexNumber[] GS50_S = new ComplexNumber[] {
        new ComplexNumber(45.0 * DEG_TO_RAD, -120.0 * DEG_TO_RAD),
        g0,
        new ComplexNumber( 0.9842990,  0.0),
        new ComplexNumber( 0.0211642,  0.0037608),
        new ComplexNumber(-0.1036018, -0.0575102),
        new ComplexNumber(-0.0329095, -0.0320119),
        new ComplexNumber( 0.0499471,  0.1223335),
        new ComplexNumber( 0.0260460,  0.0899805),
        new ComplexNumber( 0.0007388, -0.1435792),
        new ComplexNumber( 0.0075848, -0.1334108),
        new ComplexNumber(-0.0216473,  0.0776645),
        new ComplexNumber(-0.0225161,  0.0853673)
    };

    private final static ComplexNumber[] ALASKA_S = new ComplexNumber[] {
        new ComplexNumber(64.0 * DEG_TO_RAD, -152.0 * DEG_TO_RAD),
        g0,
        new ComplexNumber( 0.9972523,  0.0),
        new ComplexNumber( 0.0052513, -0.0041175),
        new ComplexNumber( 0.0074606,  0.0048125),
        new ComplexNumber(-0.0153783, -0.1968253),
        new ComplexNumber( 0.0636871, -0.1408027),
        new ComplexNumber( 0.3660976, -0.2937382)
    };

	private final static ComplexNumber[] GS48_S = new ComplexNumber[] {
        new ComplexNumber(39.0 * DEG_TO_RAD, -96.0 * DEG_TO_RAD),
        g0,
        new ComplexNumber( 0.98879,  0.0),
        new ComplexNumber( 0.0,      0.0),
        new ComplexNumber(-0.050909, 0.0),
        new ComplexNumber( 0.0,      0.0),
        new ComplexNumber( 0.085528, 0.0)
    };

    private final static ComplexNumber[] GS50_E = new ComplexNumber[] {
        new ComplexNumber(45.0 * DEG_TO_RAD, -120.0 * DEG_TO_RAD),
        g0,
        new ComplexNumber( 0.9827497,  0.0),
        new ComplexNumber( 0.0210669,  0.0053804),
        new ComplexNumber(-0.1031415, -0.0571664),
        new ComplexNumber(-0.0323337, -0.0322847),
        new ComplexNumber( 0.0502303,  0.1211983),
        new ComplexNumber( 0.0251805,  0.0895678),
        new ComplexNumber(-0.0012315, -0.1416121),
        new ComplexNumber( 0.0072202, -0.1317091),
        new ComplexNumber(-0.0194029,  0.0759677),
        new ComplexNumber(-0.0210072,  0.0834037)
    };

    private final static ComplexNumber[] ALASKA_E = new ComplexNumber[] {
        new ComplexNumber(64.0 * DEG_TO_RAD, -152.0 * DEG_TO_RAD),
        g0,
        new ComplexNumber( 0.9945303,  0.0),
        new ComplexNumber( 0.0052083, -0.0027404),
        new ComplexNumber( 0.0072721,  0.0048181),
        new ComplexNumber(-0.0151089, -0.1932526),
        new ComplexNumber( 0.0642675, -0.1381226),
        new ComplexNumber( 0.3582802, -0.2884586)
    };

    public String getName() {
        return "Modified Stereographic (ALT)";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double a = ellip.getProperty("a");
        double e = ellip.getProperty("e");

        ComplexNumber[] consts = ALASKA_E;

        if     (datum.getProperty("useMiller")   == 1.0) consts = MILLER_S;
        else if(datum.getProperty("useLee")      == 1.0) consts = LEE_S;
        else if(datum.getProperty("useGS50_S")   == 1.0) consts = GS50_S;
        else if(datum.getProperty("useAlaska_S") == 1.0) consts = ALASKA_S;
        else if(datum.getProperty("useGS48")     == 1.0) consts = GS48_S;
        else if(datum.getProperty("useGS50_E")   == 1.0) consts = GS50_E;

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double lat1 = consts[0].re();
        double lon0 = consts[0].im();

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

        double r, sp, rho, c, cosc, sinc, chi, xp, yp, phi, phi_old, esinphi;

        // To convert from the rectangular coordinates (xp, yp, with Radius 1 unit)
        // to the base map (x, y) use of complex numbers are required
        int m = consts.length - 1;
        ComplexNumber[] aM = new ComplexNumber[m];
        ComplexNumber[] bM = new ComplexNumber[m];
        ComplexNumber[] cM = new ComplexNumber[m];
        ComplexNumber[] dM = new ComplexNumber[m];

        ComplexNumber loc = new ComplexNumber();
        ComplexNumber locPrime;
        ComplexNumber f, F, locPrime_old;

        for(int i = 0; i < lon.length; ++i) {
            loc.setReal(x[i] / a);
            loc.setImaginary(y[i] / a);
            locPrime = loc.clone();

            int k = 0;
            do {
                locPrime_old = locPrime.clone();
                r = 2.0 * locPrime.re();
                sp = locPrime.absSq();

                aM[0] = consts[m];
                bM[0] = consts[m - 1];
                cM[0] = aM[0].mul(m - 1);
                dM[0] = bM[0].mul(m - 2);

                for(int j = 1; j < m; j++) {
                    aM[j] = bM[j - 1].add(aM[j - 1].mul(r));

                    // In the last iteration these values are good because
                    // g0 is in that position (usually 0.0)
                    bM[j] = consts[m - 1 - j].sub(aM[j - 1].mul(sp));

                    cM[j] = dM[j - 1].add(cM[j - 1].mul(r));
                    dM[j] = consts[m - 1 - j].mul(m - 1 - j - 1).sub(cM[j - 1].mul(sp));
                }

                f = locPrime.mul(aM[m - 2]).add(bM[m - 2]).sub(loc);
                F = locPrime.mul(cM[m - 3]).add(dM[m - 3]);
                locPrime = locPrime.sub(f.div(F));
            } while(++k < SERIES_EXPANSION_LIMIT && NEAR_ZERO_DEG < locPrime.sub(locPrime_old).absSq());

            rho = locPrime.abs();
            xp = locPrime.re();
            yp = locPrime.im();

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

            k = 0;
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
            } while (++k < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(phi - phi_old));

            lon[i] = normalizeLonRad(lon[i]);
            lat[i] = phi;
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double a = ellip.getProperty("a");
        double e = ellip.getProperty("e");

        ComplexNumber[] consts = ALASKA_E;

        if     (datum.getProperty("useMiller")   == 1.0) consts = MILLER_S;
        else if(datum.getProperty("useLee")      == 1.0) consts = LEE_S;
        else if(datum.getProperty("useGS50_S")   == 1.0) consts = GS50_S;
        else if(datum.getProperty("useAlaska_S") == 1.0) consts = ALASKA_S;
        else if(datum.getProperty("useGS48")     == 1.0) consts = GS48_S;
        else if(datum.getProperty("useGS50_E")   == 1.0) consts = GS50_E;

    	double lat1 = consts[0].re();
    	double lon0 = consts[0].im();

    	double x[] = new double[lon.length];
        double y[] = new double[lat.length];

        double hlf_e    = e * 0.5;
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

        double esinlat, chi, s, lon_M_lon0, coslon_M_lon0, coschi, sinchi;
        double r, sp;

        // To convert from the rectangular coordinates (xp, yp, with Radius 1 unit)
        // to the base map (x, y) use of complex numbers are required
        int m = consts.length - 1;
        ComplexNumber[] aM = new ComplexNumber[m];
        ComplexNumber[] bM = new ComplexNumber[m];

        ComplexNumber locPrime = new ComplexNumber();

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
        	locPrime.setReal(s * coschi * StrictMath.sin(lon_M_lon0));
        	locPrime.setImaginary(s * (coschi1 * sinchi - sinchi1 * coschi * coslon_M_lon0));

        	r = 2.0 * locPrime.re();
    		sp = locPrime.absSq();

    		aM[0] = consts[m];
    		bM[0] = consts[m - 1];

        	for(int j = 1; j < m; j++) {
        		aM[j] = bM[j - 1].add(aM[j - 1].mul(r));

        		// In the last iteration these values are good because
        		// g0 is in that position (usually 0.0)
        		bM[j] = consts[m - 1 - j].sub(aM[j - 1].mul(sp));
        	}

        	locPrime = locPrime.mul(aM[m - 2]).add(bM[m - 2]).mul(a);
        	x[i] = locPrime.re();
        	y[i] = locPrime.im();
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
                "useGS50_E",
                "useAlaska_E"}));
    }

}