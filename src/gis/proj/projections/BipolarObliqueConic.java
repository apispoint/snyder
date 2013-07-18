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
import gis.proj.Spherical;

import java.util.Collections;
import java.util.Set;


/**
 * All formulas adopted from USGS Professional Paper 1395.
 *
 * References:
 *
 * Map Projections - A Working Manual, USGS Professional Paper 1395
 * John P. Snyder
 * Pages 116-123
 *
 */

public final class BipolarObliqueConic implements Conic, Spherical {

    private static final double D_N_110  = -110.0 * DEG_TO_RAD;
    private static final double D_N_020  = -020.0 * DEG_TO_RAD;
    private static final double D_P_020  =  020.0 * DEG_TO_RAD;
    private static final double D_P_026  =  026.0 * DEG_TO_RAD;
    private static final double D_P_031  =  031.0 * DEG_TO_RAD;
    private static final double D_P_045  =  045.0 * DEG_TO_RAD;
    private static final double D_P_052  =  052.0 * DEG_TO_RAD;
    private static final double D_P_073  =  073.0 * DEG_TO_RAD;
    private static final double D_P_104  =  104.0 * DEG_TO_RAD;
    private static final double D_P_110  =  110.0 * DEG_TO_RAD;

    private static final double D_P_031H =  031.0 * 0.5 * DEG_TO_RAD;
    private static final double D_P_073H =  073.0 * 0.5 * DEG_TO_RAD;

    private static final double cos_N_020 = StrictMath.cos(D_N_020);
    private static final double cos_P_020 = StrictMath.cos(D_P_020);
    private static final double cos_P_045 = StrictMath.cos(D_P_045);
    private static final double cos_P_104 = StrictMath.cos(D_P_104);

    private static final double sin_N_020 = StrictMath.sin(D_N_020);
    private static final double sin_P_031 = StrictMath.sin(D_P_031);
    private static final double sin_P_045 = StrictMath.sin(D_P_045);
    private static final double sin_P_052 = StrictMath.sin(D_P_052);
    private static final double sin_P_073 = StrictMath.sin(D_P_073);
    private static final double sin_P_104 = StrictMath.sin(D_P_104);

    private static final double tan_P_026 = StrictMath.tan(D_P_026);

    private static final double tan_P_031H = StrictMath.tan(D_P_031H);
    private static final double tan_P_073H = StrictMath.tan(D_P_073H);

    private static final double tan_P_026_pown;
    private static final double tan_P_031H_pown;

    private static final double lonB, AzAB, AzBA, T, zc, latc, AzC, cosAzC, sinAzC, n, one_OVER_n;

	public String getName() {
        return "Bipolar Oblique Conic";
    }

    public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
        double R   = ellip.getProperty("R");

        double[] lon = new double[x.length];
        double[] lat = new double[y.length];

        double F0   = R * sin_P_031 / (n * tan_P_031H_pown);
        double k0   = 2.0 / (1.0 + n * F0 * tan_P_026_pown / (R * sin_P_052));
        double F    = k0 * F0;
        double rhoc = 0.5 * F  * T;

        double ZB, AzB, rhoB, alpha, rhoBp, xp, yp, ZA, AzA, rhoA, rhoAp, AzBp, AzAp;
        double cosAzX;

        double rhoX_old;
        int j;

        for(int i = 0; i < lon.length; ++i) {
        	xp = -x[i] * cosAzC + y[i] * sinAzC;
        	yp = -x[i] * sinAzC - y[i] * cosAzC;

        	if(xp < 0.0) {
        		// (17-45)
                rhoAp = StrictMath.hypot(xp, rhoc - yp);
                AzAp  = StrictMath.atan2(xp, rhoc - yp);

                rhoA = rhoAp;
                j = 0;
                do {
                    rhoX_old = rhoA;
                    ZA = 2.0 * StrictMath.atan(StrictMath.pow(rhoA / F, one_OVER_n));
                    alpha = StrictMath.acos(
                            (StrictMath.pow(StrictMath.tan(0.5 * ZA), n) + StrictMath.pow(StrictMath.tan(0.5 * (D_P_104 - ZA)), n))
                            / T);

                    if(j == 0 && alpha < StrictMath.abs(AzAp))
                        break;

                    rhoA = rhoAp * StrictMath.cos(alpha + AzAp);
                } while(++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(rhoA - rhoX_old));

                AzA = AzAB - AzAp * one_OVER_n;
                cosAzX = StrictMath.cos(AzA);

                lon[i] =
                		normalizeLonRad(StrictMath.atan2(
                		StrictMath.sin(AzA),
                		cos_N_020 / StrictMath.tan(ZA) - sin_N_020 * cosAzX) + D_N_110);

                lat[i] = StrictMath.asin(sin_N_020 * StrictMath.cos(ZA) + cos_P_020 * StrictMath.sin(ZA) * cosAzX);
        	} else {
        		// (17-36)
        		rhoBp = StrictMath.hypot(xp, rhoc + yp);
        		AzBp = StrictMath.atan2(xp, rhoc + yp);

        		rhoB = rhoBp;
        		j = 0;
        		do {
        			rhoX_old = rhoB;
        			ZB = 2.0 * StrictMath.atan(StrictMath.pow(rhoB / F, one_OVER_n));
        			alpha = StrictMath.acos(
        					(StrictMath.pow(StrictMath.tan(0.5 * ZB), n) + StrictMath.pow(StrictMath.tan(0.5 * (D_P_104 - ZB)), n))
        					/ T);

        			if(j == 0 && alpha < StrictMath.abs(AzBp))
                        break;

            		rhoB = rhoBp * StrictMath.cos(alpha - AzBp);
        		} while(++j < SERIES_EXPANSION_LIMIT && NEAR_ZERO_RAD < StrictMath.abs(rhoB - rhoX_old));

        		AzB = AzBA - AzBp / n;
        		cosAzX = StrictMath.cos(AzB);

        		lon[i] =
        				normalizeLonRad(lonB - StrictMath.atan2(
                		StrictMath.sin(AzB),
                		cos_P_045 / StrictMath.tan(ZB) - sin_P_045 * cosAzX));

                lat[i] = StrictMath.asin(sin_P_045 * StrictMath.cos(ZB) + cos_P_045 * StrictMath.sin(ZB) * cosAzX);
        	}
        }

        return new double[][] {lon, lat};
    }

    public double[][] forward(double[] lon, double[] lat, Ellipsoid ellip, Datum datum) {
        double R   = ellip.getProperty("R");

        double[] x = new double[lon.length];
        double[] y = new double[lat.length];

        double F0   = R * sin_P_031 / (n * tan_P_031H_pown);
        double k0   = 2.0 / (1.0 + n * F0 * tan_P_026_pown / (R * sin_P_052));
        double F    = k0 * F0;
        double rhoc = 0.5 * F  * T;

        double ZB, AzB, rhoB, alpha, rhoBp, xp, yp, ZA, AzA, rhoA, rhoAp;
        double lonB_M_lon, coslonB_M_lon, tan_half_ZX_pown, tan_half_104_M_ZX_pown, n_AzXX_M_AzX;

        double coslat, sinlat, tanlat, lon_110, coslon_110;

        for(int i = 0; i < lon.length; ++i) {
        	lonB_M_lon = normalizeLonRad(lonB - lon[i]);
        	coslonB_M_lon = StrictMath.cos(lonB_M_lon);

        	coslat = StrictMath.cos(lat[i]);
        	sinlat = StrictMath.sin(lat[i]);
        	tanlat = StrictMath.tan(lat[i]);

        	ZB = StrictMath.acos(sin_P_045 * sinlat + cos_P_045 * coslat * coslonB_M_lon);

        	AzB = StrictMath.atan2(
        			StrictMath.sin(lonB_M_lon) ,
        			(cos_P_045 * tanlat - sin_P_045 * coslonB_M_lon));

            if(AzB > AzBA) {
            	// (17-23)
            	lon_110 = lon[i] + D_P_110;
            	coslon_110 = StrictMath.cos(lon_110);

            	ZA = StrictMath.acos(sin_N_020 * sinlat + cos_N_020 * coslat * coslon_110);

            	AzA = StrictMath.atan2(
            			StrictMath.sin(lon_110),
            			cos_N_020 * tanlat - sin_N_020 * coslon_110
            			);

            	tan_half_ZX_pown       = StrictMath.pow(StrictMath.tan(0.5 * ZA), n);
            	tan_half_104_M_ZX_pown = StrictMath.pow(StrictMath.tan(0.5 * (D_P_104 - ZA)), n);

            	rhoA  = F * tan_half_ZX_pown;
            	alpha = StrictMath.acos((tan_half_ZX_pown + tan_half_104_M_ZX_pown) / T);

            	rhoAp = rhoA;

            	n_AzXX_M_AzX = n * (AzAB - AzA);

            	if(StrictMath.abs(n_AzXX_M_AzX) < alpha)
            		rhoAp = rhoA / StrictMath.cos(alpha + n_AzXX_M_AzX);

            	xp =  rhoAp * StrictMath.sin(n_AzXX_M_AzX);
            	yp = -rhoAp * StrictMath.cos(n_AzXX_M_AzX) + rhoc;
            }
            else {
            	// (17-16)
            	tan_half_ZX_pown       = StrictMath.pow(StrictMath.tan(0.5 * ZB), n);
            	tan_half_104_M_ZX_pown = StrictMath.pow(StrictMath.tan(0.5 * (D_P_104 - ZB)), n);

            	rhoB = F * tan_half_ZX_pown;
            	alpha = StrictMath.acos((tan_half_ZX_pown + tan_half_104_M_ZX_pown) / T);

            	rhoBp = rhoB;

            	n_AzXX_M_AzX = n * (AzBA - AzB);

            	if(StrictMath.abs(n_AzXX_M_AzX) < alpha)
            		rhoBp = rhoB / StrictMath.cos(alpha - n_AzXX_M_AzX);

            	xp = rhoBp * StrictMath.sin(n_AzXX_M_AzX);
            	yp = rhoBp * StrictMath.cos(n_AzXX_M_AzX) - rhoc;
            }

            x[i] = -xp * cosAzC - yp * sinAzC;
            y[i] = -yp * cosAzC + xp * sinAzC;
        }

        return new double[][] {x, y};
    }

    public Set<String> getDatumProperties() {
    	return Collections.emptySet();
    }
 
    static {
    	lonB =
        		normalizeLonRad(D_N_110 + StrictMath.acos(
        		(cos_P_104 - sin_N_020 * sin_P_045) / (cos_N_020 * cos_P_045)
        		));

        n = 
        		(StrictMath.log(sin_P_031) - StrictMath.log(sin_P_073)) /
        		(StrictMath.log(tan_P_031H) - StrictMath.log(tan_P_073H));

    	double cos_lonb_P_D_P_110 = StrictMath.cos(normalizeLonRad(lonB + D_P_110));

        AzAB = StrictMath.acos(
        		(
        				cos_N_020 * sin_P_045 -
        				sin_N_020 * cos_P_045 * cos_lonb_P_D_P_110
        		) / sin_P_104
        		);

        AzBA = StrictMath.acos(
        		(
        				cos_P_045 * sin_N_020 -
        				sin_P_045 * cos_N_020 * cos_lonb_P_D_P_110
        		) / sin_P_104
        		);

        T    = 
        		(tan_P_031H_pown = StrictMath.pow(tan_P_031H, n)) + StrictMath.pow(tan_P_073H, n);

        zc   = 2.0 * StrictMath.atan(StrictMath.pow(T * 0.5, 1.0 / n));

        latc = StrictMath.asin(sin_N_020 * StrictMath.cos(zc) + cos_N_020 * StrictMath.sin(zc) * StrictMath.cos(AzAB));
        AzC  = StrictMath.asin(cos_N_020 * StrictMath.sin(AzAB) / StrictMath.cos(latc));

        cosAzC = StrictMath.cos(AzC);
        sinAzC = StrictMath.sin(AzC);

        tan_P_026_pown = StrictMath.pow(tan_P_026, n);

        one_OVER_n = 1.0 / n;
    }

}