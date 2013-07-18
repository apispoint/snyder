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
 * Pages 76-85
 *
 */

public final class TransverseCylindricalEA implements Cylindrical {

	public String getName() {
		return "Transverse Cylindrical Equal-Area";
	}

	public double[][] inverse(double[] x, double[] y, Ellipsoid ellip, Datum datum) {
		double a   = ellip.getProperty("a");
		double e   = ellip.getProperty("e");
		double esq = ellip.getProperty("e^2");

		double h0   = datum.getProperty("h0");
		double lon0 = datum.getProperty("lon0");
		double lat0 = datum.getProperty("lat0");

		double[] lon = new double[x.length];
		double[] lat = new double[y.length];

		if(ellip.isSphere()) {
			double ah0 = a * h0;
			double sqrt_1_M_factorsq, D, factor;

			for(int i = 0; i < lon.length; ++i) {
				D = y[i] / ah0 + lat0;
				factor = (h0 * x[i]);
				sqrt_1_M_factorsq = StrictMath.sqrt(1.0 - factor * factor);

				lon[i] = normalizeLonRad(lon0 + StrictMath.atan2(h0 * x[i] / a, sqrt_1_M_factorsq * StrictMath.cos(D)));
				lat[i] = StrictMath.asin(sqrt_1_M_factorsq * StrictMath.sin(D));
			}
		}
		else {
			double one_M_esq = 1.0 - esq;
			double one_over_2e = 1.0 / (2.0 * e);

			double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
			double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
			double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
			double M_factor6 = (esq * esq * esq) * 0.011393229167;

			double
			M0  = lat0 * M_factor0;
			M0 -= StrictMath.sin(2.0 * lat0) * M_factor2;
			M0 += StrictMath.sin(4.0 * lat0) * M_factor4;
			M0 -= StrictMath.sin(6.0 * lat0) * M_factor6;
			M0 *= a;

			double sin_hlf_pi = StrictMath.sin(PI_DIV_2);
			double esin_hlf_pi = e * sin_hlf_pi;

			double qp = one_M_esq * (
					sin_hlf_pi / (1.0 - esq * sin_hlf_pi * sin_hlf_pi) - one_over_2e * StrictMath.log((1.0 - esin_hlf_pi) / (1.0 + esin_hlf_pi))
					);

			double sqrt_1_M_esq = StrictMath.sqrt(1.0 - esq);
			double E1 = (1.0 - sqrt_1_M_esq) / (1.0 + sqrt_1_M_esq);
			double MU_DIV = a * (1.0 - esq * 0.25 - esq * esq * (3.0 / 64.0) - esq * esq * esq * (5.0 / 256.0));

            double Mc, muc, beta, betac, betap;
            double phic, sinphic, esinphic, esqsinphicsinphic, cosbetac, q;

            for(int i = 0; i < lon.length; ++i) {

                Mc = M0 + y[i] / h0;
                muc = Mc / MU_DIV;

				// (3-26)
                phic = muc + ((3.0/2.0) * E1 - (27.0 / 32.0)*(E1 * E1 * E1)) * StrictMath.sin(2.0 * muc);
                phic += (((21.0/16.0) * E1 * E1) - ((55.0/32.0) * E1 * E1 * E1 * E1)) * StrictMath.sin(4.0 * muc);
                phic += ((151.0/96.0) * E1 * E1 * E1) * StrictMath.sin(6.0 * muc);
                phic += ((1097.0/512.0) * E1 * E1 * E1 * E1) * StrictMath.sin(8.0 * muc);

                sinphic = StrictMath.sin(phic);
                esinphic = e * sinphic;
                esqsinphicsinphic = esq * sinphic * sinphic;

                // (3-12)
                q = one_M_esq * (
                        sinphic / (1.0 - esqsinphicsinphic) - one_over_2e * StrictMath.log((1.0 - esinphic) / (1.0 + esinphic))
                        );

                // (3-11)
                betac = StrictMath.asin(q / qp);
                cosbetac = StrictMath.cos(betac);

                betap = -StrictMath.asin(h0 * x[i] * cosbetac * StrictMath.sqrt(1.0 - esqsinphicsinphic) / (a * StrictMath.cos(phic)));

                beta = StrictMath.asin(StrictMath.cos(betap) * StrictMath.sin(betac));

                lon[i] = lon0 - StrictMath.atan(StrictMath.tan(betap) / cosbetac);

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
		double lon0 = datum.getProperty("lon0");
		double lat0 = datum.getProperty("lat0");

		double x[] = new double[lon.length];
		double y[] = new double[lat.length];

		if(ellip.isSphere()) {
			double a_over_h0 = a / h0;
			double ah0 = a * h0;
			double lon_M_lon0;

			for(int i = 0; i < lon.length; ++i) {
				lon_M_lon0 = normalizeLonRad(lon[i] - lon0);
				x[i] = a_over_h0 * StrictMath.cos(lat[i]) * StrictMath.sin(lon_M_lon0);
				y[i] = ah0 * (StrictMath.atan2(StrictMath.tan(lat[i]), StrictMath.cos(lon_M_lon0)) - lat0);
			}
		}
		else {
			double one_M_esq = 1.0 - esq;
			double one_over_2e = 1.0 / (2.0 * e);

			double sin_hlf_pi = StrictMath.sin(PI_DIV_2);
			double esin_hlf_pi = e * sin_hlf_pi;

			double qp = one_M_esq * (
					sin_hlf_pi / (1.0 - esq * sin_hlf_pi * sin_hlf_pi) - one_over_2e * StrictMath.log((1.0 - esin_hlf_pi) / (1.0 + esin_hlf_pi))
					);

			double M_factor0 = 1.0 - esq * 0.25 - (esq * esq) * 0.046875 - (esq * esq * esq) * 0.01953125;
			double M_factor2 = esq * 0.375 + (esq * esq) * 0.09375 + (esq * esq * esq) * 0.0439453125;
			double M_factor4 = (esq * esq) * 0.05859375 + (esq * esq * esq) * 0.0439453125;
			double M_factor6 = (esq * esq * esq) * 0.011393229167;

			double
			M0  = lat0 * M_factor0;
			M0 -= StrictMath.sin(2.0 * lat0) * M_factor2;
			M0 += StrictMath.sin(4.0 * lat0) * M_factor4;
			M0 -= StrictMath.sin(6.0 * lat0) * M_factor6;
			M0 *= a;

			double Mc, phic, lon_M_lon0, q, esqsinlatsinlat, esinlat, beta, betac, sinlat, sinphic;

			for(int i = 0; i < lon.length; ++i) {
				sinlat = StrictMath.sin(lat[i]);
				esinlat = e * sinlat;
				esqsinlatsinlat = esq * sinlat * sinlat;

                lon_M_lon0 = normalizeLonRad(lon[i] - lon0);

				//(3-12)
				q = one_M_esq * (
						sinlat / (1.0 - esqsinlatsinlat) - one_over_2e * StrictMath.log((1.0 - esinlat) / (1.0 + esinlat))
						);

                // Authalic latitude

				// (3-11)
				beta = StrictMath.asin(q / qp);

				betac = StrictMath.atan2(StrictMath.tan(beta), StrictMath.cos(lon_M_lon0));

				// (3-18)
				phic = betac + (esq / 3.0 + (esq *esq) * (31.0/180.0) + (esq*esq*esq) * (517.0/5040.0)) * StrictMath.sin(2.0 * betac) +
						((esq *esq) * (23.0/360.0) + (esq*esq*esq) * (251.0/3780.0)) * StrictMath.sin(4.0 * betac) +
						((esq*esq*esq) * (761.0/45360.0)) * StrictMath.sin(6.0 * betac);

				Mc  = phic * M_factor0;
				Mc -= StrictMath.sin(2.0 * phic) * M_factor2;
				Mc += StrictMath.sin(4.0 * phic) * M_factor4;
				Mc -= StrictMath.sin(6.0 * phic) * M_factor6;
				Mc *= a;

				sinphic = StrictMath.sin(phic);
				x[i] = a * StrictMath.cos(beta) * StrictMath.cos(phic) * StrictMath.sin(lon_M_lon0) /
						(h0 * StrictMath.cos(betac) * StrictMath.sqrt(1.0 - esq * sinphic * sinphic));

				y[i] = h0 * (Mc - M0);
			}
		}

		return new double[][] {x, y};
	}

	public Set<String> getDatumProperties() {
        return new HashSet<String>(Arrays.asList(new String[]{
                "h0",
                "lon0",
                "lat0"}));
	}

}