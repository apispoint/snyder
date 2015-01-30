/*
Snyder Projection Implementation

Copyright (c) 2012-2015, APIS Point, LLC

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
package gis.proj.drivers;

import java.util.Map;

public class Fiducial {

    private String ellipsoid;
    private Map<String, String> datum;
    private String projection;

    private String lon;
    private String lat;

    private String elon;
    private String elat;

    private double ex;
    private double ey;

    private double clon;
    private double clat;

    private double cx;
    private double cy;

    private double dx;
    private double dy;
    private double dlat;
    private double dlon;

    /**
     * Default tolerance is 1.0E-7 (Degrees, ~1.11 cm)
     * Default tolerance is 1.0E-2 (Statue,  ~0.01 meters)
     * 
     * @see #setPassed(
     * 		double cx,    double cy,
     * 		double clon,  double clat,
     * 		double toleranceDegree, double toleranceStatute)
     */
    protected void setPassed(
            double cx,    double cy,
            double clon,  double clat) {
        setPassed(
                cx,     cy,
                clon,   clat,
                1.0E-7, 1.0E-2);
    }

    protected void setPassed(
            double cx,    double cy,
            double clon,  double clat,
            double toleranceDegrees,
            double toleranceStatute) {

        double e_x = getEX(); // Fiducial (expected) x
        double e_y = getEY(); // Fiducial (expected) y

        double e_lon = getELon(); // Fiducial (expected) lon
        double e_lat = getELat(); // Fiducial (expected) lat

        dx = StrictMath.abs(e_x - cx);
        dy = StrictMath.abs(e_y - cy);

        dlon = StrictMath.abs(e_lon - clon);
        dlat = StrictMath.abs(e_lat - clat);

        this.cx   = cx;
        this.clon = clon;

        this.cy   = cy;
        this.clat = clat;
    }

    public String getProjection() {
        return projection;
    }

    public void setProjection(String projection) {
        this.projection = projection;
    }

    public String getEllipsoid() {
        return ellipsoid;
    }

    public void setEllipsoid(String ellipsoid) {
        this.ellipsoid = ellipsoid;
    }

    public Map<String, String> getDatum() {
        return datum;
    }

    public void setDatum(Map<String, String> datum) {
        this.datum = datum;
    }

    public double getLon() {
        return Snyder.toDD(lon);
    }

    public void setLon(String lon) {
        this.lon = lon;
    }

    public double getELat() {
        return Snyder.toDD(elat);
    }

    public void setELat(String elat) {
        this.elat = elat;
    }

    public double getELon() {
        return Snyder.toDD(elon);
    }

    public void setELon(String elon) {
        this.elon = elon;
    }

    public double getLat() {
        return Snyder.toDD(lat);
    }

    public void setLat(String lat) {
        this.lat = lat;
    }

    public double getEX() {
        return ex;
    }

    public void setEX(double ex) {
        this.ex = ex;
    }

    public double getEY() {
        return ey;
    }

    public void setEY(double ey) {
        this.ey = ey;
    }

    public double getClon() {
        return clon;
    }

    public double getclat() {
        return clat;
    }

    public double getCx() {
        return cx;
    }

    public double getCy() {
        return cy;
    }

    public double getDx() {
        return dx;
    }

    public double getDy() {
        return dy;
    }

    public double getDlat() {
        return dlat;
    }

    public double getDlon() {
        return dlon;
    }

}