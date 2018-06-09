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

import gis.proj.Datum;
import gis.proj.Ellipsoid;
import gis.proj.Projection;
import gis.proj.SnyderMath;

import java.io.File;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;

public final class Snyder {

    public static final String RESOURCE_DIR = "src/test/resources/";
    public static final String INFO_DIR     = RESOURCE_DIR + "json/";
    public static final String WMAP_DIR     = RESOURCE_DIR + "gshhg/";

    private static ObjectMapper om = new ObjectMapper();

    @Parameter(names = "-ellip", description = "Ellipsoid name", required = true)
    private String ellpStr;

    @Parameter(names = "-proj", description = "Projection classname", required = true)
    private String projStr;

    @Parameter(names = "-inv", description = "Inverse projection")
    private boolean inverse = false;

    private Snyder() {}

    protected static String licenseInformation(String prog) {
        return 
                prog +
                "  Copyright (c) 2012-2015 APIS Point, LLC.  " +
                "\nThis program comes with ABSOLUTELY NO WARRANTY. " +
                "\nThis is free software, and you are welcome to redistribute it " +
                "\nunder certain conditions.\n";
    }

    protected static void printLicenseInformation(String prog) {
        System.out.println(licenseInformation(prog));
    }

    protected static <T> T fromJSON(final TypeReference<T> type, final String filename) {
        T data = null;

        try {
            data = om.readValue(new File(filename), type);
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Bad JSON file.");
            System.exit(100);
        }

        return data;
    }

    protected static double parseDatumVal(String dVal) {
        double value = Double.NaN;

        try {
            String valueStr = dVal.trim();
            boolean doRadians = !(valueStr.lastIndexOf("r") < 0);
            double scalar = 1.0;

            //
            // Assume all values should be passed in as is, unless they
            // end with an 'r' which means the following:
            //         * Convert to radians
            //
            // All values are parsed using DMS to DD
            if(doRadians) {
                scalar = SnyderMath.DEG_TO_RAD;
                valueStr = valueStr.substring(0, valueStr.length() - 1);
            }

            value = toDD(valueStr) * scalar;
        } catch(Exception e) { }

        return value;
    }

    protected static double toDD(String dms) {
        double dd;

        dms = dms.toLowerCase();
        String[] deg = dms.split("d");

        if(deg.length > 1) {
            String[] min = deg[1].split("m");

            double sec = 0;
            if(min.length > 1)
                sec = Double.parseDouble(min[1].split("s")[0]);

            double sign = Double.parseDouble(deg[0]) < 0.0 ? -1.0 : 1.0;

            dd = (StrictMath.abs(Double.parseDouble(deg[0])) +
                    Double.parseDouble(min[0]) / 60.0 +
                    sec                        / 3600.0) * sign;
        }
        else
            dd = Double.parseDouble(deg[0]);

        return dd;
    }

    public static void main(String... args) {
        Snyder snyder = new Snyder();
        JCommander jc = new JCommander(snyder);

        try {
            jc.parse(args);
        } catch(Exception e) {
            jc.usage();
            System.exit(-10);
        }

        String fFormat = "(forward)   X: %18.9f,   Y: %18.9f%n";
        String iFormat = "(inverse) Lon: %18.9f, Lat: %18.9f%n%n";

        double[][] xy, ll = null;

        java.util.regex.Pattern pq = java.util.regex.Pattern.compile("quit", java.util.regex.Pattern.CASE_INSENSITIVE);
        java.util.Scanner s = new java.util.Scanner(System.in);

        Projection proj = null;

        printLicenseInformation("Snyder");

        try {
            System.out.println("Loading projection: " + snyder.projStr);
            Class<?> cls = Class.forName(snyder.projStr);
            proj = (Projection) cls.getConstructor().newInstance();
        } catch (Exception ex) {
            ex.printStackTrace();
        }

        double lon = 0.0, lat = 0.0, x = 0.0, y = 0.0;

        Datum d = new Datum();
        Ellipsoid e = EllipsoidFactory.getInstance().getEllipsoid(snyder.ellpStr);
        System.out.println("\nEllipsoid: " + snyder.ellpStr + ", " + e.getName() + ", " + e.getId() + ", " + e.getDescription());

        for(String prop : e.getPropertyNames()) {
            System.out.println("\t" + prop + "\t" + e.getProperty(prop));
        }

        String cmdEntryLine =
                (snyder.inverse ? "\nx y " : "\nlon lat") + ": ";

        for(String dProp : proj.getDatumProperties()) {
            d.setUserOverrideProperty(dProp, 0.0);
        }
        System.out.print(cmdEntryLine);

        while(s.hasNext(pq) == false) {
            if(snyder.inverse == false) {
                lon  = parseDatumVal(s.next());
                lat  = parseDatumVal(s.next());
            }
            else {
                x  = parseDatumVal(s.next());
                y  = parseDatumVal(s.next());
            }

            for(String dp : d.getPropertyNames()) {
                System.out.print(dp + ": ");
                d.setUserOverrideProperty(dp, parseDatumVal(s.next()));
            }

            System.out.println();

            if(snyder.inverse == false) {
                xy = proj.forward(
                        new double[] {lon},
                        new double[] {lat},
                        e, d
                        );

                System.out.printf(fFormat, xy[0][0], xy[1][0]);

                ll = proj.inverse(
                        new double[] {xy[0][0]},
                        new double[] {xy[1][0]},
                        e, d
                        );

                System.out.printf(
                        iFormat,
                        StrictMath.toDegrees(ll[0][0]),
                        StrictMath.toDegrees(ll[1][0]));
            }
            else {
                ll = proj.inverse(
                        new double[] {x},
                        new double[] {y},
                        e, d
                        );

                System.out.printf(
                        iFormat,
                        StrictMath.toDegrees(ll[0][0]),
                        StrictMath.toDegrees(ll[1][0]));

                xy = proj.forward(
                        new double[] {ll[0][0]},
                        new double[] {ll[1][0]},
                        e, d
                        );

                System.out.printf(fFormat, xy[0][0], xy[1][0]);
            }

            System.out.print(cmdEntryLine);
        }

        s.close();
    }

}