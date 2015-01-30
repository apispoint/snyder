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
import java.util.ArrayList;
import java.util.Map.Entry;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;

public final class SnyderTest {

    @Parameter(
            names = "-if", description = "Fiducial JSON file", required = true)
    private String ifName;

    @Parameter(
            names = "-of", description = "Fiducial JSON result file",
            required = true)
    private String ofName;

    private ObjectMapper om = new ObjectMapper();

    private SnyderTest() {}

    public static void main(String... args) {
        SnyderTest sfTest = new SnyderTest();
    	JCommander jc = new JCommander(sfTest);

        try {
        	jc.parse(args);
        } catch(Exception e) {
        	jc.usage();
        	System.exit(-10);
        }

        Snyder.printLicenseInformation("SnyderTest");

        double[][]  fVal = null; // forward values
        double[][]  iVal = null; // inverse values

        Ellipsoid  ellip;
        Datum      datum;
        Projection proj;

        ArrayList<Fiducial> testData = null;

        try {
            testData = Snyder.fromJSON(
                    new TypeReference<ArrayList<Fiducial>>() {},
                    sfTest.ifName);
        } catch(Exception e) {}

        if(testData != null) {
            for(Fiducial sf : testData) {
                ellip = null;
                datum = null;
                proj  = null;

                ellip = EllipsoidFactory.getInstance().getEllipsoid(
                        sf.getEllipsoid());

                datum = new Datum();

                for(Entry<String, String> entry : sf.getDatum().entrySet()) {
                    datum.setUserOverrideProperty(
                    		entry.getKey(),
                    		Snyder.parseDatumVal(entry.getValue().toLowerCase()));
                }

                try {
                    Class<?> cls = Class.forName(sf.getProjection());
                    proj = (Projection) cls.newInstance();
                } catch (Exception ex) {
                    ex.printStackTrace();
                }

                fVal = proj.forward(
                		new double[] { sf.getLon() * SnyderMath.DEG_TO_RAD },
                		new double[] { sf.getLat() * SnyderMath.DEG_TO_RAD },
                		ellip, datum);

                iVal = proj.inverse(
                		new double[] { fVal[0][0] },
                		new double[] { fVal[1][0] },
                		ellip, datum);

                iVal[0][0] *= SnyderMath.RAD_TO_DEG;
                iVal[1][0] *= SnyderMath.RAD_TO_DEG;

                sf.setPassed(
                		fVal[0][0], fVal[1][0],
                		iVal[0][0], iVal[1][0]);
            }

            try {

            	if(sfTest.ofName != null)
            		sfTest.om.writeValue(new File(sfTest.ofName), testData);

            } catch(Exception e) {
                System.out.println(
                        "Couldn't write the fiducial data results, of=" +
                        sfTest.ofName);
                System.exit(-1);
            }

            System.out.print("\n**********************************************************************\n");
            System.out.print("SNYDER TEST SUMMARY");
            System.out.print("\n**********************************************************************\n");
            System.out.println("\n\t  Total test cases : " + testData.size());
        }
        else {
            System.out.println(
                    "Couldn't load the fiducial data, if=" + sfTest.ifName);
            System.exit(-1);
        }
    }

}