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

import gis.proj.Ellipsoid;

import java.io.File;
import java.util.Map;
import java.util.Set;

import com.fasterxml.jackson.databind.ObjectMapper;

final class EllipsoidFactory {

    private final static EllipsoidFactory EF;
    private final static ObjectMapper     OM;

    private Map<String, Ellipsoid> ellipsoidMap;

    private EllipsoidFactory() {}

    public static EllipsoidFactory getInstance() {
        return EF;
    }

    public Ellipsoid getEllipsoid(String name) {
        Ellipsoid ell = ellipsoidMap.get(name);
        if(ell == null && ellipsoidMap.containsKey(name)) {
            try {
                ell = (Ellipsoid) OM.readValue(
                        new File(Snyder.INFO_DIR + name),
                        gis.proj.drivers.Ellipsoid.class);
                ellipsoidMap.put(name, ell);
            } catch(Exception e) {
                e.printStackTrace();
            }
        }
        return ell;
    }

    public Set<String> getEllipsoidNames() {
        return ellipsoidMap.keySet();
    }

    protected void setEllipsoids(Map<String, Ellipsoid> map) {
        ellipsoidMap = map;
    }

    static {
        EllipsoidFactory tef = null;
        ObjectMapper     tom = null;

        try{
            tom = new ObjectMapper();
            tef = tom.readValue(
                    new File(Snyder.INFO_DIR + "ellipsoids"),
                    EllipsoidFactory.class);
        } catch(Exception e) {
            e.printStackTrace();
        }

        EF = tef;
        OM = tom;
    }

}