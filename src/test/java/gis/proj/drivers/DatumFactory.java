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

import java.io.File;
import java.util.Map;
import java.util.Set;

import com.fasterxml.jackson.databind.ObjectMapper;

final class DatumFactory {

    private final static DatumFactory DF;
    private final static ObjectMapper OM;

    private Map<String, Datum> datumMap;

    private DatumFactory() {}

    public static DatumFactory getInstance() {
        return DF;
    }

    public Datum getDatum(String name) {
        Datum dat = datumMap.get(name);
        if(dat == null && datumMap.containsKey(name)) {
            try {
                dat = (Datum) OM.readValue(
                        new File(Snyder.INFO_DIR + name),
                        gis.proj.drivers.Datum.class);
                datumMap.put(name, dat);
            } catch(Exception e) {
//                e.printStackTrace();
            }
        }
        return dat;
    }

    public Set<String> getDatumNames() {
        return datumMap.keySet();
    }

    protected void setDatums(Map<String, Datum> map) {
        datumMap = map;
    }

    static {
        DatumFactory tDF = null;
        ObjectMapper tom = null;

        try{
            tom = new ObjectMapper();
            tDF = tom.readValue(
                    new File(Snyder.INFO_DIR + "datums"),
                    DatumFactory.class);
        } catch(Exception e) {
            e.printStackTrace();
        }

        DF = tDF;
        OM = tom;
    }

}