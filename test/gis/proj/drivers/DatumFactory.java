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
package gis.proj.drivers;

import gis.proj.Datum;
import gis.proj.drivers.Snyder;

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