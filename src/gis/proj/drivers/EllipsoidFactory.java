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

import gis.proj.Ellipsoid;
import gis.proj.drivers.Snyder;

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