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
package gis.proj;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class Datum implements Property {

    protected String datumName = null;

    // Properties
    // e.g. k0
    // e.g. lon0
    // e.g. lat1
    protected Map<String, Double> datumProperties    = new HashMap<String, Double>();
    protected Map<String, Double> overrideProperties = new HashMap<String, Double>();

    public String getName() {
        return datumName;
    }

    public double getProperty(String property) {
        Double ret =
                overrideProperties.containsKey(property) ?
                        overrideProperties.get(property) :
                        datumProperties.get(property);

        return ret == null ? Double.NaN : ret;
    }

    public Set<String> getPropertyNames() {
        Set<String> keys = new HashSet<String>(datumProperties.keySet());
        keys.addAll(overrideProperties.keySet());
        return keys;
    }

    public Set<String> getUserOverrideProperty() {
        return overrideProperties.keySet();
    }

    public void setUserOverrideProperty(String property, double value) {
        overrideProperties.put(property, value);
    }

}