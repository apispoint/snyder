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