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

import static gis.proj.SnyderMath.*;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class Ellipsoid implements Property {

    protected String  name        = null;
    protected String  description = null;
    protected String  id          = null;
    protected boolean isSphere    = false;

    // Properties
    // e.g. a
    // e.g. b
    // e.g. e
    protected Map<String, Double> properties = new HashMap<String, Double>();

    public Ellipsoid(double a, double inverseF) {
        this(a, inverseF, null, null, null);
    }

    public Ellipsoid(double a, double inverseF, String name, String id, String description) {
        double f = 1.0 / inverseF;
        double b = a;
        double n = f / (2.0 - f); // f'' 3rd flattening or (a-b)/(a+b)

        // Handle spherical cases
        if(Double.isNaN(f) == false && Double.isInfinite(f) == false)
            b = a * (1.0 - f);
        else
            f = n = 0.0;

        double e  = StrictMath.sqrt(1.0 - (b * b) / (a * a));
        double ep = e / StrictMath.sqrt(1.0 - (e * e));

        // Precalculate common values
        properties.put("a",    a);
        properties.put("1/f",  inverseF);
        properties.put("b",    b);
        properties.put("R",    a);
        properties.put("a^2",  a * a);
        properties.put("b^2",  b * b);
        properties.put("f",    f);
        properties.put("n",    n);
        properties.put("e",    e);
        properties.put("e^2",  e * e);
        properties.put("e^3",  e * e * e);
        properties.put("e^4",  e * e * e * e);
        properties.put("e^5",  e * e * e * e * e);
        properties.put("e^6",  e * e * e * e * e * e);
        properties.put("e'",   ep);
        properties.put("e'^2", ep * ep);

        isSphere         = e <= NEAR_ZERO_DEG;
        this.name        = name;
        this.id          = id;
        this.description = description;
    }

    public boolean isSphere() {
        return isSphere;
    }

    public String getName() {
        return name;
    }

    public String getId() {
        return id;
    }

    public String getDescription() {
        return description;
    }

    public double getProperty(String property) {
        Double ret = properties.get(property);

        return ret == null ? Double.NaN : ret;
    }

    public Set<String> getPropertyNames() {
        return properties.keySet();
    }

}