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