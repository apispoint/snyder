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

import java.util.Map;
import java.util.Set;

final class Ellipsoid extends gis.proj.Ellipsoid {

    protected Ellipsoid() {
        super(1.0, 1.0, null, null, null);
    }

    protected void setIsSphere(boolean isSphere) {
        this.isSphere = isSphere;
    }

    protected void setName(String ellipsoidName) {
       this.name = ellipsoidName;
    }

    protected void setId(String id) {
        this.id = id;
    }

    protected void setDescription(String desc) {
        this.description = desc;
    }

    protected void setProperties(Map<String, Double> ellipsoidProperties) {
        // Create a "complete" ellipsoid by merging the curated JSON ellipsoid properties
        // with the properties from an instantiated ellipsoid instance.
        //
        // Once the properties are merged the temporary ellipsoid instance is discarded.
        //
        // This allows for easy overriding of computed ellipsoid properties because
        // user defined properties take precedence over the computed properties.
        properties = ellipsoidProperties;

        double a        = properties.get("a");
        double inverseF = properties.get("1/f");

        gis.proj.Ellipsoid computedEllipsoid =
                new gis.proj.Ellipsoid(a, inverseF, name, id, description);

        // Set difference
        Set<String> props = computedEllipsoid.getPropertyNames();
        props.removeAll(properties.keySet());

        // Add only the missing properties
        for(String k : props)
            this.properties.put(k, computedEllipsoid.getProperty(k));
    }

}