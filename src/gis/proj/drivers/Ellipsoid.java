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