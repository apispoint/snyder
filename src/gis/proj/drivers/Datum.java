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
import java.util.Map.Entry;

final class Datum extends gis.proj.Datum {

    protected Datum() {}

    protected void setName(String datumName) {
       this.datumName = datumName;
    }

    protected void setProperties(Map<String, String> datumProperties) {
    	for(Entry<String, String> e : datumProperties.entrySet())
    		this.datumProperties.put(e.getKey(), Snyder.parseDatumVal(e.getValue()));
    }

}