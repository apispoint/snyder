Snyder Geospatial Projection System
===================================

Copyright (c) 2012-2017<br/>
APIS Point, LLC

Description
===========

Snyder is a geospatial projection library implementing all projections and the appendix fiducial values for unit testing as defined in John P. Snyder's seminal USGS professional paper, Map Projections--A Working Manual.

The library is free and open-source software licensed under The MIT License. Snyder runs on Android, is thread safe, and does not require any external dependencies.

Features
========
    + Android compatible
    + Thread safe
    + Minimum class topology
    + Small library footprint

Code Example
============

***Stereographic Projection Example***
```Java
    public static void main(String... args) {
        Projection proj = new Stereographic();
        //
        // WGS 84
        //
        Ellipsoid e = new Ellipsoid(6378137.0, 298.257223563);

        //
        // Print all available parameters
        //
        for(String prop : proj.getDatumProperties())
            System.out.println(prop);

        Datum d = new Datum();

        //
        // Core projection library assumes all long/lat values are in Radians
        //
        d.setUserOverrideProperty("lon0", 0);
        d.setUserOverrideProperty("lat1", 0);
        d.setUserOverrideProperty("k0", 1);

        double[] lon = new double[] {1.0 * SnyderMath.DEG_TO_RAD}; // 1.0 degree
        double[] lat = new double[] {2.0 * SnyderMath.DEG_TO_RAD}; // 2.0 degree

        //
        // [0][] -> x[]
        // [1][] -> y[]
        //
        double[][] forward = proj.forward(lon, lat, e, d);
        System.out.println("x: " + forward[0][0] + ", y: " + forward[1][0]);

        //
        // [0][] -> lon[]
        // [1][] -> lat[]
        //
        double[][] inverse = proj.inverse(forward[0], forward[1], e, d);
        System.out.println(
                "lon: "   + (float) (inverse[0][0] * SnyderMath.RAD_TO_DEG) +
                ", lat: " + (float) (inverse[1][0] * SnyderMath.RAD_TO_DEG));
    }
```

Test Drivers + Deployment
=========================

The **src** / **test** / **resources** directory contains the test data fiducials,
projection and ellipsoid configuration parameter files necessary to run only the test drivers.

Fiducial values file: **src** / **test** / **resources** / **json** / **Fiducials.json**

Building Snyder
=========================

```
$ gradle build
```

The resulting directory **build** contains the built artifacts including the deployment jar file:<br/>
**build** / **libs** / **ap.snyder-{VERSION}.jar**

References
==========

   Snyder, J. P. Map Projections--A Working Manual.<br/>
   U. S. Geological Survey Professional Paper 1395.<br/>
   Washington, DC: U. S. Government Printing Office, 1987, 1994 3rd Printing.<br/>
   http://pubs.er.usgs.gov/publication/pp1395

   A Global Self-consistent, Hierarchical, High-resolution Geography Database<sup>\*</sup>
   http://www.soest.hawaii.edu/pwessel/gshhg<br/>
   Extracted from version 2.3.4 on 1 January 2015

License
=======

Code licensed under the MIT License

<sub><sup>\* SnyderVis, GSHHG files will no longer be updated.</sup></sub>
