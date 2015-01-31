Snyder Map Projection Implementation
====================================

Copyright (c) 2012-2015  APIS Point, LLC

Description
===========

Snyder is a map projection library implementing Snyder's seminal USGS projection
text including all fiducial appendix values for testing.  The code is written in
Java and runs on Android.  The core library code has no external dependencies.

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

The **lib** + **resources** directories contain the libraries, projection and
ellipsoid configuration parameter files necessary to run only the test drivers.
When creating a distribution jar file (e.g. snyder-*.jar) from the build
system, it is not necessary to include anything from either of the above
directories.

Fiducial values file: **src** / **test** / **resources** / **json** / **Fiducials.json**

Building Snyder
=========================

```
$ gradle build
```

The result is the directory **build** containing the build artifacts including the deployment jar.

Deployment jar library: **build** / **libs** / **snyder-{VERSION}.jar**

References
==========

    Snyder, J. P. Map Projections--A Working Manual.
    U. S. Geological Survey Professional Paper 1395.
    Washington, DC: U. S. Government Printing Office, 1987, 1994 3rd Printing.
    http://pubs.er.usgs.gov/publication/pp1395

    SnyderVis data++ provided by:

    A Global Self-consistent, Hierarchical, High-resolution Geography Database
    http://www.soest.hawaii.edu/pwessel/gshhg
    Extracted from version 2.3.4 on 1 January 2015

**NOTES**:
> ++ SnyderVis GSHHG data files will no longer be updated. For the latest GSHHG binary updates please follow the GSHHG URI reference above


Code licensed under the The MIT License
