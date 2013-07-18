Snyder Projection Implementation

Copyright (C) 2012, 2013  Jeremy J. Gibbons

References:

    Snyder, J. P. Map Projections--A Working Manual.
    U. S. Geological Survey Professional Paper 1395.
    Washington, DC: U. S. Government Printing Office, 1987, 1994 3rd Printing.
    http://pubs.er.usgs.gov/publication/pp1395

SnyderVis data provided by:

    A Global Self-consistent, Hierarchical, High-resolution Geography Database
    http://www.soest.hawaii.edu/pwessel/gshhg
    Extracted from version 2.2.3 on 1 July 2013

REPO NOTES
==========

The lib/ + resources/ directories are only required by the classes in the following package:

    gis.proj.drivers.*

The driver classes are not required for deployment and should be used as a tutorial and usage reference
for the library.  Since the package is not required for deployment, it is excluded in the build.xml.
