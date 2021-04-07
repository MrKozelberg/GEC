/*
 *  Copyright (C) 2017 IAPRAS - All Rights Reserved
 *
 *  This file is part of the GEC calculator.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "geomagnetic.h"


extern "C"
{
    // doubles in Fortran transmitted by pointer inputs and outputs both
    void gdz_geo_(double *lati, double *longi, double *alti, double *xx, double *yy, double *zz);
    void init_dtd_(double *dyear);
    // double arrays in Fortran transmitted like in C
    void geo_mag_(double xGEO[3], double xMAG[3]);
    void car_sph_(double xMAG[3], double *altm, double *latm, double *longm);

    void dtd_(double *x, double *y, double *z, double *Bx, double *By, double *Bz);
}

void init_dtd(double year) {
    init_dtd_(&year);
}

/** geomagnetic coordinates by geodetic coordinates
 *
 * @param dyear real value of year (2018.5 for 30 june 2018)
 * @param lati latitude in degrees
 * @param longi geomagnetic longitude in degrees
 * @param alti altitude in km
 * @param latm  geomagnetic latitude in degrees
 * @param longm   geomagnetic longitude in degrees
 * @param altm  altitude in km
 */
void gdz_to_mag(double dyear,
                double lati, double longi, double alti,
                double &latm, double &longm, double &altm)
{
    double xGEO[3], xMAG[3];

    gdz_geo_(&lati, &longi, &alti, &xGEO[0], &xGEO[1], &xGEO[2]);
//    init_dtd_(&dyear);
    geo_mag_(xGEO, xMAG);

    car_sph_(xMAG, &altm, &latm, &longm);
}


/** geomagnetic latitude in degrees
 *
 * @param dyear real value of year (2018.5 for 30 june 2018)
 * @param lat   latitude    in degrees
 * @param lon   longitude   in degrees
 * @param alt   altitude    in km
 * @return geomagnetic latitude in degrees
 */
double geomag_lat_by_gzd(double dyear, double lat, double lon, double alt)
{
    double mlat, mlon, malt;
    gdz_to_mag(dyear,
               lat, lon, alt,
               mlat, mlon, malt);
    return mlat;
}


void dtd(double x, double y, double z, double &Bx, double &By, double &Bz) {
    dtd_(&x, &y, &z, &Bx, &By, &Bz);
}
