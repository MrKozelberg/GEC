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

#pragma once

void gdz_to_mag(double dyear,
                double lati, double longi, double alti,
                double &latm, double &longm, double &altm);

double geomag_lat_by_gzd(double dyear, double lat, double lon, double alt=0.0);

void dtd(double x, double y, double z, double &Bx, double &By, double &Bz);

void init_dtd(double year);
