/*
 *  Copyright 2011 Brad Parks
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.bradsbrain.simpleastronomy;

import static com.bradsbrain.simpleastronomy.BaseUtils.adjustTo360Range;
import static com.bradsbrain.simpleastronomy.BaseUtils.cosDegrees;
import static com.bradsbrain.simpleastronomy.BaseUtils.sinDegrees;
import static com.bradsbrain.simpleastronomy.BaseUtils.tanDegrees;

import java.util.Calendar;

import android.util.Log;

public class SunPosition {
	// some handy constants
	private double EPOCH = 2447891.5; // 1990 January 0.0
	private double ECLIPTIC_LONGITUDE_OF_PERIGREE = 282.768422;  //(w)
	private double ECLIPTIC_LONGITUDE_AT_EPOCH_1990 = 279.403303;  
	private double ECCENTRICITY_OF_ORBIT = 0.016713;  //(e)

	static boolean DEBUG = true;
	private double sunElevation;
	private double sunAzimuth;
	
	/**
	 * The geocentric ecliptic longitude.  <br>
	 * Calculation is good to 3 decimal places <br>
	 * me: 337.44406603442917,   book: 337.444194
	 */
	private double geoEclipticLongitude = 0;   // oft represented as a lambda with little circle+dot
	/**
	 * The mean anomaly
	 */
	private double meanAnomaly = 0; // oft represented as capital M with little circle+dot

	public SunPosition(Calendar cal) {
		Calendar myCal = BaseUtils.getSafeLocalCopy(cal.getTimeInMillis());

		double daysSince = BaseUtils.exactDaysSince(myCal, EPOCH);

		double N = (360 / 365.242191 * daysSince) % 360;
		if (N < 0) {
			N += 360;
		}

		meanAnomaly = computeMeanAnomaly(N);
		geoEclipticLongitude = computeGeoEclipticLongitude(N);
	}


	public void calcalateSunAzEl(Calendar cal) {

		//Conversion of date date and time 

		//get Julian Day-
		//double JD = JulianDate.makeJulianDateUsingMyModified(cal);
		double JD = JulianDate.getJulianDay(cal);
		//convert to Julian Ephemeris Day
		double JDE = JD;
		//find the time in Julian Centuries
		double T = (JDE - 2451545.0) / 36525.0 ;

		if(DEBUG) Log.d("SP","Time in Julian centuries is " +T);


		//calculate the angle terms (degrees 0-360)


		/* Solar Coordinates (according to: Jean Meeus: Astronomical Algorithms), accuracy of 0.01 degree 
		 * 
		 */
		
		// mean geometric longitude of the sun, degree
		// aka mean equinox of date
		double L0 = 280.46645 + 36000.76983*T + 0.0003032*T*T; 
		
		//mean anomaly of the sun, degree
		double M = 357.52910 + 35999.05030*T - 0.0001559*T*T - 0.00000048*T*T*T ;

		//eccentricity of earth's orbit
		double e = 0.016798617 - 0.000042037*T - 0.0000001236*T*T;
		
		//the sun's equation of center
		double c = (1.914600 - 0.004817*T - 0.000014*T*T)*sinDegrees(M)
				+ (0.019993 - 0.000101*T)*sinDegrees(2*M) 
				+ 0.000290*sinDegrees(3*M);

		// sun's true longitude, degree
		double L = L0 + c ;


		/* convert ecliptic longitude L to right ascension RA and declination delta
		 *  	 (the ecliptic latitude of the Sun is assumed to be zero):
		 */
		// number of Julian centuries since Jan 1, 2000, 12 UT
		// T = (JD-2451545.0) / 36525

		// obliquity eps of ecliptic:
		double eps = 23.0 + 26.0/60.0 + 21.448/3600.0 - (46.8150*T + 0.00059*T*T - 0.001813*T*T*T)/3600.;

		double X = cosDegrees(L);
		double Y = cosDegrees(eps)*sinDegrees(L);
		double Z = sinDegrees(eps)*sinDegrees(L);
		double R = Math.sqrt(1.0-Z*Z);

		double delta = (180./Math.PI)*Math.atan2(Z,R); // in degrees
		double RA = (24./180.)*(180./Math.PI)*Math.atan2(Y,(X+R)); // in hours

		
		/* compute sidereal time at Greenwich (according to: Jean Meeus: Astronomical Algorithms)
		 */

	    //T = (JD - 2451545.0 ) / 36525

	    double theta0 = 280.46061837 + 360.98564736629*(JD-2451545.0) + 0.000387933*T*T - T*T*T/38710000.0;
	    
	    double beta = 30.77; // local geo. latitude
	    double mlong = -97.0;// local geo. longitude
	    
	    double theta = theta0 + mlong;
	    double tau = theta - RA;
	    
		//local hour angle
		//tau = adjustTo360Range(tau);  
		
	    /* convert tau, delta to horizon coordinates 
	     * of the observer (altitude h, azimuth az):
	     */

	    double h = (180./Math.PI)*Math.asin(sinDegrees(beta)*sinDegrees(delta)+cosDegrees(beta)*cosDegrees(delta)*cosDegrees(tau) );
	    double az = (180./Math.PI)*Math.atan2(-sinDegrees(tau),(cosDegrees(beta)*tanDegrees(delta)-sinDegrees(beta)*cosDegrees(tau)) );

		if(DEBUG) Log.d("SP","Sun Elevation: " + h);
		if(DEBUG) Log.d("SP","Sun Azimuth: " + az);

		//stash final values
		sunElevation = h;
		sunAzimuth = az;
	}

	public double getSunElevation() {
		   return sunElevation;
	}
	public double getSunAzimuth() {
		   return sunAzimuth;
	}
	
	private double computeGeoEclipticLongitude(double nValue) {
		double Ec = (360.0 / Math.PI) * ECCENTRICITY_OF_ORBIT * Math.sin(Math.toRadians(meanAnomaly));
		double preliminaryLongitude = nValue + Ec + ECLIPTIC_LONGITUDE_AT_EPOCH_1990;
		if (preliminaryLongitude > 360) {
			preliminaryLongitude -= 360;
		}
		return preliminaryLongitude;
	}

	private double computeMeanAnomaly(double nValue) {
		double someMean = nValue + ECLIPTIC_LONGITUDE_AT_EPOCH_1990 - ECLIPTIC_LONGITUDE_OF_PERIGREE;
		return someMean < 0 ? someMean + 360 : someMean;

	}

	/**
	 * TODO: implement this someday
	 */
	public RightAscension getRightAscension() {
		return null;
	}

	/**
	 * TODO: implement this someday
	 */
	public Declination getDeclination() {
		return null;
	}

	public double getEclipticLongitude() {
		return geoEclipticLongitude;
	}

	public double getMeanAnomaly() {
		return meanAnomaly;
	}

}
