package com.bradsbrain.simpleastronomy;

import java.util.Calendar;

public class MoonAngles {

//			% Programed by Darin C. Koblick 2/14/2009
//			%
//			% Updated on 03/04/2009 to clean up code and add quadrant check to Azimuth
//			% Thank you Doug W. for your help with the test code to find the quadrant check 
//			% error.
//			%
//			% Updated on 04/13/2009 to add Lunar perturbation offsets (this will
//			% increase the accuracy of calculations)
//			%
//			% Updated on 08/17/2010 to make use of the site altitude (this will affect
//			% the elevation angle)
//
//			% External Function Call Sequence:
//			% [Az El] = LunarAzEl('1991/05/19 13:00:00',50,10,0)
//
//
//			% Function Description:
//			% LunarAzEl will ingest a Universal Time, and specific site location on earth
//			% it will then output the lunar Azimuth and Elevation angles relative to that
//			% site.
//
//			% External Source References:
//
//			% Basics of Positional Astronomy and Ephemerides
//			% http://jgiesen.de/elevazmoon/basics/index.htm
//
//			% Computing planetary positions - a tutorial with worked examples
//			% http://stjarnhimlen.se/comp/tutorial.html
//
//			%Input Description:
//			% UTC (Coordinated Universal Time YYYY/MM/DD hh:mm:ss)
//			% Lat (Site Latitude in degrees -90:90 -> S(-) N(+))
//			% Lon (Site Longitude in degrees -180:180 W(-) E(+))
//			% Altitude of the site above sea level (km)
//
//			%Output Description:
//			%Az (Azimuth location of the moon in degrees)
//			%El (Elevation location of the moon in degrees)
//
//			%Verified output by comparison with the following source data:
//			%http://aa.usno.navy.mil/data/docs/AltAz.php

			// Code Sequence
			//--------------------------------------------------------------------------
			//Do initial Longitude Latitude check

		public void calculateAngles(double Lat, double Lon, Calendar cal) {
	
			if( Lon > 180.)
			   Lon = Lon - 360.;
			else if( Lon < -180.)
			   Lon = Lon + 360.; 
			if (Lat > 90.)
			    Lat = Lat - 360.;
			else if (Lat < -90.)
			   Lat = Lat + 360.; 


			//Declare Earth Equatorial Radius Measurements in km
			double EarthRadEq = 6378.1370;

			//Convert Universal Time to Ephemeris Time
			//jd = juliandate(UTC,'yyyy/mm/dd HH:MM:SS');

			double jd = JulianDate.getJulianDay(cal);
			
			
			//Find the Day Number
			double d = jd - 2451543.5;

			//Keplerian Elements of the Moon
			//This will also account for the Sun's perturbation
			 double N = 125.1228 - 0.0529538083 * d; 	//   (Long asc. node deg)
			 double i = 5.1454; 						//    (Inclination deg)
			 double w = 318.0634 + 0.1643573223 * d; 	//    (Arg. of perigee deg)
			 double a =  60.2666;						//    (Mean distance (Earth's Equitorial Radii)
			 double e = 0.054900;						//    (Eccentricity)
			 double M = (115.3654+13.0649929509*d)%(360);//(Mean anomaly deg)
			  
			 double LMoon =  (N + w + M)%(360);        //(Moon's mean longitude deg)
			 double FMoon =  (LMoon - N)%(360);        //(Moon's argument of latitude)

			 //Keplerian Elements of the Sun
			 double wSun = (282.9404 + 4.70935E-5*d)%(360);    // (longitude of perihelion)
			 double MSun = (356.0470 + 0.9856002585*d)%(360);  // (Sun mean anomaly)
			 double LSun = (wSun + MSun) % 360;                 // (Sun's mean longitude)
			     
			 double DMoon =  LMoon - LSun;                     // (Moon's mean elongation)  

			 //Calculate Lunar perturbations in Longitude
			 LunarPLon = [ -1.274.*Math.sin((M - 2.*DMoon).*(Math.PI/180)); ...
			     .658.*Math.sin(2.*DMoon.*(Math.PI/180)); ...
			     -0.186.*Math.sin(MSun.*(Math.PI/180)); ...
			     -0.059.*Math.sin((2.*M-2.*DMoon).*(Math.PI/180)); ...
			     -0.057.*Math.sin((M-2.*DMoon + MSun).*(Math.PI/180)); ...
			     .053.*Math.sin((M+2.*DMoon).*(Math.PI/180)); ...
			     .046.*Math.sin((2.*DMoon-MSun).*(Math.PI/180)); ...
			     .041.*Math.sin((M-MSun).*(Math.PI/180)); ...
			    -0.035.*Math.sin(DMoon.*(Math.PI/180)); ...           
			    -0.031.*Math.sin((M+MSun).*(Math.PI/180)); ...
			    -0.015.*Math.sin((2.*FMoon-2.*DMoon).*(Math.PI/180)); ...
			    .011.*Math.sin((M-4.*DMoon).*(Math.PI/180))];
			 
			 //Calculate Lunar perturbations in Latitude 
			 LunarPLat = [ -0.173.*Math.sin((FMoon-2.*DMoon).*(Math.PI/180)); ...
			    -0.055.*Math.sin((M-FMoon-2.*DMoon).*(Math.PI/180)); ...
			    -0.046.*Math.sin((M+FMoon-2.*DMoon).*(Math.PI/180)); ...
			    +0.033.*Math.sin((FMoon+2.*DMoon).*(Math.PI/180)); ...
			    +0.017.*Math.sin((2.*M+FMoon).*(Math.PI/180))];

			//Calculate perturbations in Distance
			 LunarPDist = [ -0.58*Math.cos((M-2.*DMoon).*(Math.PI/180)); ...
			    -0.46.*Math.cos(2.*DMoon.*(Math.PI/180))];

			// Compute E, the eccentric anomaly

			//E0 is the eccentric anomaly approximation estimate 
			//(this will initially have a relatively high error)
			double E0 = M+(180/Math.PI)*e*Math.sin(M*(Math.PI/180))*(1+e*Math.Math.cos(M*(Math.PI/180)));

			//Compute E1 and set it to E0 until the E1 == E0
			double E1 = E0-(E0-(180/Math.PI)*e*Math.sin(E0*(Math.PI/180))-M)/(1-e*Math.cos(E0*(Math.PI/180)));
			if ( (E1-E0) > .000005 ){
			    E0 = E1;
			    E1 = E0-(E0-(180/Math.PI)*e*Math.sin(E0*(Math.PI/180))-M)/(1-e*Math.cos(E0*(Math.PI/180)));    
			}
			double E = E1;

			//Compute rectangular coordinates (x,y) in the plane of the lunar orbit
			double x = a*(Math.cos(E*(Math.PI/180))-e);
			double y = a*Math.sqrt(1-e*e)*Math.sin(E*(Math.PI/180));

			//convert this to distance and true anomaly
			double r = Math.sqrt(x*x + y*y);
			double v = Math.atan2(y*(Math.PI/180),x*(Math.PI/180))*(180/Math.PI);

			//Compute moon's position in ecliptic coordinates
			double xeclip = r*(Math.cos(N*(Math.PI/180))*Math.cos((v+w)*(Math.PI/180))-Math.sin(N*(Math.PI/180))*Math.sin((v+w)*(Math.PI/180))*Math.cos(i*(Math.PI/180)));
			double yeclip = r*(Math.sin(N*(Math.PI/180))*Math.cos((v+w)*(Math.PI/180))+Math.cos(N*(Math.PI/180))*Math.sin(((v+w)*(Math.PI/180)))*Math.cos(i*(Math.PI/180)));
			double zeclip = r*Math.sin((v+w)*(Math.PI/180))*Math.sin(i*(Math.PI/180));

			//Add the calculated lunar perturbation terms to increase model fidelity
			[eLon eLat eDist] = cart2sph(xeclip,yeclip,zeclip);
			[xeclip yeclip zeclip] = sph2cart(eLon + sum(LunarPLon).*(Math.PI/180), ...
			                                  eLat + sum(LunarPLat).*(Math.PI/180), ...
			                                  eDist + sum(LunarPDist));
			clear eLon eLat eDist;
			                              
			//convert the latitude and longitude to right ascension RA and declination
			//delta
			T = (jd-2451545.0)/36525.0;

			//Generate a rotation matrix for ecliptic to equitorial
			//RotM=rotm_coo('E',jd);
			//See rotm_coo.m for obl and rotational matrix transformation
			Obl = 23.439291 - 0.0130042.*T - 0.00000016.*T.*T + 0.000000504.*T.*T.*T;
			Obl = Obl.*(Math.PI/180);
			RotM = [1 0 0; 0 Math.cos(Obl) Math.sin(Obl); 0 -Math.sin(Obl) Math.cos(Obl)]';

			//Apply the rotational matrix to the ecliptic rectangular coordinates
			//Also, convert units to km instead of earth equatorial radii
			sol = RotM*[xeclip yeclip zeclip]'.*EarthRadEq;

			//Find the equatorial rectangular coordinates of the location specified
			[xel yel zel] = sph2cart(Lon.*(Math.PI/180),Lat.*(Math.PI/180),Alt+EarthRadEq);

			//Find the equatorial rectangular coordinates of the location @ sea level
			[xsl ysl zsl] = sph2cart(Lon.*(Math.PI/180),Lat.*(Math.PI/180),EarthRadEq);

			//Find the Angle Between sea level coordinate vector and the moon vector
			theta1 = 180 - aMath.cosd(dot([xsl ysl zsl],[sol(1)-xsl sol(2)-ysl sol(3)-zsl]) ...
			        ./(sqrt(xsl.^2 + ysl.^2 + zsl.^2) ...
			         .*sqrt((sol(1)-xsl).^2 + (sol(2)-ysl).^2 + (sol(3)-zsl).^2)));

			//Find the Angle Between the same coordinates but at the specified elevation
			theta2 = 180 - aMath.cosd(dot([xel yel zel],[sol(1)-xel sol(2)-yel sol(3)-zel]) ...
			    ./(sqrt(xel.^2 + yel.^2 + zel.^2) ...
			         .*sqrt((sol(1)-xel).^2 + (sol(2)-yel).^2 + (sol(3)-zel).^2)));
			     
			//Find the Difference Between the two angles (+|-) is important
			thetaDiff = theta2 - theta1;

			// equatorial to horizon coordinate transformation
			 [RA,delta] = cart2sph(sol(1),sol(2),sol(3));
			 delta = delta.*(180/Math.PI);
			 RA = RA.*(180/Math.PI);
			 
			//Following the RA DEC to Az Alt conversion sequence explained here:
			//http://www.stargazing.net/kepler/altaz.html

			//Find the J2000 value
			J2000 = jd - 2451545.0;
			hourvec = datevec(UTC,'yyyy/mm/dd HH:MM:SS');
			UTH = hourvec(4) + hourvec(5)/60 + hourvec(6)/3600;

			//Calculate local siderial time
			LST = mod(100.46+0.985647.*J2000+Lon+15*UTH,360);

			//Replace RA with hour angle HA
			HA = LST-RA;

			//Find the h and AZ at the current LST
			h = aMath.sin(Math.sin(delta.*(Math.PI/180)).*Math.sin(Lat.*(Math.PI/180)) + Math.cos(delta.*(Math.PI/180)).*Math.cos(Lat.*(Math.PI/180)).*Math.cos(HA.*(Math.PI/180))).*(180/Math.PI);
			Az = aMath.cos((Math.sin(delta.*(Math.PI/180)) - Math.sin(h.*(Math.PI/180)).*Math.sin(Lat.*(Math.PI/180)))./(Math.cos(h.*(Math.PI/180)).*Math.cos(Lat.*(Math.PI/180)))).*(180/Math.PI);

			//Add in the angle offset due to the specified site elevation
			h = h + thetaDiff;

			if Math.sin(HA.*(Math.PI/180)) >= 0
			   Az = 360-Az; 
			end

			//Apply Paralax Correction if we are still on earth
			if Alt < 100
			    horParal = 8.794/(r*6379.14/149.59787e6);
			    p = aMath.sin(Math.cos(h.*(Math.PI/180))*Math.sin((horParal/3600).*(Math.PI/180))).*(180/Math.PI);
			    h = h-p;
			end

			
//			//This sub function is provided in case juliandate does not come with your 
//			//distribution of Matlab
//
//			[year month day hour min sec] = datevec(datenum(varargin{:}));
//
//			for k = length(month):-1:1
//			    if ( month(k) <= 2 ) % january & february
//			        year(k)  = year(k) - 1.0;
//			        month(k) = month(k) + 12.0;
//			    end
//			end
//
//			jd = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
//			    floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5 + ...
//			    (hour + min/60 + sec/3600)/24;
			    
			    
		}
	
}
