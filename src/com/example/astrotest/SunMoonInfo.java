package com.example.astrotest;

import android.util.Log;

public class SunMoonInfo {



	static double pi = Math.PI; 
	static double r2d = (180./Math.PI);
	static double d2r  = (Math.PI/180.);

	static double DEFLAT = 9999;
	static double DEFLON = 9999;

	/* function declarations */
/*
	int timescale(int y, int m, int d);
	double eccentric_anomaly(double M, double e);
	void sunposition(double d, double* Ms, double* ws, double* RAs, double* Decs);
	void moonposition(double d,
						double Ms,
						double ws,
						double* rm,
						double* RAm,
						double* Decm);
	double topocentric_correction(double r, double alt_geoc, double lat);
	void altitude(double UT,
					double lon,
					double lat,
					double RA,
					double Decl,
					double Ms,
					double ws,
					double* AZ,
					double* Alt);
	double modpi(double a);
	static void usage();

*/

	/*-----------------------------------------------------------------------
	Main orbital elements
	    N = longitude of the ascending node
	    i = inclination to the ecliptic (plane of the Earth's orbit)
	    w = argument of perihelion
	    a = semi-major axis, or mean distance from Sun
	    e = eccentricity (0=circle, 0-1=ellipse, 1=parabola)
	    M = mean anomaly (0 at perihelion; increases uniformly with time)

	Related orbital elements
	    w1 = N + w = longitude of perihelion
	    L  = M + w1 = mean longitude
	    q  = a*(1-e) = perihelion distance
	    Q  = a*(1+e) = aphelion distance
	    P  = a ^ 1.5 = orbital period (years if a is in AU, astronomical units)
	    T  = Epoch_of_M - (M(deg)/360) / P  = time of perihelion
	    v  = true anomaly (angle between position and perihelion)
	    E  = eccentric anomaly
	--------------------------------------------------------------------------*/

	int main(int argc, char* argv[])
	{
		int	i;

		double Ms;			/* Mean Anomaly of the Sun */
		double ws;			/* Argument of perihelion for the Sun */
		double RAs;
		double Decs;
		double RAm;
		double Decm;
		double rm;			/* Moon distance in Earth radii */

		double lon;
		double lat;
		int year;
		int month;
		int day;

		char*	datearg;

		if (argc != 4) {
			usage();
			return 1;
		}

		lon = DEFLON;
		lat = DEFLAT;
		datearg = 0;

		for (i = 1; i < 4; ++i) {
			if (strncmp(argv[i], "lat=", 4) == 0) {
				lat = atof(argv[i] + 4);
			} else if (strncmp(argv[i], "lon=", 4) == 0) {
				lon = atof(argv[i] + 4);
			} else if (strncmp(argv[i], "date=", 5) == 0) {
				datearg = argv[i] + 5;
			}
		}

		if (lon == DEFLON || lat == DEFLAT || datearg == 0) {
			usage();
			return 1;
		}

		year = 2000;
		month = 1;
		day = 1;

		switch (sscanf(datearg, "%d-%d-%d", &year, &month, &day)) {
		case 3:
			break;
		case 2:
			day = 1;
			break;
		case 1:
			month = 1;
			day = 1;
			break;
		case 0:
		default:
			usage();
			return 1;
		}

		if (month < 1 || month > 12) {
			fprintf(stderr, "bad month value\n");
			usage();
			return 1;
		}

		/*---
		note: allow arbitrary day, including forays into next or previous months

		if (day < 1 || day > 31) {
			fprintf(stderr, "bad day value\n");
			usage();
			return 1;
		}
		----*/

		printf("longitude %f\n", lon);
		printf("latitude %f\n", lat);
		printf("yyyy-mm-dd %04d-%02d-%02d\n", year, month, day);
		printf("\n");

		lon *= d2r;
		lat *= d2r;
		double d = timescale(year, month, day);


		printf("            Sun            Moon\n");
		printf("UTC      Alt      AZ    Alt      AZ\n");

		for (int i = 0; i < 24*6; ++i) {
			double UT = 0 + ((double)i)/(24*6);
			sunposition(d + UT, &Ms, &ws, &RAs, &Decs);
			moonposition(d + UT, Ms, ws, &rm, &RAm, &Decm);
			double AZs, Alts;
			double AZm, Altm;
			altitude(UT*24, lon, lat, RAs, Decs, Ms, ws, &AZs, &Alts);
			altitude(UT*24, lon, lat, RAm, Decm, Ms, ws, &AZm, &Altm);
			Altm = topocentric_correction(rm, Altm, lat);
			AZs = modpi(AZs);
			AZm = modpi(AZm);

			/* express the day into the correct month.
				if the day is already in the right month, this code does nothing */
			struct tm stm;
			memset(&stm, 0, sizeof(stm));
			stm.tm_year = year - 1900;
			stm.tm_mon = month - 1;
			stm.tm_mday = day;
			time_t mkt = mktime(&stm);
			struct tm* gmt = gmtime(&mkt);

			printf("%04d-%02d-%02d ", gmt->tm_year + 1900, gmt->tm_mon + 1, gmt->tm_mday);
			printf("%02d:%02d  %5.1f  %6.1f  %5.1f  %6.1f\n",
					i/6,
					10*(i%6),
					Alts*r2d,
					AZs*r2d,
					Altm*r2d,
					AZm*r2d);

		}
	}

	static void usage()
	{
//		fprintf(stderr, "usage: altaz lat={lat} lon={lon} date={yyyy-mm-dd}\n");
		Log.d("INFO","usage: altaz lat={lat} lon={lon} date={yyyy-mm-dd}\n");

	}

	int timescale(int y, int m, int d)
	{
		/*-----------
		The time scale used here is a "day number" from 2000 Jan 0.0 TDT, which
		is the same as 1999 Dec 31.0 TDT, i.e. precisely at midnight TDT at the
		start of the last day of this century.  With the modest accuracy we
		strive for here, one can usually disregard the difference between
		TDT (formerly canned ET) and UT.
	    d  =  JD - 2451543.5  =  MJD - 51543.0
		We can also compute d directly from the calendar date like this:
	    d = 367*Y - (7*(Y + ((M+9)/12)))/4 + (275*M)/9 + D - 730530
		------------*/

		/* days since jan 01 2000 */
		return 367*y - 7*(y + (m+9)/12)/4 + 275*m/9 + d - 730530;
	}

	double eccentric_anomaly(double M, double e)
	{ 
		/* E and M in radians */
	    double E0 = M + e * sin(M) * ( 1.0 + e * cos(M) );
	    double E1 = E0 - ( E0 - e * sin(E0) - M ) / ( 1 - e * cos(E0) );
		int iteration = 0;
		if (e > 0.05) {
			while (fabs(E1 - E0) > 0.001*d2r) {
	    		E0 = E1;
	    		E1 = E0 - ( E0 - e * sin(E0) - M ) / ( 1 - e * cos(E0) );
				if (++iteration >= 10) {
//					fprintf(stderr, "eccentric anomaly does not converge\n");
//					exit(1);
					Log.e("ERR","eccentric anomaly does not converge");
					break;
				}
			}
		}
		return E1;
	}

	void sunposition(double d, double Ms, double ws, double RAs, double Decs)
	{
		/* obliquity of the ecliptic */

	    double ecl = 23.4393 - 3.563E-7 * d;
		ecl *= d2r;

		/* Sun orbital elements (angles in degrees)
			N = longitude of the ascending node
	    	i = inclination to the ecliptic (plane of the Earth's orbit)
	    	w = argument of perihelion
	    	a = semi-major axis, or mean distance from Earth
	    	e = eccentricity (0=circle, 0-1=ellipse, 1=parabola)
	    	M = mean anomaly (0 at perigee; increases uniformly with time)
		*/

	    double N = 0.0;
	    double i = 0.0;
	    double w = 282.9404 + 4.70935E-5 * d;
	    double a = 1.000000;		  /* AU */
	    double e = 0.016709 - 1.151E-9 * d;
	    double M = 356.0470 + 0.9856002585 * d;

		N *= d2r;
		i *= d2r;
		w *= d2r;
		M *= d2r;

		/* Excentric anomaly from the mean anomaly */

	    double E = M + e*sin(M)*(1.0 + e*cos(M));

		/*
			Compute Sun's distance r and its true anomaly v using
	    		xv = r * cos(v) = a * (cos(E) - e)
	    		yv = r * sin(v) = a * sqrt(1.0 - e*e) * sin(E)
		*/

	    double xv = cos(E) - e;
	    double yv = sqrt(1.0 - e*e) * sin(E);
	    double v = atan2(yv, xv);
	    double r = sqrt(xv*xv + yv*yv);

		/* the Sun's true longitude: */

	    double lonsun = v + w;

		/* Convert lonsun,r to ecliptic rectangular geocentric coordinates
		   xs,ys: */

	    double xs = r*cos(lonsun);
	    double ys = r*sin(lonsun);
		/* zs is zero */

		/* convert this to equatorial, rectangular, geocentric coordinates: */

	    double xe = xs;
	    double ye = ys*cos(ecl);
	    double ze = ys*sin(ecl);

		/* compute the Sun's Right Ascension (RA) and Declination (Dec): */

	    double RA  = atan2(ye, xe);
	    double Dec = atan2(ze, sqrt(xe*xe+ye*ye));

	/*--
		printf("Sun's RA %g\n", RA*r2d);
		printf("Sun's Dec %g\n", Dec*r2d);
	--*/

		Ms = M;
		ws = w;
		RAs = RA;
		Decs = Dec;
	}

	void moonposition(double d,
						double Ms,
						double ws,
						double rm,
						double RAm,
						double Decm)
	{
		/* I. The first order calculation is similar for all planets */

		/* obliquity of the ecliptic */

	    double ecl = 23.4393 - 3.563E-7 * d;
		ecl *= d2r;

		/* Moon orbital elements (angles in degrees)
			N = longitude of the ascending node
	    	i = inclination to the ecliptic (plane of the Earth's orbit)
	    	w = argument of perihelion
	    	a = semi-major axis, or mean distance from Earth
	    	e = eccentricity (0=circle, 0-1=ellipse, 1=parabola)
	    	M = mean anomaly (0 at perigee; increases uniformly with time)
		*/

	    double N = 125.1228 - 0.0529538083 * d;
	    double i = 5.1454;
	    double w = 318.0634 + 0.1643573223 * d;
	    double a = 60.2666;				  /* Earth radii */
	    double e = 0.054900;
	    double M = 115.3654 + 13.0649929509 * d;

		N *= d2r;
		i *= d2r;
		w *= d2r;
		M *= d2r;

		/*
			First, compute the eccentric anomaly, E, from M, the mean anomaly,
			and e, the eccentricity.
		*/

		double E = eccentric_anomaly(M, e);

		/*
			Now compute the Moon's distance and true anomaly, using
	    	xv = r * cos(v) = a * ( cos(E) - e )
	    	yv = r * sin(v) = a * ( sqrt(1.0 - e*e) * sin(E) )
		*/
	    double xv = a * ( cos(E) - e );
	    double yv = a * ( sqrt(1.0 - e*e) * sin(E) );
	    double v = atan2( yv, xv );
	    double r = sqrt( xv*xv + yv*yv );

		/* Compute the Moon's position in 3-dimensional space: */

	    double xg = r * ( cos(N) * cos(v+w) - sin(N) * sin(v+w) * cos(i) );
	    double yg = r * ( sin(N) * cos(v+w) + cos(N) * sin(v+w) * cos(i) );
	    double zg = r * ( sin(v+w) * sin(i) );

		/*
			As a check one can compute sqrt(xh*xh+yh*yh+zh*zh), which of course
			should equal r (except for small round-off errors).
		*/
	/*--
			printf("\tcheck 1: r %f = sqrt(xg*xg+yg*yg+zg*zg) %f\n",
						r,
						sqrt(xg*xg+yg*yg+zg*zg));
	--*/

		/*
			This is the geocentric (Earth-centered) position in the ecliptic
			coordinate system.
			If one wishes, one can compute the ecliptic longitude and
			latitude (this must be done if one wishes to correct for
			perturbations, or if one wants to precess the position to a standard
			epoch.
		*/

		/* the ecliptic longitude and latitude: */

	    double lonecl = atan2( yg, xg );
	    double latecl = atan2( zg, sqrt(xg*xg + yg*yg) );

		/* II. Perturbations of the Moon */

		/*
			If the position of the Moon is computed, and one wishes a better
			accuracy than about 2 degrees, the most important perturbations has
			to be taken into account.  If one wishes 2 arc minute accuracy, all
			the following terms should be accounted for.  If less accuracy is
			needed, some of the smaller terms can be omitted.
		*/

		/* First compute:

	    	Ms, Mm             Mean Anomaly of the Sun and the Moon
	    	Nm                 Longitude of the Moon's node
	    	ws, wm             Argument of perihelion for the Sun and the Moon
		*/

		double Mm = M;
		double Nm = N;
		double wm = w;

	    double Ls = Ms + ws;       /* Mean Longitude of the Sun  (Ns=0) */
	    double Lm = Mm + wm + Nm;  /* Mean longitude of the Moon */
	    double D = Lm - Ls;        /* Mean elongation of the Moon */
	    double F = Lm - Nm;        /* Argument of latitude for the Moon */

		/* Add these terms to the Moon's longitude (degrees): */

		lonecl += (
	    	-1.274 * sin(Mm - 2*D)          /* the Evection */
	    	+0.658 * sin(2*D)               /* the Variation */
	    	-0.186 * sin(Ms)                /* the Yearly Equation */
	    	-0.059 * sin(2*Mm - 2*D)
	    	-0.057 * sin(Mm - 2*D + Ms)
	    	+0.053 * sin(Mm + 2*D)
	    	+0.046 * sin(2*D - Ms)
	    	+0.041 * sin(Mm - Ms)
	    	-0.035 * sin(D)                 /* the Parallactic Equation */
	    	-0.031 * sin(Mm + Ms)
	    	-0.015 * sin(2*F - 2*D)
	    	+0.011 * sin(Mm - 4*D))*d2r;

		/* Add these terms to the Moon's latitude (degrees): */

		latecl += (
	    	-0.173 * sin(F - 2*D)
	    	-0.055 * sin(Mm - F - 2*D)
	    	-0.046 * sin(Mm + F - 2*D)
	    	+0.033 * sin(F + 2*D)
	    	+0.017 * sin(2*Mm + F))*d2r;

		/* Add these terms to the Moon's distance (Earth radii): */

		r +=
	    	-0.58 * cos(Mm - 2*D)
	    	-0.46 * cos(2*D);

		/*
			All perturbation terms that are smaller than 0.01 degrees in
			longitude or latitude and smaller than 0.1 Earth radii in distance
			have been omitted here.  A few of the largest perturbation terms even
			have their own names!  The Evection (the largest perturbation) was
			discovered already by Ptolemy a few thousand years ago (the Evection
			was one of Ptolemy's epicycles). The Variation and the Yearly
			Equation were both discovered by Tycho Brahe in the 16'th
			century.
		*/

		/*
			Now we have computed the geocentric (Earth-centered) coordinate of
			the Moon, and we have included the most important perturbations.
			We should convert the perturbed lonecl, latecl, r to (perturbed)
			xg, yg, zg:
		*/

	    xg = r * cos(lonecl) * cos(latecl);
	    yg = r * sin(lonecl) * cos(latecl);
	    zg = r               * sin(latecl);

		/*
			We now have the Moons's geocentric (Earth centered) position in
			rectangular, ecliptic coordinates.
		*/

		/* Equatorial coordinates */

		/*
			Let's convert our rectangular, ecliptic coordinates to rectangular,
			equatorial coordinates: simply rotate the y-z-plane by ecl, the angle
			of the obliquity of the ecliptic:
		*/

	    double xe = xg;
	    double ye = yg * cos(ecl) - zg * sin(ecl);
	    double ze = yg * sin(ecl) + zg * cos(ecl);

		/*
			Finally, compute the planet's Right Ascension (RA)
			and Declination (Dec):
		*/

	    double RA  = atan2( ye, xe );
	    double Dec = atan2( ze, sqrt(xe*xe+ye*ye) );

	/*--
		printf("Moon's RA %g\n", RA*r2d);
		printf("Moon's Dec %g\n", Dec*r2d);
	--*/

		/* Compute the geocentric distance:
			rg = sqrt(xg*xg+yg*yg+zg*zg) = sqrt(xe*xe+ye*ye+ze*ze)
		*/

	/*--
		printf("\tcheck 2: r %f = sqrt(xg*xg+yg*yg+zg*zg) %f\n",
					r,
					sqrt(xg*xg+yg*yg+zg*zg));
		printf("\tcheck 3: r %f = sqrt(xe*xe+ye*ye+ze*ze) %f\n",
					r,
					sqrt(xe*xe+ye*ye+ze*ze));
	--*/

		/* Topocentric coordinates */
		/* I ignore this small effect */

		rm = r;
		RAm = RA;
		Decm = Dec;
	}

	double topocentric_correction(double r, double alt_geoc, double lat)
	{
		/* The Moon's topocentric position */

		/*
			The Moon's position, as computed earlier, is geocentric, i.e. as seen
			by an imaginary observer at the center of the Earth.  Real observers
			dwell on the surface of the Earth, though, and they will see a
			different position - the topocentric position.  This position can
			differ by more than one degree from the geocentric position.  To
			compute the topocentric positions, we must add a correction to the
			geocentric position.
		*/

		/*
			Let's start by computing the Moon's parallax, i.e. the apparent
			size of the (equatorial) radius of the Earth, as seen from the Moon:
		*/

	    double mpar = asin( 1/r );

		/*
			where r is the Moon's distance in Earth radii.  It's simplest to apply
			the correction in horizontal coordinates (azimuth and altitude):
			within our accuracy aim of 1-2 arc minutes, no correction need to be
			applied to the azimuth.  One need only apply a correction to the
			altitude above the horizon:
		*/

	    double alt_topoc = alt_geoc - mpar * cos(alt_geoc);
		return alt_topoc;
	}

	void altitude(double UT,
				double lon,
				double lat,
				double RA,
				double Decl,
				double Ms,
				double ws,
				double AZ,
				double Alt)
	{
		/*
			First we compute the Sun's RA and Decl for this moment.  Now we need
			to know our Local Sidereal Time.  We start by computing the sidereal
			time at Greenwich at 00:00 Universal Time, let's call this quantity
			GMST0:
		*/

	    double L = Ms + ws;		/* the Sun's mean longitude */
	    double GMST0 = L*r2d + 180;

		/*
			Note that we express GMST in degrees here to simplify the
			computations.  360 degrees of course corresponds to 24 hours, i.e.
			each hour corresponds to 15 degrees.
		*/

		/*
			Now we can compute our Local Sidereal Time (LST):
		*/

	    double LST = GMST0 + UT*15.0 + lon*r2d;

		/*
			UT is the Universal Time, expressed in hours+decimals, the remaining
			quantities are expressed in degrees.  To convert UT to degrees we
			must multiply it by 15 above.  long is our local longitude in
			degrees, where east longitude counts as positive, and west longitude
			as negative.  (this is according to the geographic standard, and
			the recent astronomical standard; if you prefer to use the older
			astronomical standard where west longitude counts as positive, then
			you must change the '+' in front of 'lon' to a '-' above).
		*/

		/*
			Next let's compute the Sun's Local Hour Angle (LHA), i.e. the angle
			the Earth has turned since the Sun last was in the south:
		*/

	    double LHA = LST - RA*r2d;

		/*
			A negative hour angle means the Sun hasn't been in the south yet, this
			day.  The angle -10 degrees is of course the same as 350 degrees,
			i.e. adding or subtracting even multiples of 360 degrees does not
			change the angle.
		*/

		/*
			We also need to know our latitude (lat), where north latitude
			counts as positive and south latitude as negative.  Now we
			can compute the Sun's altitude above the horizon:
		*/

		/*
			We compute sin(h), and then take the arcsine of this to get h, the
			Sun's altitude above the horizon.
		*/

	    double sin_h = sin(lat) * sin(Decl) + cos(lat) * cos(Decl) * cos(LHA*d2r);
		double h = asin(sin_h);
		Alt = h;

		/* formula for the azimuth, reckoned eastward from north: */
		double azim = atan2(-sin(LHA*d2r)*cos(Decl),
							cos(lat)*sin(Decl) - sin(lat)*cos(Decl)*cos(LHA*d2r)); 
		AZ = azim;
	}

	double modpi(double a)
	{
		while (a < 0) {
			a += 2*pi;
		}
		while (a > 2*pi) {
			a -= 2*pi;
		}
		return a;
	}

	/*math functions*/
	double sin(double d) {
		return Math.sin(d);
	}
	double cos(double d) {
		return Math.cos(d);
	}
	double tan(double d) {
		return Math.tan(d);
	}
	double atan2(double n, double d) {
		return Math.atan2(n,d);
	}
	double asin(double d) {
		return Math.asin(d);
	}
	double sqrt(double d) {
		return Math.sqrt(d);
	}
	double fabs(double d) {
		return Math.abs(d);
	}
	
}
