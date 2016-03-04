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

import static com.bradsbrain.simpleastronomy.BaseUtils.sinDegrees;
import static com.bradsbrain.simpleastronomy.BaseUtils.cosDegrees;
import static com.bradsbrain.simpleastronomy.BaseUtils.tanDegrees;
import static com.bradsbrain.simpleastronomy.BaseUtils.adjustTo360Range;

import java.text.DecimalFormat;
import java.util.Calendar;

import android.util.Log;

public class MoonPosition {
    // some handy constants
    private double EPOCH = 2447891.5; // 1990 January 0.0
    private double MEAN_LONGITUDE_AT_EPOCH = 318.351648;
    private double MEAN_LONGITUDE_OF_PERIGREE_AT_EPOCH = 36.340410;
    
    private boolean DEBUG = true;
    
    double moonLambda;
    double moonBeta;
    double moonDelta;
    double moonElevationApparent;
    double moonAzimuth;
    
    /**
     * The True Longitude
     */
    private double trueOrbitalLongitude;

    /**
     * This is from section 65, page 144
     *
     * @param cal the calendar date for which to compute the moon position
     */
    public MoonPosition(Calendar cal) {
        Calendar myCal = BaseUtils.getSafeLocalCopy(cal.getTimeInMillis());
        double daysSince = BaseUtils.exactDaysSince(myCal, EPOCH);

        // l
        double moonMeanLongitude = computeMeanLongitude(daysSince);
        // M m
        double moonMeanAnomaly = computeMeanAnomaly(daysSince, moonMeanLongitude);

        SunPosition sunPos = new SunPosition(myCal);

        MoonCorrections corrections = new MoonCorrections(moonMeanLongitude, moonMeanAnomaly, sunPos.getEclipticLongitude(), sunPos.getMeanAnomaly());
        trueOrbitalLongitude = corrections.getCorrectedLongitude() - corrections.getVariationCorrection();
    }

    public double getTrueLongitude() {
        return trueOrbitalLongitude;
    }

    /**
     * Compute the Moon Mean Longitude	l
     *
     * @param daysSince
     * @return
     */
    private double computeMeanLongitude(double daysSince) {
        double moonMeanLongitude = 13.1763966 * daysSince + MEAN_LONGITUDE_AT_EPOCH;
        return BaseUtils.adjustTo360Range(moonMeanLongitude);
    }

    /**
     * Compute the Moon Mean Anomaly	M m
     *
     * @param daysSince
     * @param moonMeanLongitude
     * @return
     */
    private double computeMeanAnomaly(double daysSince, double moonMeanLongitude) {
        double moonMeanAnomaly = moonMeanLongitude
                - (0.1114041 * daysSince)
                - MEAN_LONGITUDE_OF_PERIGREE_AT_EPOCH;
        return BaseUtils.adjustTo360Range(moonMeanAnomaly);
    }

    /**
     * Private internal class used for lots of computations and corrections
     */
    private static class MoonCorrections {
        private double moonMeanLongitude;
        private double moonMeanAnomaly;
        private double sunLongitude;
        private double sunMeanAnomaly;

        MoonCorrections(double moonMeanLongitude, double moonMeanAnomaly,
                        double sunLongitude, double sunMeanAnomaly) {
            this.moonMeanAnomaly = moonMeanAnomaly;
            this.sunMeanAnomaly = sunMeanAnomaly;
            this.moonMeanLongitude = moonMeanLongitude;
            this.sunLongitude = sunLongitude;
        }

        /**
         * V
         *
         * @return
         */
        public double getVariationCorrection() {
            return 0.6583 * sinDegrees(2 * (getCorrectedLongitude() - sunLongitude));
        }

        /**
         * l'
         *
         * @return
         */
        public double getCorrectedLongitude() {
            return moonMeanLongitude
                    + getEvictionCorrection()
                    + getCorrectionForEquationCentre()
                    - getAnnualEquationCorrection()
                    + getYetAnotherCorrectionTerm();
        }

        /**
         * A 4
         *
         * @return
         */
        private double getYetAnotherCorrectionTerm() {
            return 0.214 * sinDegrees(2 * getMoonCorrectedAnomaly());
        }

        /**
         * E c
         *
         * @return
         */
        private double getCorrectionForEquationCentre() {
            return 6.2886 * sinDegrees(getMoonCorrectedAnomaly());
        }

        /**
         * M' m
         *
         * @return
         */
        private double getMoonCorrectedAnomaly() {
            return moonMeanAnomaly + getEvictionCorrection() - getAnnualEquationCorrection() - getUnnamedThirdCorrection();
        }

        /**
         * E v
         *
         * @return
         */
        public double getEvictionCorrection() {
            double C = moonMeanLongitude - sunLongitude;
            return 1.2739 * sinDegrees(2.0 * C - moonMeanAnomaly);
        }

        /**
         * A e
         *
         * @return
         */
        public double getAnnualEquationCorrection() {
            return 0.1858 * sinDegrees(sunMeanAnomaly);
        }

        /**
         * A 3
         *
         * @return
         */
        public double getUnnamedThirdCorrection() {
            return 0.37 * sinDegrees(sunMeanAnomaly);
        }
    }

   public void calcalateMoonAzEl(Calendar cal) {
	   
	   //Conversion of date date and time 
	   
	   //get Julian Day-
	   //double JD = JulianDate.makeJulianDateUsingMyModified(cal);
	   double JD = JulianDate.getJulianDay(cal);
	   //convert to Julian Ephemeris Day
	   //double JDE = JD;
	   //find the time in Julian Centuries
	   double T = (JD - 2451545.0) / 36525.0 ;
	   
	   //testing
	   //T =-0.077221081451; //1992 April 12 at 0h
	   
	   if(DEBUG) Log.d("MP","Time in Julian centuries is " +T);


	   //calculate the angle terms (degrees 0-360)
	   
	   /* calculate Lprime
	   * The Moons mean longitude referred to the mean equinox of the date, and including the constant term of the effect of light-time
	   */
	   double Lp = 218.3164591 
			   + 481267.88134236*T 
			   - 0.0013268*T*T
			   + T*T*T/538841.
			   - T*T*T*T/65194000.;
	   
	   Lp = adjustTo360Range(Lp);
	   if(DEBUG) Log.d("MP","Angle Lp is " +new DecimalFormat("000.0").format(Lp)+ "\u00b0\n");
	   
	   
	   /* calculate D
	    * the Mean elongation of the Moon
	    */
	   double D = 297.8502042 
			   + 445267.1115168*T
			   -0.0016300*T*T
			   +T*T*T/545868.
			   -T*T*T*T/113065000.;
	   
	   D = adjustTo360Range(D);
	   if(DEBUG) Log.d("MP","Angle D is " +new DecimalFormat("000.0").format(D)+ "\u00b0\n");
	   
	   /* calculate M
	    * The Sun's mean anomaly
	    */
	   double M = 357.5291092 
			   + 35999.0502909*T
			   - 0.0001536*T*T
			   + T*T*T/24490000.;
	   
	   M = adjustTo360Range(M);
	   if(DEBUG) Log.d("MP","Angle M is " +new DecimalFormat("000.0").format(M)+ "\u00b0\n");
	   
	   /* calculate Mp
	    * The Moon's mean anomaly
	    */
	   double Mp = 134.9634114 
			   + 477198.8676313*T
			   + 0.0089970*T*T
			   + T*T*T/69699.
			   - T*T*T*T/14712000.;
	   
	   Mp = adjustTo360Range(Mp);
	   if(DEBUG) Log.d("MP","Angle Mp is " +new DecimalFormat("000.0").format(Mp)+ "\u00b0\n");
	   
	   /* calculate F
	    * The Moon's argument of latitude (mean distance of the Moon from its ascending node)
	    */
	   double F = 93.2720993
			   + 483202.0175273*T
			   - 0.0034029*T*T
			   - T*T*T/3526000.
			   + T*T*T*T/863310000.;
	   
	   F = adjustTo360Range(F);
	   if(DEBUG) Log.d("MP","Angle F is " +new DecimalFormat("000.0").format(F)+ "\u00b0\n");
	   
	   /* Three further arguments: A1, A2, A3
	    */
	   double A1 = 119.75 + 131.849*T;  //due to the action of Venus
	   double A2 = 53.09 + 479264.290*T;  //due to the action of Jupiter
	   double A3 = 313.45 + 481266.484*T;  //due to flattening of the Earth
	   
	   
	   /* calculate E
	    * to take into account the decreasing eccentricity of the earth's orbit over time
	    */
	   double E = 1.0 - 0.002516*T - 0.0000074*T*T;
	   
	   /* Now, calculate the sums sigmaL and sigmaR and sigmaB
	    * using the first few terms (plenty more available to add later)
	    */
	   double sigmaL = 0.0;
	   sigmaL = 6288774.0 * sinDegrees(Mp);
	   sigmaL = sigmaL + 127402.*sinDegrees(2*D -Mp);
	   sigmaL = sigmaL + 658314.*sinDegrees(2*D);
	   sigmaL = sigmaL + 213618.*sinDegrees(2*Mp);
	   sigmaL = sigmaL + -185116.*sinDegrees(M) * E;
	   sigmaL = sigmaL + -114332.*sinDegrees(2*F);
	   sigmaL = sigmaL + 58793.*sinDegrees(2*D - 2*Mp);
	   sigmaL = sigmaL + 57066.*sinDegrees(2*D -M -Mp) * E;
	   sigmaL = sigmaL + 53322.*sinDegrees(2*D -Mp) ;
	   sigmaL = sigmaL + 45758.*sinDegrees(2*D -M ) * E;
	   if(DEBUG) Log.d("MP","sigma L (10 terms):  " + new DecimalFormat("00.0000").format(sigmaL) );
	   sigmaL = sigmaL + -40923.*sinDegrees(M -Mp) * E;
	   sigmaL = sigmaL + -34720.*sinDegrees(D );
	   sigmaL = sigmaL + -30383.*sinDegrees(M +Mp ) * E;
	   sigmaL = sigmaL + 15327.*sinDegrees(2*D +2*F) ;
	   sigmaL = sigmaL + -12528.*sinDegrees(Mp +2*F);
	   sigmaL = sigmaL + 10980.*sinDegrees(Mp -2*F) ;
	   sigmaL = sigmaL + 10675.*sinDegrees(4*D -Mp) ;
	   sigmaL = sigmaL + 10034.*sinDegrees(3*Mp) ;
	   sigmaL = sigmaL + 8548.*sinDegrees(4*D -2*Mp) ;
	   sigmaL = sigmaL + -7888.*sinDegrees(2*D +M -Mp) * E;
	   sigmaL = sigmaL + -6766.*sinDegrees(2*D +M ) * E;
	   sigmaL = sigmaL + -5163.*sinDegrees(D -Mp) ;
	   sigmaL = sigmaL + 4987.*sinDegrees(D +M ) ;
	   sigmaL = sigmaL + 4036.*sinDegrees(2*D -M +Mp) * E;
	   sigmaL = sigmaL + 3994.*sinDegrees(2*D -2*Mp) ;
	   sigmaL = sigmaL + 3861.*sinDegrees(4*D ) ;
	   sigmaL = sigmaL + 3665.*sinDegrees(2*D -3*Mp) ;
	   sigmaL = sigmaL + -2689.*sinDegrees(M -2*Mp) * E;
	   sigmaL = sigmaL + -2602.*sinDegrees(2*D -Mp +2*F);
	   sigmaL = sigmaL + 2390.*sinDegrees(2*D -M -2*Mp) * E;
	   sigmaL = sigmaL + -2348.*sinDegrees(D +Mp) * E;
	   sigmaL = sigmaL + 2236.*sinDegrees(2*D -2*M ) * E*E;
	   sigmaL = sigmaL + -2120.*sinDegrees(M -2*Mp) * E;
	   sigmaL = sigmaL + -2069.*sinDegrees(2*M )* E* E;
	   sigmaL = sigmaL + 2048.*sinDegrees(2*D -2*M - Mp )* E* E;
	   sigmaL = sigmaL + -1773.*sinDegrees(2*D -Mp - 2*F );
	   sigmaL = sigmaL + -1595.*sinDegrees(2*D -2*F );
	   sigmaL = sigmaL + 1215.*sinDegrees(4*D -M -Mp )* E;
	   sigmaL = sigmaL + -1110.*sinDegrees(2*Mp * 2*F );
	   sigmaL = sigmaL + -892.*sinDegrees(3*D -Mp );
	   sigmaL = sigmaL + -810.*sinDegrees(2*D +M +Mp )* E;
	   sigmaL = sigmaL + 759.*sinDegrees(4*D -M -2*Mp )* E;
	   sigmaL = sigmaL + -713.*sinDegrees(2*M -Mp )* E* E;
	   sigmaL = sigmaL + -700.*sinDegrees(2*D +2*M -Mp )* E* E;
	   sigmaL = sigmaL + 691.*sinDegrees(2*D +M -2*M )* E;
	   sigmaL = sigmaL + 596.*sinDegrees(2*D -M -2*F )* E;
	   sigmaL = sigmaL + 549.*sinDegrees(4*D +Mp);
	   sigmaL = sigmaL + 537.*sinDegrees(4*Mp );
	   sigmaL = sigmaL + 520.*sinDegrees(4*D -M )* E;
	   sigmaL = sigmaL + -487.*sinDegrees(D -2*Mp );
	   sigmaL = sigmaL + -399.*sinDegrees(2*D -M -2*F )* E;
	   sigmaL = sigmaL + -381.*sinDegrees(2*Mp -2*F );
	   sigmaL = sigmaL + 351.*sinDegrees(D +M +Mp )* E;
	   sigmaL = sigmaL + -340.*sinDegrees(3*D -2*M )* E* E;
	   sigmaL = sigmaL + 330.*sinDegrees(4*D -3*Mp );
	   sigmaL = sigmaL + 327.*sinDegrees(2*D -M +2*Mp )* E;
	   sigmaL = sigmaL + -323.*sinDegrees(2*M +Mp )* E* E;
	   sigmaL = sigmaL + 299.*sinDegrees(D +M -Mp )* E;
	   sigmaL = sigmaL + 294.*sinDegrees(2*D +3*Mp );
	   sigmaL = sigmaL + 0.*sinDegrees(2*D -Mp -2*F );
	   if(DEBUG) Log.d("MP","sigma L (all terms):  " + new DecimalFormat("00.0000").format(sigmaL) ); 


	   //add in the A factors
	   sigmaL = sigmaL 
			   + 3958. * sinDegrees(A1)
	   		   + 1962. * sinDegrees(Lp - F)
	   		   + 318. * sinDegrees(A2);
	   
	   if(DEBUG) Log.d("MP","sigma L with A terms: " +sigmaL);

	   
	   //sigma R is used in calculating the distance to the Moon. 
	   double sigmaR = 0.0;
	   sigmaR = -20905355. * cosDegrees(Mp);
	   sigmaR = sigmaR + -3699111.*cosDegrees(2*D -Mp);
	   sigmaR = sigmaR + -2955968.*cosDegrees(2*D);
	   sigmaR = sigmaR + -569925.*cosDegrees(2*Mp);
	   sigmaR = sigmaR + 48888.*cosDegrees(M) * E;
	   sigmaR = sigmaR + -3149.*cosDegrees(2*F);
	   sigmaR = sigmaR + 246158.*cosDegrees(2*D - 2*Mp);
	   sigmaR = sigmaR + -152138.*cosDegrees(2*D -M -Mp) * E;
	   sigmaR = sigmaR + -170733.*cosDegrees(2*D -Mp) ;
	   sigmaR = sigmaR + -204586.*cosDegrees(2*D -M ) * E;
	   if(DEBUG) Log.d("MP","sigma R (10 terms):  " + new DecimalFormat("00.0000").format(sigmaR) );
	   sigmaR = sigmaR + -129620.*cosDegrees(M -Mp) * E;
	   sigmaR = sigmaR + 108743.*cosDegrees(D );
	   sigmaR = sigmaR + 104755.*cosDegrees(M +Mp ) * E;
	   sigmaR = sigmaR + 10321.*cosDegrees(2*D +2*F) ;
	   sigmaR = sigmaR + 0.*cosDegrees(Mp +2*F);
	   sigmaR = sigmaR + 79661.*cosDegrees(Mp -2*F) ;
	   sigmaR = sigmaR + -34782.*cosDegrees(4*D -Mp) ;
	   sigmaR = sigmaR + -23210.*cosDegrees(3*Mp) ;
	   sigmaR = sigmaR + -21636.*cosDegrees(4*D -2*Mp) ;
	   sigmaR = sigmaR + 24208.*cosDegrees(2*D +M -Mp) * E;
	   sigmaR = sigmaR + 30824.*cosDegrees(2*D +M ) * E;
	   sigmaR = sigmaR + -8379.*cosDegrees(D -Mp) ;
	   sigmaR = sigmaR + -16675.*cosDegrees(D +M ) ;
	   sigmaR = sigmaR + -12831.*cosDegrees(2*D -M +Mp) * E;
	   sigmaR = sigmaR + -10445.*cosDegrees(2*D -2*Mp) ;
	   sigmaR = sigmaR + -11650.*cosDegrees(4*D ) ;
	   sigmaR = sigmaR + 14403.*cosDegrees(2*D -3*Mp) ;
	   sigmaR = sigmaR + -7003.*cosDegrees(M -2*Mp) * E;
	   sigmaR = sigmaR + 0.*cosDegrees(2*D -Mp +2*F);
	   sigmaR = sigmaR + 10056.*cosDegrees(2*D -M -2*Mp) * E;
	   sigmaR = sigmaR + 6322.*cosDegrees(D +Mp) * E;
	   sigmaR = sigmaR + -9884.*cosDegrees(2*D -2*M ) * E*E;
	   sigmaR = sigmaR + 5751.*cosDegrees(M -2*Mp) * E;
	   sigmaR = sigmaR + 0.*cosDegrees(2*M )* E* E;
	   sigmaR = sigmaR + -4950.*cosDegrees(2*D -2*M - Mp )* E* E;
	   sigmaR = sigmaR + 4130.*cosDegrees(2*D -Mp - 2*F );
	   sigmaR = sigmaR + 0.*cosDegrees(2*D -2*F );
	   sigmaR = sigmaR + -3958.*cosDegrees(4*D -M -Mp )* E;
	   sigmaR = sigmaR + 0.*cosDegrees(2*Mp * 2*F );
	   sigmaR = sigmaR + 3258.*cosDegrees(3*D -Mp );
	   sigmaR = sigmaR + 2616.*cosDegrees(2*D +M +Mp )* E;
	   sigmaR = sigmaR + -1897.*cosDegrees(4*D -M -2*Mp )* E;
	   sigmaR = sigmaR + -2117.*cosDegrees(2*M -Mp )* E* E;
	   sigmaR = sigmaR + 2354.*cosDegrees(2*D +2*M -Mp )* E* E;
	   sigmaR = sigmaR + 0.*cosDegrees(2*D +M -2*M )* E;
	   sigmaR = sigmaR + 0.*cosDegrees(2*D -M -2*F )* E;
	   sigmaR = sigmaR + -1423.*cosDegrees(4*D +Mp);
	   sigmaR = sigmaR + -1117.*cosDegrees(4*Mp );
	   sigmaR = sigmaR + -1571.*cosDegrees(4*D -M )* E;
	   sigmaR = sigmaR + -1739.*cosDegrees(D -2*Mp );
	   sigmaR = sigmaR + 0.*cosDegrees(2*D -M -2*F )* E;
	   sigmaR = sigmaR + -4421.*cosDegrees(2*Mp -2*F );
	   sigmaR = sigmaR + 0.*cosDegrees(D +M +Mp )* E;
	   sigmaR = sigmaR + 0.*cosDegrees(3*D -2*M )* E* E;
	   sigmaR = sigmaR + 0.*cosDegrees(4*D -3*Mp );
	   sigmaR = sigmaR + 0.*cosDegrees(2*D -M +2*Mp )* E;
	   sigmaR = sigmaR + 1165.*cosDegrees(2*M +Mp )* E* E;
	   sigmaR = sigmaR + 0.*cosDegrees(D +M -Mp )* E;
	   sigmaR = sigmaR + 0.*cosDegrees(2*D +3*Mp );
	   sigmaR = sigmaR + 8752.*cosDegrees(2*D -Mp -2*F );
	   if(DEBUG) Log.d("MP","sigma R (all terms):  " + new DecimalFormat("00.0000").format(sigmaR) ); 	   
	   
	   
	   double sigmaB = 0.0;
	   sigmaB = 5128122.0 * sinDegrees(F);
	   sigmaB = sigmaB + 280602.*sinDegrees(Mp +F);
	   sigmaB = sigmaB + 277693.*sinDegrees(Mp -F);
	   sigmaB = sigmaB + 173237.*sinDegrees(2*D - F);
	   sigmaB = sigmaB + 55413.*sinDegrees(2*D -Mp +F);
	   sigmaB = sigmaB + 46271.*sinDegrees(2*D -Mp -F);
	   sigmaB = sigmaB + 32573.*sinDegrees(2*D);
	   sigmaB = sigmaB + 17198.*sinDegrees(2*Mp );
	   sigmaB = sigmaB + 9266.*sinDegrees(2*D + Mp -F);
	   sigmaB = sigmaB + 8822.*sinDegrees(2*Mp -F);
	   sigmaB = sigmaB + 8216.*sinDegrees(2*D -M -F)* E;
	   sigmaB = sigmaB + 4324.*sinDegrees(2*D -2*Mp -F);
	   sigmaB = sigmaB + 4200.*sinDegrees(2*D +Mp +F);
	   sigmaB = sigmaB + -3359.*sinDegrees(2*D +M -F)* E;
	   sigmaB = sigmaB + 2463.*sinDegrees(2*D -M -Mp +F)* E;
	   sigmaB = sigmaB + 2211.*sinDegrees(2*D -M +F)* E;
	   sigmaB = sigmaB + 2065.*sinDegrees(2*D -M -Mp -F)* E;
	   sigmaB = sigmaB + -1870.*sinDegrees(M -Mp +F);
	   sigmaB = sigmaB + 1828.*sinDegrees(4*D -Mp -F);
	   sigmaB = sigmaB + -1794.*sinDegrees(M+F)* E;
	   sigmaB = sigmaB + -1749.*sinDegrees(3*F);
	   sigmaB = sigmaB + -1565.*sinDegrees(M-Mp+F)* E;
	   sigmaB = sigmaB + -1491.*sinDegrees(D+F);
	   sigmaB = sigmaB + -1475.*sinDegrees(M +Mp +F)* E;
	   sigmaB = sigmaB + -1410.*sinDegrees(M +Mp -F)* E;
	   sigmaB = sigmaB + -1344.*sinDegrees(M -F)* E;
	   sigmaB = sigmaB + -1335.*sinDegrees(D -F);
	   sigmaB = sigmaB + 1107.*sinDegrees(3*Mp +F);
	   sigmaB = sigmaB + 1021.*sinDegrees(4*D -F);
	   sigmaB = sigmaB + 833.*sinDegrees(4*D -Mp +F);
	   sigmaB = sigmaB + 777*sinDegrees(Mp-3*F);
	   sigmaB = sigmaB + 671.*sinDegrees(4*D -2*Mp +F);
	   sigmaB = sigmaB + 607.*sinDegrees(2*D -3*F);
	   sigmaB = sigmaB + 596.*sinDegrees(2*D +2*Mp - F);
	   sigmaB = sigmaB + 491.*sinDegrees(2*D -M +Mp -F)* E;
	   sigmaB = sigmaB  -451.*sinDegrees(2*D -2*Mp +F);
	   sigmaB = sigmaB + 439.*sinDegrees(3*Mp -F);
	   sigmaB = sigmaB + 422.*sinDegrees(2*D +2*Mp +F);
	   sigmaB = sigmaB + 421.*sinDegrees(2*D -3*Mp -F);
	   sigmaB = sigmaB  -366.*sinDegrees(2*D +M -Mp +F)* E;
	   sigmaB = sigmaB  -351.*sinDegrees(2*D +M +F)* E;
	   sigmaB = sigmaB + 331.*sinDegrees(4*D +F);
	   sigmaB = sigmaB + 315.*sinDegrees(2*D -M +Mp +F)* E;
	   sigmaB = sigmaB + 302.*sinDegrees(2*D -2*M -F)* E* E;
	   sigmaB = sigmaB  -283.*sinDegrees(Mp +3*F);
	   sigmaB = sigmaB  -229.*sinDegrees(2*D +M +Mp -F)* E;
	   sigmaB = sigmaB + 223.*sinDegrees(D +M -F)* E;
	   sigmaB = sigmaB + 223.*sinDegrees(D +M +F)* E;
	   sigmaB = sigmaB  -220.*sinDegrees(M -2*Mp -F)* E;
	   sigmaB = sigmaB  -220.*sinDegrees(2*D +M -Mp -F)* E;
	   sigmaB = sigmaB  -185.*sinDegrees(D +Mp +F);
	   sigmaB = sigmaB + 181.*sinDegrees(2*D -M -2*Mp -F)* E;
	   sigmaB = sigmaB  -177.*sinDegrees(M +2*Mp +F)* E;
	   sigmaB = sigmaB + 176.*sinDegrees(4*D -2*Mp -F);
	   sigmaB = sigmaB + 166.*sinDegrees(4*D -M -Mp -F)* E;
	   sigmaB = sigmaB  -164.*sinDegrees(D -Mp -F);
	   sigmaB = sigmaB + 132.*sinDegrees(4*D +Mp -F);
	   sigmaB = sigmaB  -119.*sinDegrees(D -Mp -F);
	   sigmaB = sigmaB + 115.*sinDegrees(4*D -M -F)* E;
	   sigmaB = sigmaB + 107.*sinDegrees(2*D -2*M +F)* E* E;


	   //add in the A factors
	   sigmaB = sigmaB 
			   - 2235. * sinDegrees(Lp)
	   		   + 382. * sinDegrees(A3)
	   		   + 175. * sinDegrees(A1 - F)
	   		   + 175. * sinDegrees(A1 + F)
	   		   + 127. * sinDegrees(Lp - Mp)
	   		   - 115. * sinDegrees(Lp + Mp);
	   
	   
	   if(DEBUG) Log.d("MP","sigma B with A terms: " +sigmaB);

	   
	   /* The coordinates of the Moon are given by 
	    * geocentric longitude lambda (degrees)
	    * geocentric latitude beta (degrees)
	    * distance between earth and moon centers delta (km)
	    */
	   moonLambda = Lp + sigmaL/1000000.0;   	//in deg
	   moonBeta = sigmaB/1000000.0;   			//in deg
	   moonDelta = 385000.56 + sigmaR/1000.0;  	//in km
	   
	   
	   /* Next Step 
	    * convert ecliptic latitude and longitude to right ascension RA and declination delta 
	    */
	   
	   double lambda = moonLambda; //ecliptic longitude
	   double beta = moonBeta; //ecliptic latitude

	   //obliquity of ecliptic:
	   double eps = 23.0 + 26.0/60.0 + 21.448/3600.0;
	   eps = eps -46.8150*T/3600.;
	   eps = eps - 0.00059*T*T/3600.;
	   eps = eps + 0.001813*T*T*T/3600.;
	   if(DEBUG) Log.d("MP","obliquity of ecliptic: " +eps);

	   
	   double X = cosDegrees(beta)*cosDegrees(lambda);
	   double Y = cosDegrees(eps)*cosDegrees(beta)*sinDegrees(lambda) - sinDegrees(eps)*sinDegrees(beta);
	   double Z = sinDegrees(eps)*cosDegrees(beta)*sinDegrees(lambda) - cosDegrees(eps)*sinDegrees(beta);
	   double R = Math.sqrt(1.0-Z*Z);

	   //right ascension RA and declination delta
	   double ddelta = (180./Math.PI)*Math.atan2(Z,R); // in degrees
	   //double RA = (24./Math.PI)*Math.atan(Y/(X+R)); // in hours	   
	   double RA = (24./360.)*adjustTo360Range((180./Math.PI)*(Math.atan2(Y,(X+R)))); // in hours	   

//	   if(DEBUG) Log.d("MP","Moon Right Ascension: " + new DecimalFormat("00.00").format(RA) +" hrs");
//	   if(DEBUG) Log.d("MP","Moon Declination: " + new DecimalFormat("#00.00").format(ddelta) +" deg");
//

	   /* Next Step
	    * compute sidereal time (degrees) at Greenwich 
	    */
//	   double mlong = -97.5;  //my geographic longitude
//	   double beta = 30.17;  //my geographic latitude
	   
	   //double T = (JD - 2451545.0 ) / 36525;

//	   double theta = theta0 + mlong;  //(eastern longitudes positive, western negative)
//	   double tau = theta - RA;
	   
	   //local hour angle
	   //tau = adjustTo360Range(tau);  
	   
	   /* Final Result
	    * convert (tau, delta) to horizon coordinates (h, az) of the observer (geo. latitude beta, geo. logitude mlong) 
	    * 	   sin(h) = sin(beta )*sin(delta) + cos(beta)*cos(delta)*cos(tau)
	    *      tan(az) = - sin(tau) / (cos(beta)*tan(delta) - sin(beta)*cos(tau))
	    */
	   
	   
//	   //ecliptic longitude
//	   //eq. (12.1) 
//	   double lambda = 180./Math.PI*Math.atan( 
//			   (sinDegrees(alpha)*cosDegrees(epsilon) + tanDegrees(delta)*sinDegrees(epsilon))
//			   / cosDegrees(alpha) );
//	 
//	   //ecliptic latitude
//	   //eq. (12.2)
//	   double beta = 180./Math.PI*Math.asin( 
//			   sinDegrees(delta)*cosDegrees(epsilon) 
//			   - cosDegrees(delta)*sinDegrees(epsilon)*sinDegrees(alpha));
	   
	   
	   
	   //30.2500° N, 97.7500° W
	   double phi = 30.25 ; //observer's latitude positive in Northern hemi, neg in Southern hemi
	   double L = -90.75; //observer's longitude, west-negative 

	   //epsilon = obliquity of the ecliptic
	   //eq 21.2 (mean obliquity of ecliptic), T is in julian centuries from the epoch J2000
	   //double epsilon = 23.439291 - 0.013004*T;
	   double epsilon = 23.0 + 26.0/60.0 + 21.448/3600.0 - (46.8150*T+ 0.00059*T*T- 0.001813*T*T*T)/3600.;

	   //the sidereal time at greenwich
	   double theta0 = 280.46061837 ;
	   theta0 = theta0	+ 360.98564736629*(JD-2451545.0) ;
	   theta0 = theta0	+ 0.000387933*T*T ;
	   theta0 = theta0	- T*T*T/38710000.0; // degrees 

	   //alpha = the right ascension
	   //eq. (12.3)
	   double alpha = Math.atan2( 
			   ( sinDegrees(lambda)*cosDegrees(epsilon) - tanDegrees(beta)*sinDegrees(epsilon)) 
			   , (cosDegrees(lambda)) );
	   alpha = adjustTo360Range(alpha*180./Math.PI);
	   if(DEBUG) Log.d("MP","Moon's right ascension: " +new DecimalFormat("#00.00").format(alpha/360.*24.)+ "hrs");

	   
	   
	   //H: local hour angle
	   //double H = theta - alpha;
	   double H = theta0 - L + alpha;
	   
	   //delta = declination, positive if north of the celestial equator
	   //eq. (12.4)
	   double delta = Math.asin(
			   sinDegrees(beta)*cosDegrees(epsilon) 
			   + cosDegrees(beta)*sinDegrees(epsilon)*sinDegrees(lambda) );
	   //delta = adjustTo360Range(delta*180./Math.PI);
	   delta = delta*180./Math.PI;
	   if(DEBUG) Log.d("MP","Moon's declination: " +new DecimalFormat("#00.00").format(delta)+ "\u00b0");

	   
//	   //eq. (12.5) //A=azimuth, h=attitude(elevation)
	   double A = Math.atan2(sinDegrees(H),cosDegrees(H)*sinDegrees(phi) - tanDegrees(delta)*cosDegrees(phi));
	   A = A * 180./Math.PI + 180.; //convert to degrees, +180 to reference CW from N instead of south
	   //eq. (12.6)
	   double h = Math.asin(sinDegrees(phi)*sinDegrees(delta) + cosDegrees(phi)*cosDegrees(delta)*cosDegrees(H));
	   h = h * 180./Math.PI; //convert to degrees
	   
//		//Find the h and AZ at the current LST
//		double A = Math.acos( (-sinDegrees(delta) + sinDegrees(h)*sinDegrees(L) )/ (cosDegrees(h)*cosDegrees(L)) );
//		A = A*(180/Math.PI);
//	   
	   /**/
	   
	   //double h = ( 180./Math.PI*Math.asin(sinDegrees(B )*sinDegrees(delta) + cosDegrees(beta)*cosDegrees(delta)*cosDegrees(tau) ) );
	   //double az = adjustTo360Range( 180./Math.PI*Math.atan2( - sinDegrees(tau), (cosDegrees(beta)*tanDegrees(delta) - sinDegrees(B)*cosDegrees(tau)) ) );
	   
	   if(DEBUG) Log.d("MP","H calculated: " + h);

	   /* Calculate for parallax in altitude
	    */
	   double horParal = 8.794 / (moonDelta/149.59787E6); // horizontal parallax (arcseconds), Meeus S. 263
	   double p = 180./Math.PI*Math.asin(cosDegrees(h)*sinDegrees(horParal/3600.0)); // parallax in altitude (degrees)
	   
	   /*The refraction 'ref' is calculated by Saemundsson's formula (Meeus, Astronomical Algorithms):
	    *    h is the true (airless) altitude in degrees, 
	    *    'ref' is in minutes of arc. 
	    *    The apparent (measured) altitude is h+R. 
	    */
	   double ref = 1.02 / (tanDegrees(h + 10.3/(h + 5.11)));
	   
	   /* apparent altitude h (corrected for parallax and refraction) 
	    */
	   double happarent = h + ref + p;
	   if(DEBUG) Log.d("MP","H apparent: " + happarent);
	   
	   //stash final values
	   moonElevationApparent = happarent;
	   moonAzimuth = A;
	   
	   if(DEBUG) Log.d("MP","TEST lambda: " + lambda);
	   if(DEBUG) Log.d("MP","TEST beta: " + beta);
	   if(DEBUG) Log.d("MP","TEST Delta: " + moonDelta +"km");
	   if(DEBUG) Log.d("MP","TEST Lp: " + Lp);
	   if(DEBUG) Log.d("MP","TEST D: " + D);
	   if(DEBUG) Log.d("MP","TEST M: " + M);
	   if(DEBUG) Log.d("MP","TEST Mp: " + Mp);
	   if(DEBUG) Log.d("MP","TEST F: " + F);
	   if(DEBUG) Log.d("MP","TEST epsilon: " + epsilon);
	   if(DEBUG) Log.d("MP","TEST alpha: " + alpha);
	   if(DEBUG) Log.d("MP","TEST delta: " + delta);
	   
/* for T=-0.077221081451
 * should see -
 * Lp = 134
 * D = 113
 * M = 97
 * Mp = 5
 * F = 219
 * lambda = 133deg
 * beta = -3
 * Ddelta = 368405km
 * delta 13deg
 * alpha = 134deg
 * 
 */

	   
   }
   
   
   public void calcalateMoonAzEl2(Calendar cal) {
	   
	   //Conversion of date date and time 
	   
	   //get Julian Day-
	   //double JD = JulianDate.makeJulianDateUsingMyModified(cal);
	   double JD = JulianDate.getJulianDay(cal);
	   //convert to Julian Ephemeris Day
	   double JDE = JD;
	   //find the time in Julian Centuries
	   double T = (JDE - 2451545.0) / 36525.0 ;
	   
	   if(DEBUG) Log.d("MP","Time in Julian centuries is " +T);


	   //calculate the angle terms (degrees 0-360)
	   
	   /* calculate Lprime
	   * The Moons mean longitude referred to the mean equinox of the date, and including the constant term of the effect of light-time
	   */
	   double Lp = 218.3164591 
			   + 481267.88134236*T 
			   - 0.0013268*T*T
			   + T*T*T/65194000.
			   - T*T*T*T/65194000.;
	   
	   Lp = adjustTo360Range(Lp);
	   if(DEBUG) Log.d("MP","Angle Lp is " +Lp);
	   
	   
	   /* calculate D
	    * the Mean elongation of the Moon
	    */
	   double D = 297.8502042 
			   + 445267.1115168*T
			   -0.0016300*T*T
			   +T*T*T/545868.
			   -T*T*T*T/113065000.;
	   
	   D = adjustTo360Range(D);
	   if(DEBUG) Log.d("MP","Angle D is " +D);
	   
	   /* calculate M
	    * The Sun's mean anomaly
	    */
	   double M = 357.5291092 
			   + 35999.0502909*T
			   - 0.0001536*T*T
			   + T*T*T/24490000.;
	   
	   M = adjustTo360Range(M);
	   if(DEBUG) Log.d("MP","Angle M is " +M);
	   
	   /* calculate Mp
	    * The Moon's mean anomaly
	    */
	   double Mp = 134.9634114 
			   + 477198.8676313*T
			   + 0.0089970*T*T
			   + T*T*T/69699.
			   - T*T*T*T/14712000.;
	   
	   Mp = adjustTo360Range(Mp);
	   if(DEBUG) Log.d("MP","Angle Mp is " +Mp);
	   
	   /* calculate F
	    * The Moon's argument of latitude (mean distance of the Moon from its ascending node)
	    */
	   double F = 93.2720993
			   + 483202.0175273*T
			   - 0.0034029*T*T
			   - T*T*T/3526000.
			   + T*T*T*T/863310000.;
	   
	   F = adjustTo360Range(F);
	   if(DEBUG) Log.d("MP","Angle F is " +F);
	   
	   /* Three further arguments: A1, A2, A3
	    */
	   double A1 = 119.75 + 131.849*T;  //due to the action of Venus
	   double A2 = 53.09 + 479264.290*T;  //due to the action of Jupiter
	   double A3 = 313.45 + 481266.484*T;  //due to flattening of the Earth
	   
	   
	   /* calculate E
	    * to take into account the descreasing eccentricity of the earth's orbit over time
	    */
	   double E = 1.0 - 0.002516*T - 0.0000074*T*T;
	   
	   /* Now, calculate the sums sigmaL and sigmaR and sigmaB
	    * using the first few terms (plenty more available to add later)
	    */
	   double sigmaL = 0.0;
	   sigmaL = 6288774.0 * sinDegrees(Mp);
	   sigmaL = sigmaL + 127402.*sinDegrees(2*D -Mp);
	   sigmaL = sigmaL + 658314.*sinDegrees(2*D);
	   sigmaL = sigmaL + 213618.*sinDegrees(2*Mp);
	   sigmaL = sigmaL + -185116.*sinDegrees(M) * E;
	   sigmaL = sigmaL + -114332.*sinDegrees(2*F);
	   sigmaL = sigmaL + 58793.*sinDegrees(2*D - 2*Mp);
	   sigmaL = sigmaL + 57066.*sinDegrees(2*D -M -Mp) * E;

	   if(DEBUG) Log.d("MP","sigmaL is " +sigmaL);

	   //add in the A factors
	   sigmaL = sigmaL 
			   + 3958. * sinDegrees(A1)
	   		   + 1962. * sinDegrees(Lp - F)
	   		   + 318. * sinDegrees(A2);
	   
	   if(DEBUG) Log.d("MP","sigmaL is " +sigmaL);

	   
	   double sigmaR = 0.0;
	   sigmaR = -20905355. * cosDegrees(Mp);
	   sigmaR = sigmaR + -3699111.*cosDegrees(2*D -Mp);
	   sigmaR = sigmaR + -2955968.*cosDegrees(2*D);
	   sigmaR = sigmaR + -569925.*cosDegrees(2*Mp);
	   sigmaR = sigmaR + 48888.*cosDegrees(M) * E;
	   sigmaR = sigmaR + -3149.*cosDegrees(2*F);
	   sigmaR = sigmaR + 246158.*cosDegrees(2*D - 2*Mp);
	   sigmaR = sigmaR + -152138.*cosDegrees(2*D -M -Mp) * E;
	   
	   
	   //add in the A factors
	   sigmaR = sigmaR 
			   - 2235. * sinDegrees(Lp)
	   		   + 382. * sinDegrees(A3)
	   		   + 175. * sinDegrees(A1 - F)
	   		   + 175. * sinDegrees(A1 + F)
	   		   + 127. * sinDegrees(Lp - Mp)
	   		   - 115. * sinDegrees(Lp + Mp);
	   
	   double sigmaB = 0.0;
	   sigmaB = 5128122.0 * sinDegrees(F);
	   sigmaB = sigmaB + 280602.*sinDegrees(Mp +F);
	   sigmaB = sigmaB + 277693.*sinDegrees(Mp -F);
	   sigmaB = sigmaB + 173237.*sinDegrees(2*D - F);
	   sigmaB = sigmaB + 55413.*sinDegrees(2*D -Mp +F);
	   sigmaB = sigmaB + 46271.*sinDegrees(2*D -Mp -F);
	   sigmaB = sigmaB + 32573.*sinDegrees(2*D);
	   sigmaB = sigmaB + 17198.*sinDegrees(2*Mp +F);
	   
	   
	   
	   /* The coordinates of the Moon are given by 
	    * geocentric longitude lambda (degrees)
	    * geocentric latitude beta (degrees)
	    * distance between earth and moon centers delta (km)
	    */
	   moonLambda = Lp + sigmaL/100000.0;
	   moonBeta = sigmaB/1000000.0;
	   moonDelta = 385000.56 + sigmaR/1000.0;
	   
	   
	   /* Next Step 
	    * convert ecliptic latitude and longitude to right ascension RA and declination delta 
	    * 
	    */
	   //double T = (JD-2451545.0)/36525.0;

	   double L = moonLambda; //ecliptic longitude
	   double B = moonBeta; //ecliptic latitude

	   //obliquity of ecliptic:
	   double eps = 23.0 + 26.0/60.0 + 21.448/3600.0 - (46.8150*T+ 0.00059*T*T- 0.001813*T*T*T)/3600;

	   double X = cosDegrees(B)*cosDegrees(L);
	   double Y = cosDegrees(eps)*cosDegrees(B)*sinDegrees(L) - sinDegrees(eps)*sinDegrees(B);
	   double Z = sinDegrees(eps)*cosDegrees(B)*sinDegrees(L) - cosDegrees(eps)*sinDegrees(B);
	   double R = Math.sqrt(1.0-Z*Z);

	   //right ascension RA and declination delta
	   double delta = (180/Math.PI)*Math.atan(Z/R); // in degrees
	   double RA = (24/Math.PI)*Math.atan(Y/(X+R)); // in hours	   
	   

	   /* Next Step
	    * compute sidereal time (degrees) at Greenwich 
	    */
	   double mlong = -97.0;  //my geographic longitude
	   double beta = 30.17;  //my geographic latitude
	   
	   //double T = (JD - 2451545.0 ) / 36525;
	   double theta0 = 280.46061837 
			   + 360.98564736629*(JD-2451545.0) 
			   + 0.000387933*T*T 
			   - T*T*T/38710000.0; // degrees 
	   double theta = theta0 + mlong;  //(eastern longitudes positive, western negative)
	   double tau = theta - RA;
	   
	   //local hour angle
	   tau = adjustTo360Range(tau);  
	   
	   /* Final Result
	    * convert (tau, delta) to horizon coordinates (h, az) of the observer (geo. latitude beta, geo. logitude mlong) 
	    * 	   sin(h) = sin(beta )*sin(delta) + cos(beta)*cos(delta)*cos(tau)
	    *      tan(az) = - sin(tau) / (cos(beta)*tan(delta) - sin(beta)*cos(tau))
	    */
	   double h = ( 180./Math.PI*Math.asin(sinDegrees(B )*sinDegrees(delta) + cosDegrees(beta)*cosDegrees(delta)*cosDegrees(tau) ) );
	   double az = adjustTo360Range( 180./Math.PI*Math.atan2( - sinDegrees(tau), (cosDegrees(beta)*tanDegrees(delta) - sinDegrees(B)*cosDegrees(tau)) ) );
	   
	   if(DEBUG) Log.d("MP","H calculated: " + h);

	   /* Calculate for parallax in altitude
	    */
	   double horParal = 8.794 / (moonDelta/149.59787E6); // horizontal parallax (arcseconds), Meeus S. 263
	   double p = 180./Math.PI*Math.asin(cosDegrees(h)*sinDegrees(horParal/3600.0)); // parallax in altitude (degrees)
	   
	   /*The refraction 'ref' is calculated by Saemundsson's formula (Meeus, Astronomical Algorithms):
	    *    h is the true (airless) altitude in degrees, 
	    *    'ref' is in minutes of arc. 
	    *    The apparent (measured) altitude is h+R. 
	    */
	   double ref = 1.02 / (tanDegrees(h + 10.3/(h + 5.11)));
	   
	   /* apparent altitude h (corrected for parallax and refraction) 
	    */
	   double happarent = h + ref + p;
	   if(DEBUG) Log.d("MP","H apparent: " + happarent);
	   
	   //stash final values
	   moonElevationApparent = happarent;
	   moonAzimuth = az;
   }
   
   public double getMoonLambda() {
	   return moonLambda;
   }
   public double getMoonBeta() {
	   return moonBeta;
   }
   public double getMoonDelta() {
	   return moonDelta;
   }
   public double getMoonAzimuth() {
	   return moonAzimuth;
   }
   public double getMoonElevationApparent() {
	   return moonElevationApparent;
   }
   
}
