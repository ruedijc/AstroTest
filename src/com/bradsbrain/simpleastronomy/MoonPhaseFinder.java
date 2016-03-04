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

import java.util.Calendar;
import java.util.Date;

public class MoonPhaseFinder {

    private static final int NUM_STEP_MINUTES = 1;

    private static MoonChecker newMoonChecker = new NewMoonChecker();
    private static MoonChecker fullMoonChecker = new FullMoonChecker();

    public enum MoonPhase {
        NEW,
        WAXINGCRESCENT,
        FIRSTQUARTER,
        WAXINGGIBBOUS,
        FULL,
        WANINGGIBBOUS,
        LASTQUARTER,
        WANINGCRESCENT;


        static MoonPhase finder(double percent) {
            return FULL;
        }

    }

    /**
     * Someday this will return a descriptive MoonPhase enum when handed a cal/date
     *
     * @param cal
     * @return
     */
    public static MoonPhase findMoonPhaseAt(Calendar cal) {

        return null;
    }

    public static String getMoonPhaseName(Calendar cal) {
    	double angle = getMoonAngle(cal);
    	String phaseName = "";
    	
    	if(angle>-7 && angle <= 7)
    		phaseName = "New Moon";
    	else if(angle>7 && angle <= 67.5)
    		phaseName = "Waxing Cresent";
    	else if(angle>67.5 && angle <= 112.5)
    		phaseName = "First Quarter";
    	else if(angle>112.5 && angle <= 173)
    		phaseName = "Waxing Gibbous";
    	else if(angle>173 && angle <= 187)
    		phaseName = "Full Moon";
    	else if(angle>187 && angle <= 247.5)
    		phaseName = "Waning Gibbous";
    	else if(angle>247.5 && angle <= 292.5)
    		phaseName = "Last Quarter";
    	else if(angle>292.5 && angle <= 353)
    		phaseName = "Waning Crescent";
    	else if(angle>353 && angle <= 367)
    		phaseName = "New Moon";

    	return phaseName;
    }
    public static Date findFullMoonFollowing(Calendar cal) {
        return findDatePassingBounds(cal, fullMoonChecker);
    }

    public static Date findNewMoonFollowing(Calendar cal) {
        return findDatePassingBounds(cal, newMoonChecker);
    }

    public static Date findFullMoonFollowing(Date date) {
        long time = date.getTime();
        Calendar cal = Calendar.getInstance();
        cal.setTimeInMillis(time);
        return findDatePassingBounds(cal, fullMoonChecker);
    }

    public static Date findNewMoonFollowing(Date date) {
        long time = date.getTime();
        Calendar cal = Calendar.getInstance();
        cal.setTimeInMillis(time);
        return findDatePassingBounds(cal, newMoonChecker);
    }

    /**
     * TODO: This loop isn't very efficient, come up with a better implementation
     *
     * @param cal         the calendar date for which to compute the moon position
     * @param moonChecker the NewMoon or FullMoon checker
     * @return the forward date which passes the given bounds provided
     */
    private static Date findDatePassingBounds(Calendar cal, MoonChecker moonChecker) {
        Calendar myCal = BaseUtils.getSafeLocalCopy(cal.getTimeInMillis());
        Calendar thirtyOneDaysLater = Calendar.getInstance();
        thirtyOneDaysLater.setTimeInMillis(myCal.getTimeInMillis());
        thirtyOneDaysLater.add(Calendar.DAY_OF_MONTH, 31);
        // if we don't find a new moon after 31 days we're not going to find it. days between phases is ~29.5
        while (myCal.before(thirtyOneDaysLater)) {
            myCal.add(Calendar.MINUTE, NUM_STEP_MINUTES);
            double percent = 100 * MoonPhaseFinder.getMoonVisiblePercent(myCal);
            double angle = MoonPhaseFinder.getMoonAngle(myCal);
            if (moonChecker.isCorrectAngle(angle) && moonChecker.isCorrectPercent(percent)) {
                return myCal.getTime();
            }
        }
        return null;
    }

    /**
     * Returns a (much-too-)high-precision value for the amount of moon visible.
     * Value will be somewhere in the range 0% to 100%  (i.e. 0.00 to 1.00)
     *
     * @param cal
     * @return percent of moon which is visible
     */
    public static double getMoonVisiblePercent(Calendar cal) {
        double moonAngle = getMoonAngle(cal);
        return BaseUtils.useLessPrecision(0.5 * (1 - BaseUtils.cosDegrees(moonAngle)), 3);
    }

    /**
     * The moon angle.  For that we need the sun's position and the moon's position. </br>
     * The moon angle will be in the range 0 to 360.  <br/>
     * 0 or 360 is NEW, 180 is FULL
     *
     * @param cal
     * @return the angle of the moon in relation to the earth
     */
    public static double getMoonAngle(Calendar cal) {
        Calendar myCal = BaseUtils.getSafeLocalCopy(cal.getTimeInMillis());
        SunPosition sunPos = new SunPosition(myCal);
        MoonPosition moonPos = new MoonPosition(myCal);

        double angleAge = moonPos.getTrueLongitude() - sunPos.getEclipticLongitude();
        if (angleAge < 0) {
            return 360 + angleAge;
        }
        return angleAge;
    }

}
