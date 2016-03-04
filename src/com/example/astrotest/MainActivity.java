package com.example.astrotest;

import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.Locale;
import java.util.TimeZone;

import com.bradsbrain.simpleastronomy.MoonPhaseFinder;
import com.bradsbrain.simpleastronomy.MoonPosition;
import com.bradsbrain.simpleastronomy.SunPosition;



import android.app.Activity;
import android.content.Context;
import android.location.Location;
import android.location.LocationManager;
import android.os.Bundle;
import android.util.Log;
import android.view.Menu;
import android.view.View;
import android.widget.Button;
import android.widget.TextView;

public class MainActivity extends Activity {

	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_main);
		
		final TextView tvOutput = (TextView) findViewById(R.id.textView1);

		
		// Wire up the Go button 
		final Button okButton = (Button) findViewById(R.id.button1);
		okButton.setOnClickListener(new View.OnClickListener() {
			public void onClick(View v) {
				Log.d("++++","Clicked the GO button");

	        	// Get time in Zulu
	        	//final Calendar now = new Calendar.getInstance() ;
				final Date now = new Date();
				SimpleDateFormat dateFormatGmt = new SimpleDateFormat("HHmm");
	        	dateFormatGmt.setTimeZone(TimeZone.getTimeZone("GMT"));
	        	String timeZulu = dateFormatGmt.format(now);
	        	SimpleDateFormat hoursGmt = new SimpleDateFormat("HH");
	        	hoursGmt.setTimeZone(TimeZone.getTimeZone("GMT"));
	        	int zuluHour = Integer.parseInt(hoursGmt.format(now));
	        	Log.d("++++","Current GMT hour is " + zuluHour);
	        	
	        	SimpleDateFormat timeUTC = new SimpleDateFormat("HH:mm:ss MM/dd/yy");
	        	timeUTC.setTimeZone(TimeZone.getTimeZone("GMT"));
	        	String strtimeUTC = timeUTC.format(now);
	        	Log.d("++++","Current UTC time is " + strtimeUTC);
	        	
	        	
	        	Calendar nowCal =  Calendar.getInstance(TimeZone.getTimeZone("GMT"));
	        	nowCal.setTime(now);
				MoonPosition mMoonPosition = new MoonPosition(nowCal);
				mMoonPosition.calcalateMoonAzEl(nowCal);
				//MoonPhaseFinder mMoonPhaseFinder = new MoonPhaseFinder();

				SunPosition mSunPosition = new SunPosition(nowCal);
				mSunPosition.calcalateSunAzEl(nowCal);
				
				String moonLongitude = new DecimalFormat("###").format( mMoonPosition.getTrueLongitude()) ;
				String sunLongitude =  new DecimalFormat("###").format( mSunPosition.getEclipticLongitude()) ;
				//String moonElevation= String.valueOf(mMoonPosition.());
				//String sunLongitude = String.valueOf(mSunPosition.getDeclination());
				String moonPhase = new DecimalFormat("###").format( MoonPhaseFinder.getMoonAngle(nowCal));
				
				
				tvOutput.setText(strtimeUTC + "\n" 
						+ "Moon Longitude: " + moonLongitude + "\u00b0\n"
						+ "Sun Longitude: " + sunLongitude + "\u00b0\n"
						+ "Moon Phase Angle, 0-360\u00b0 (180\u00b0=full): " + moonPhase  + "\u00b0\n"
						+ "Phase Desc: " + MoonPhaseFinder.getMoonPhaseName(nowCal) + "\n"
						+ "Moon Altitude: " + new DecimalFormat("##0.0").format( mMoonPosition.getMoonElevationApparent()) + "\u00b0\n"
						+ "Moon Azimuth: " + new DecimalFormat("##0.0").format( mMoonPosition.getMoonAzimuth()) + "\u00b0\n"
						+ "Sun Az: " + new DecimalFormat("##0.00").format( mSunPosition.getSunAzimuth()) + "\u00b0\n"
						+ "Sun El: " + new DecimalFormat("##0.00").format( mSunPosition.getSunElevation()) + "\u00b0\n"
						+ "Distance to Moon: " + new DecimalFormat("#,###,###").format( mMoonPosition.getMoonDelta()) + " km\n");
			
			
			} 
		});
		
	}

	@Override
	public boolean onCreateOptionsMenu(Menu menu) {
		// Inflate the menu; this adds items to the action bar if it is present.
		getMenuInflater().inflate(R.menu.main, menu);
		return true;
	}
	
	private double[] getGPS() {
		 LocationManager lm = (LocationManager) getSystemService(
		  Context.LOCATION_SERVICE);
		 List<String> providers = lm.getProviders(true);

		 Location l = null;

		 for (int i=providers.size()-1; i>=0; i--) {
		  l = lm.getLastKnownLocation(providers.get(i));
		  if (l != null) break;
		 }

		 double[] gps = new double[2];
		 if (l != null) {
		  gps[0] = l.getLatitude();
		  gps[1] = l.getLongitude();
		 }

		 return gps;
		}

}
