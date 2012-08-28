import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileReader;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.DeviationRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * This is a simple Java program to model the behaviour of a RepRap hot end (or 
 * anything else that gets heated) to allow the parameters of a PID for controlling
 * it to be set.
 * 
 * The Hot End is modelled as follows:
 * 
 *   dT/dt = a*p + b*(T - T0)    ----- (1)
 *   
 *   t is time in seconds
 *   T is temperature in C
 *   p is the dimensionless input power.
 *   T0 is room temperature
 *   a and b are parameters that get set experimentally (see below).
 *   
 * For constant p, ODE (1) has the solution:
 *   
 *   T = T0 + a*p*[1 - e^(b(tLag - t))]/b  -----(2)
 *   
 * where tLag is the time between applying power and any change in temperature occurring.
 *   
 * p is changed to a PWM value in [0, range] internally.  
 * 
 * Generally, the unknown values are a, b and tLag.
 * 
 * So the procedure is to fix a value of p, run the hot end up from
 * cold using it, and record the temperatures every dt seconds into
 * a data file.
 * 
 * This program then reads that data file (the first temperature in which is taken to be T0)
 * and compares it with Equation (2).  Arbitrary, but sensible, values of a, b and tLag are set
 * for the initial comparison.  The program works out a residual sum of squares for the errors, 
 * then iterates to find the values of a, b and tLag that minimise that RSS.  The iteration 
 * computes numerical values for the partial derivatives of the RSS with respect to a, b and tLag
 * then uses a downhill pathfinder to get to a minimum RSS.
 * 
 * @author Adrian Bowyer
 * RepRapPro Ltd
 * http://reprappro.com
 * 
 * Licence: GPL
 *
 */

public class HotEndModel 
{

	
	double dt = 0.4;        // Time increment for Euler integration and experiment
	double endTime = 240;   // How long to run the simulation (secs)
	double T0 = 25;         // Room temperature
	double gradEst = 0.01;  // Used in numerical estimation of derivatives
	double pTest = 0.4;		// Power used in the experimental test
	double range = 255;		// Actual power values to multiply p by
	
	// Hot end parameters
	
	static final int a = 0;
	static final int b = 1;
	static final int lag = 2;
	double aDefault = 0.05;      // Hot end parameter - see above
	double bDefault = 0.02;      // Hot end parameter - see above
	double lagDefault = 5;       // Time lag between power change and effect (secs)	
	
	// PID variables
	
	static final int Kp = 0;
	static final int Ki = 1;
	static final int Kd = 2;
	static final int clamp = 3;
	double e = 0;
	double eLast = 0;
	double eIntegral = 0;
	double target = 205;
	double heatingTime = 40;
	double KpDefault = 0.0015;
	double KiDefault = 0.005;	
	double KdDefault = 0.02;
	double clampDefault = 70;		

	String experimentalData = "t04";
	
	double experiment[];
	double ideal[];
	
	/**
	 * Ring buffer to hold power values to achieve time lag
	 * @author ensab
	 *
	 */
	class Ring
	{
		int l;
		double buffer[];
		int in;
		int out;
		
		Ring(double v, double lagTime)
		{
			l = (int) (1 + lagTime/dt);
			buffer = new double[l];
			in = 0;
			out = l;
			for(int i = 0; i < l; i++) buffer[i] = v;
		}
		
		void add(double v)
		{
			buffer[in] = v;
			in++;
			if (in >= l) in = 0;
		}
		
		double get()
		{
			out++;
			if (out >= l) out = 0;
			return buffer[out];
		}
	}
	
	public HotEndModel()
	{
		experiment = new double[5 + (int)(endTime/dt)];
		ideal = new double[5 + (int)(endTime/dt)];
		int i = 0;
		try 
		{
		    BufferedReader in = new BufferedReader(new FileReader(experimentalData));
		    String str;
		    while ((str = in.readLine()) != null)
		    {
		    	if(!str.isEmpty()) // Data files sometimes have extra newlines
		    	{
		    		experiment[i] = Double.parseDouble(str);
		    		i++;
		    	}
		    }
		    in.close();
		} catch (Exception e) 
		{
			e.printStackTrace();
		}
		T0 = experiment[0];
	}

	/**
	 * Set up the ideal heating profile
	 * @param pr
	 */
	public void setIdeal()
	{
		double t = 0;
		int i = 0;
		while(t < endTime)
		{
			if(t < heatingTime)
			{
				ideal[i] = T0 + (target - T0)*t/heatingTime;
			} else
				ideal[i] = target;
			i++;
			t += dt;
		}
	}
	
    public void buildChart(XYSeries dataA, XYSeries dataB, XYSeries dataC, XYSeries dataD) throws Exception 
    {  
       XYSeriesCollection dataset = new XYSeriesCollection();
       if(dataA != null) dataset.addSeries( dataA );  
       if(dataB != null) dataset.addSeries( dataB );
       if(dataC != null) dataset.addSeries( dataC );
       if(dataD != null) dataset.addSeries( dataD );
       JFreeChart chart = ChartFactory.createXYLineChart( "Hot-end graph",  
                                                          "Time (secs)",  
                                                          "Temp (C)",  
                                                          dataset,  
                                                          PlotOrientation.VERTICAL,  
                                                          true,  
                                                          false,  
                                                          false );
       XYPlot plot = (XYPlot) chart.getPlot();
//       plot.setBackgroundPaint(Color.lightGray);
//       plot.setDomainGridlinePaint(Color.white);
//       plot.setRangeGridlinePaint(Color.white);
       DeviationRenderer renderer = new DeviationRenderer(true, false);
       renderer.setSeriesPaint(0, Color.red);
       renderer.setSeriesPaint(1, Color.blue);
       renderer.setSeriesPaint(2, Color.yellow);
       renderer.setSeriesPaint(3, Color.black);
       plot.setRenderer(renderer); 
       
       ChartFrame chartFrame = new ChartFrame("XYArea Chart", chart);
       chartFrame.setVisible(true);
       chartFrame.setSize(800, 600);
    }  
	
    /**
     * Take one Euler step in the ODE
     * TODO: Change to Runge-Kutta?
     * @param T
     * @param p
     * @param hotEnd
     * @param r
     * @return
     */
	public double nextT(double T, double p, double[] hotEnd, Ring r)
	{
		double result = T + ( hotEnd[a]*r.get()*range - hotEnd[b]*(T-T0) )*dt;
		r.add(p);
		return result;
	}
	
//	public XYSeries constantPower(Parameters parm)
//	{
//		XYSeries result = new XYSeries( "Constant power: " + (int)(pTest*range) );
//		double t = 0;
//		double T = T0;
//		Ring r = new Ring(0, parm);
//		while(t < endTime)
//		{
//			result.add(t, T); 
//			T = nextT(T, pTest, parm, r);
//			t += dt;
//		}
//		return result;
//	}
	
	/**
	 * Return the power demanded by the PID
	 */
	public double PID(double T, double[] pid)
	{
		eLast = e;
		e = target - T;
		double result = pid[Kp]*e + pid[Ki]*eIntegral + pid[Kd]*(e - eLast)/dt;
		eIntegral += e*dt;
		if(eIntegral < -pid[clamp]) eIntegral = -pid[clamp];
		if(eIntegral > pid[clamp]) eIntegral = pid[clamp];
		if(result < 0) result = 0;
		if(result > 1) result = 1;
		return result;
	}
	
	/**
	 * Calculate the PID response
	 * @param pid
	 * @return
	 */
	public double[] pidPredict(double[] pid, double[] hotEnd)
	{
		double result[] = new double[experiment.length];
		double t = 0;
		double T = T0;
		int i = 0;
		Ring r = new Ring(0, hotEnd[lag]);
		while(t < endTime)
		{
			result[i] = T;
			double p = PID(T, pid);
			T = nextT(T, p, hotEnd, r);
			t += dt;
			i++;
		}
		return result;
	}
	
	/**
	 * Plot curve for the PID response
	 * @param pid
	 * @return
	 */
	public XYSeries PIDPower(double[] pid, double[] hotEnd)
	{
		XYSeries result = new XYSeries("PID");
		double r[] = pidPredict(pid, hotEnd);
		int i = 0;
		double t = 0;
		while(t < endTime)
		{
			result.add(t, r[i]);
			t += dt;
			i++;
		}
		return result;
	}
	
	/**
	 * Return the analytical solution at time t for constant power pTest
	 * @param t
	 * @param hotEnd
	 * @return
	 */
	public double odeSolution(double t, double[] hotEnd)
	{
		if(t <= hotEnd[lag]) return T0;
		return T0 + hotEnd[a]*pTest*range*( 1 - Math.exp(hotEnd[b]*(hotEnd[lag] - t)) )/hotEnd[b];
	}
	
	/**
	 * Plot curve for the analytical solution at constant power pTest
	 * @param hotEnd
	 * @return
	 */
	public XYSeries analyticalSolution(double[] hotEnd)
	{
		XYSeries result = new XYSeries( "Theory");
		double t = 0;
		while(t < endTime)
		{
			result.add(t, odeSolution(t, hotEnd));
			t += dt;
		}
		return result;
	}
	
	/**
	 * Compute the residual sum of squares between the analytical solution for
	 * parm and the experiment.
	 * @param hotEnd
	 * @return
	 */
	public double heRSS(double[] hotEnd)
	{
		double r = 0;
		double t = 0;
		int i = 0;
		while(t < endTime)
		{
			double d = odeSolution(t, hotEnd) - experiment[i];
			r += d*d;
			t += dt;
			i++;
		}
		return r;
	}
	
	/**
	 * Compute the residual sum of squares between the PID 
	 * and the ideal.
	 * @param pid
	 * @return
	 */
	public double pidRSS(double[] pid, double[] hotEnd)
	{
		double v[] = pidPredict(pid, hotEnd);
		double r = 0;
		double t = 0;
		int i = 0;
		while(t < endTime)
		{
			double d = v[i] - ideal[i];
			r += d*d;
			t += dt;
			i++;
		}
		return r;
	}
	
	/**
	 * Compute the partial derivative of the experiment RSS w.r.t. parm[i]
	 * @param r
	 * @param hotEnd
	 * @param i
	 * @return
	 */
	public double heDrByDp(double r, double[] hotEnd, int i)
	{
		double delta = Math.abs(hotEnd[i]*gradEst);
		double[] temp = new double[hotEnd.length];
		for(int j = 0; j < hotEnd.length; j++) temp[j] = hotEnd[j];
		temp[i] = hotEnd[i] + delta;
		double newR = heRSS(temp);
		return (newR - r)/delta;
	}
	
	/**
	 * Compute the partial derivative of the PID RSS w.r.t. parm[i]
	 * @param r
	 * @param pid
	 * @param i
	 * @return
	 */
	public double pidDrByDp(double r, double[] pid, double[] hotEnd, int i)
	{
		double delta = Math.abs(pid[i]*gradEst);
		double[] temp = new double[pid.length];
		for(int j = 0; j < pid.length; j++) temp[j] = pid[j];
		temp[i] = pid[i] + delta;
		double newR = pidRSS(temp, hotEnd);
		return (newR - r)/delta;
	}

	
	/**
	 * Dumb descend-the-slope optimiser to fit the parameters to the
	 * experiment.
	 * @return
	 */
	public double[] heFit()
	{
		double[] hotEnd = {aDefault, bDefault, lagDefault};
		double r = heRSS(hotEnd);
		double[] newHotEnd;
		double newR;
		int count = 0;
		double g;
		while(count < 200 && Math.sqrt(r)*dt/endTime > 0.1)
		{
			for(int i = 0; i < hotEnd.length; i++)
			{
				g = heDrByDp(r, hotEnd, i);
				newHotEnd = new double[hotEnd.length];
				for(int j = 0; j < hotEnd.length; j++) newHotEnd[j] = hotEnd[j];
				if(g > 0)
					newHotEnd[i] = hotEnd[i] - Math.abs(hotEnd[i]*gradEst);
				else
					newHotEnd[i] = hotEnd[i] + Math.abs(hotEnd[i]*gradEst);
				newR = heRSS(newHotEnd);
				if(newR < r)
				{
					r = newR;
					hotEnd = newHotEnd;
				}
			}
			count++;
		}
		System.out.println(" Hot end RMS error: " + Math.sqrt(r)*dt/endTime + ", iterations: " + count);
		System.out.println(" Hot end parameters: a = " + hotEnd[a] + ", b = " + hotEnd[b] + ", lag = " + hotEnd[lag]);
		return hotEnd;
	}
	
	/**
	 * Dumb descend-the-slope optimiser to fit the pid to the
	 * ideal.
	 * @return
	 */
	public double[] pidFit(double[] hotEnd)
	{
		double[] pid = { KpDefault, KiDefault, KdDefault, clampDefault };
		double r = pidRSS(pid, hotEnd);
		double[] newPid;
		double newR;
		int count = 0;
		double g;
		while(count < 200 && Math.sqrt(r)*dt/endTime > 0.1)
		{
			for(int i = 0; i < pid.length; i++)
			{
				g = pidDrByDp(r, pid,  hotEnd, i);
				newPid = new double[pid.length];
				for(int j = 0; j < pid.length; j++) newPid[j] = pid[j];
				if(g > 0)
					newPid[i] = pid[i] - Math.abs(pid[i]*gradEst);
				else
					newPid[i] = pid[i] + Math.abs(pid[i]*gradEst);
				newR = pidRSS(newPid, hotEnd);
				if(newR < r)
				{
					r = newR;
					pid = newPid;
				}
			}
			count++;
		}
		System.out.println(" PID RMS error: " + Math.sqrt(r)*dt/endTime + ", iterations: " + count);
		System.out.println(" PID parameters: Kp = " + pid[Kp] + ", Ki = " + pid[Ki] + ", Kd = " + pid[Kd]  + ", clamp = " + pid[clamp]);
		return pid;
	}
	
	/**
	 * Plot curve for the experimental results
	 * @return
	 */
	public XYSeries experimentalObservations()
	{
		XYSeries result = new XYSeries( "Experiment");
		double t = 0;
		int i = 0;
		while(t < endTime)
		{
			result.add(t, experiment[i]);
			t += dt;
			i++;
		}
		return result;
	}
	
	/**
	 * Plot curve for the ideal profile
	 * @return
	 */
	public XYSeries idealProfile()
	{
		XYSeries result = new XYSeries( "Ideal");
		double t = 0;
		int i = 0;
		while(t < endTime)
		{
			result.add(t, ideal[i]);
			t += dt;
			i++;
		}
		return result;
	}
	
	public static void main(String[] args) 
	{
		HotEndModel hem = new HotEndModel();
		double[] hotEnd = hem.heFit();
		hem.setIdeal();
		double[] pid = { hem.KpDefault, hem.KiDefault, hem.KdDefault, hem.clampDefault };
		XYSeries pd = hem.PIDPower(pid, hotEnd);
		XYSeries analytical = hem.analyticalSolution(hotEnd);
		XYSeries expt = hem.experimentalObservations();
		XYSeries idl = hem.idealProfile();
		try 
		{
			hem.buildChart(expt, analytical, idl, pd);
		} catch (Exception e) 
		{
			e.printStackTrace();
		}  
	}
}
