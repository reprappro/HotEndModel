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
 * This is a simple program to model the behaviour of a RepRap hot end (or 
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
 * For constant p, this ODE has the solution:
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

	
	double dt = 0.2;        // Time increment for Euler integration
	double endTime = 120;    // How long to run the simulation (secs)
	double T0 = 25;         // Room temperature
	double a = 0.15;        // Hot end parameter - see above
	double b = 0.15;        // Hot end parameter - see above
	double lag = 2;         // Time lag between power change and effect (secs)
	
	// PID variables
	
	double e = 0;
	double eLast = 0;
	double eIntegral = 0;
	double pTest = 0.4;
	double Kp = 0.005;
	double Ki = 0.014;	
	double Kd = 0.005;
	double clamp = 50;
	double target = 205;
	double range = 255;
	
	String dataFile = "/home/ensab/Pro/Git/Marlin/Slave/xx.pde";
	
	double experiment[];
	
	class Ring
	{
		int l;
		double buffer[];
		int in;
		int out;
		
		Ring(double v)
		{
			l = (int) (1 + lag/dt);
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
		int i = 0;
		try 
		{
		    BufferedReader in = new BufferedReader(new FileReader(dataFile));
		    String str;
		    while ((str = in.readLine()) != null)
		    {
		    	if(!str.isEmpty())
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

	}
	
    public void buildChart(XYSeries dataA, XYSeries dataB, XYSeries dataC) throws Exception 
    {  
       XYSeriesCollection dataset = new XYSeriesCollection();
       if(dataA != null) dataset.addSeries( dataA );  
       if(dataB != null) dataset.addSeries( dataB );
       if(dataC != null) dataset.addSeries( dataC );
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
       plot.setRenderer(renderer); 
       
       ChartFrame chartFrame = new ChartFrame("XYArea Chart", chart);
       chartFrame.setVisible(true);
       chartFrame.setSize(800, 600);
    }  
	
	public double nextT(double T, double p, Ring r)
	{
		double result = T + (a*r.get()*range - b*(T-T0))*dt;
		r.add(p);
		return result;
	}
	
//	public XYSeries constantPower()
//	{
//		XYSeries result = new XYSeries( "Constant power: " + (int)(pTest*range) );
//		double t = 0;
//		double T = T0;
//		Ring r = new Ring(0);
//		while(t < endTime)
//		{
//			result.add(t, T); 
//			T = nextT(T, pTest, r);
//			t += dt;
//		}
//		return result;
//	}
	
	public double PID(double T)
	{
		eLast = e;
		e = target - T;
		double result = Kp*e + Ki*eIntegral + Kd*(e - eLast)/dt;
		eIntegral += e*dt;
		if(eIntegral < -clamp) eIntegral = -clamp;
		if(eIntegral > clamp) eIntegral = clamp;
		return result;
	}
	
	public XYSeries PIDPower()
	{
		XYSeries result = new XYSeries( "PID power");
		double t = 0;
		double T = T0;
		Ring r = new Ring(0);
		while(t < endTime)
		{
			result.add(t, T);
			double p = PID(T);
			T = nextT(T, p, r);
			t += dt;
		}
		return result;
	}
	
	public double odeSolution(double t)
	{
		if(t <= lag) return T0;
		return T0 + a*pTest*range*( 1 - Math.exp(b*(lag - t)) )/b;
	}
	
	public XYSeries analyticalSolution()
	{
		XYSeries result = new XYSeries( "Theory");
		double t = 0;
		while(t < endTime)
		{
			result.add(t, odeSolution(t));
			t += dt;
		}
		return result;
	}
	
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
	
	public static void main(String[] args) 
	{
	  HotEndModel hem = new HotEndModel();
      XYSeries pid = hem.PIDPower();
      XYSeries analytical = hem.analyticalSolution();
      XYSeries expt = hem.experimentalObservations();
		try {
			hem.buildChart(expt, analytical, pid);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}  
	}
}
