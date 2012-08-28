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

	
	double dt = 0.4;        // Time increment for Euler integration
	double endTime = 240;   // How long to run the simulation (secs)
	double T0 = 25;         // Room temperature
	double gradEst = 0.01;  // Used in numerical estimation of derivatives
	
	// PID variables
	
	double e = 0;
	double eLast = 0;
	double eIntegral = 0;
	double pTest = 0.4;
	double Kp = 0.001;
	double Ki = 0.0014;	
	double Kd = 0;
	double clamp = 50;
	double target = 205;
	double range = 255;
	
	String dataFile = "t04";
	
	double experiment[];
	
	class Parameters
	{
		// Initial guess values
		
		double a = 0.05;      // Hot end parameter - see above
		double b = 0.02;      // Hot end parameter - see above
		double lag = 5;       // Time lag between power change and effect (secs)
		
		Parameters() {}
		
		Parameters(double aa, double bb, double ll)
		{
			a = aa;
			b = bb;
			lag = ll;
		}
		
		public String toString()
		{
			String r = "a: " + a + ", b: " + b + ", lag: " + lag;
			return r;
		}
	}
	
	class Ring
	{
		int l;
		double buffer[];
		int in;
		int out;
		
		Ring(double v, Parameters parm)
		{
			l = (int) (1 + parm.lag/dt);
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
		T0 = experiment[0];
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
	
    /**
     * Take one Euler step in the ODE
     * TODO: Change to Runge-Kutta?
     * @param T
     * @param p
     * @param parm
     * @param r
     * @return
     */
	public double nextT(double T, double p, Parameters parm, Ring r)
	{
		double result = T + (parm.a*r.get()*range - parm.b*(T-T0))*dt;
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
	
	/**
	 * Plot curve for the PID response
	 * @param parm
	 * @return
	 */
	public XYSeries PIDPower(Parameters parm)
	{
		XYSeries result = new XYSeries( "PID power");
		double t = 0;
		double T = T0;
		Ring r = new Ring(0, parm);
		while(t < endTime)
		{
			result.add(t, T);
			double p = PID(T);
			T = nextT(T, p, parm, r);
			t += dt;
		}
		return result;
	}
	
	/**
	 * Return the analytical solution at time t for constant power pTest
	 * @param t
	 * @param parm
	 * @return
	 */
	public double odeSolution(double t, Parameters parm)
	{
		if(t <= parm.lag) return T0;
		return T0 + parm.a*pTest*range*( 1 - Math.exp(parm.b*(parm.lag - t)) )/parm.b;
	}
	
	/**
	 * Plot curve for the analytical solution at constant power pTest
	 * @param parm
	 * @return
	 */
	public XYSeries analyticalSolution(Parameters parm)
	{
		XYSeries result = new XYSeries( "Theory");
		double t = 0;
		while(t < endTime)
		{
			result.add(t, odeSolution(t, parm));
			t += dt;
		}
		return result;
	}
	
	/**
	 * Compute the residual sum of squares between the analytical solution for
	 * parm and the experiment.
	 * @param parm
	 * @return
	 */
	public double rss(Parameters parm)
	{
		double r = 0;
		double t = 0;
		int i = 0;
		while(t < endTime)
		{
			double d = odeSolution(t, parm) - experiment[i];
			r += d*d;
			t += dt;
			i++;
		}
		return r;
	}
	
	/**
	 * Compute the partial derivative of the RSS w.r.t. parm.a
	 * @param r
	 * @param parm
	 * @return
	 */
	public double drByDa(double r, Parameters parm)
	{
		double delta = Math.abs(parm.a*gradEst);
		Parameters temp = new Parameters(parm.a + delta, parm.b, parm.lag);
		double newR = rss(temp);
		return (newR - r)/delta;
	}
	
	/**
	 * Compute the partial derivative of the RSS w.r.t. parm.b
	 * @param r
	 * @param parm
	 * @return
	 */
	public double drByDb(double r, Parameters parm)
	{
		double delta = Math.abs(parm.b*gradEst);
		Parameters temp = new Parameters(parm.a, parm.b + delta, parm.lag);
		double newR = rss(temp);
		return (newR - r)/delta;
	}
	
	/**
	 * Compute the partial derivative of the RSS w.r.t. parm.lag
	 * @param r
	 * @param parm
	 * @return
	 */
	public double drByDlag(double r, Parameters parm)
	{
		double delta = Math.abs(parm.lag*gradEst);
		Parameters temp = new Parameters(parm.a, parm.b, parm.lag + delta);
		double newR = rss(temp);
		return (newR - r)/delta;
	}
	
	/**
	 * Dumb descend-the-slope optimiser to fit the parameters to the
	 * experiment.
	 * @return
	 */
	public Parameters fit()
	{
		Parameters pr = new Parameters();
		double r = rss(pr);
		double newR;
		Parameters newPr;
		int count = 0;
		double g;
		while(count < 200 && Math.sqrt(r)*dt/endTime > 0.1)
		{
			g = drByDa(r, pr);
			if(g > 0)
				newPr = new Parameters(pr.a - Math.abs(pr.a*gradEst), pr.b, pr.lag);
			else
				newPr = new Parameters(pr.a + Math.abs(pr.a*gradEst), pr.b, pr.lag);
			newR = rss(newPr);
			if(newR < r)
			{
				r = newR;
				pr = newPr;
			}
			
			g = drByDb(r, pr);
			if(g > 0)
				newPr = new Parameters(pr.a, pr.b - Math.abs(pr.b*gradEst), pr.lag);
			else
				newPr = new Parameters(pr.a, pr.b + Math.abs(pr.b*gradEst), pr.lag);
			newR = rss(newPr);
			if(newR < r)
			{
				r = newR;
				pr = newPr;
			}
			
			g = drByDlag(r, pr);
			if(g > 0)
				newPr = new Parameters(pr.a, pr.b, pr.lag - Math.abs(pr.lag*gradEst));
			else
				newPr = new Parameters(pr.a, pr.b, pr.lag + Math.abs(pr.lag*gradEst));
			newR = rss(newPr);
			if(newR < r)
			{
				r = newR;
				pr = newPr;
			}
			
			count++;
		}
		System.out.println("RMS error: " + Math.sqrt(r)*dt/endTime + ", iterations: " + count);
		System.out.println("Parameters: " + pr);
		return pr;
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
	
	public static void main(String[] args) 
	{
		HotEndModel hem = new HotEndModel();
		Parameters pr = hem.fit();
		XYSeries pid = hem.PIDPower(pr);
		XYSeries analytical = hem.analyticalSolution(pr);
		XYSeries expt = hem.experimentalObservations();
		try 
		{
			hem.buildChart(expt, analytical, pid);
		} catch (Exception e) 
		{
			e.printStackTrace();
		}  
	}
}
