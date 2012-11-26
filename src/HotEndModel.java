import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.FileReader;

import javax.swing.BorderFactory;
//import javax.swing.ButtonGroup;
import javax.swing.JButton;
//import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
//import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

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
 * anything else that gets heated) to allow the parameters of a PID for
 * controlling it to be set.
 * 
 * The Hot End is modelled as follows:
 * 
 * dT/dt = a*p + b*(T - T0) ----- (1)
 * 
 * t is time in seconds T is temperature in C p is the dimensionless input
 * power. T0 is room temperature a and b are parameters that get set
 * experimentally (see below).
 * 
 * For constant p, ODE (1) has the solution:
 * 
 * T = T0 + a*p*[1 - e^(b(tLag - t))]/b -----(2)
 * 
 * where tLag is the time between applying power and any change in temperature
 * occurring.
 * 
 * p is changed to a PWM value in [0, range] internally.
 * 
 * Generally, the unknown values are a, b and tLag.
 * 
 * So the procedure is to fix a value of p, run the hot end up from cold using
 * it, and record the temperatures every dt seconds into a data file.
 * 
 * Referring to the graph, you run a simple experiment setting the PWM
 * into the heater constant (around 0.4*255-ish) and logging the results every
 * 0.4s (or whatever). This gives the red curve.
 * 
 * The program then does an iterative least-squares fit (the blue curve) of the
 * analytical solution of the differential equation that describes the hot end
 * to the experimental results. This gives all the hot end parameters (response
 * to input power, cooling, and time lag) automatically.
 * 
 * You then specify your ideal heating characteristic (the yellow curve).
 * 
 * The program then takes the hot end parameters it found plus the PID
 * controller terms and does another least-squares fit, this time keeping the
 * found hot-end parameters constant and optimising the PID parameters to get
 * the response as close as possible to the yellow curve - the black curve.
 * 
 * Then the final step (which I hope to do later today...): you put those PID
 * parameters into your firmware and see if it really does do what the model
 * thinks...
 * 
 * @author Adrian Bowyer RepRapPro Ltd http://reprappro.com
 * 
 *         Licence: GPL
 * 
 */


/**
 * Class to put a selection panel up for the user to type numbers into
 */
class UserInteraction extends JPanel 
{
	static String experimentalData = null; // The name of the file with the experimental data

	private HotEndModel hem;
	
	private double dt = 1.0;      // Time increment for Euler integration and experiment
	private double endTime = 240; // How long to run the simulation (secs)
	//private double T0 = 25;       // Room temperature
	//private double gradEst = 0.01;// Used in numerical estimation of derivatives
	private double pTest = 0.25;   // Power used in the experimental test
	private double range = 255;   // Actual power values to multiply p by

	// Hot end parameters

	private double aDefault = 0.05;   // Guessed hot end parameter - see above
	private double bDefault = 0.02;   // Guessed hot end parameter - see above
	private double lagDefault = 5;    // Guessed time lag between power change and effect (secs)

	// PID variables

	private double target = 205;      // Target hot end temperature
	private double heatingTime = 40;  // Time to get to temp ideally
	private double KpDefault = 0.0015;// Guessed PID parameter
	private double KiDefault = 0.005; // Guessed PID parameter
	private double KdDefault = 0.02;  // Guessed PID parameter
	private double clampDefault = 70; // Guessed PID parameter (The integral clamps to +/- this value)


	
	private final long serialVersionUID = 1L;
	private JDialog dialog;
	private JTextField dtBox;
	private JTextField endTimeBox;          // How long to run the simulation (secs)

	private JTextField pTestBox;            // Power used in the experimental test
	private JTextField rangeBox;            // Actual power values to multiply p by

//	Hot end parameters
	
	private JTextField aDefaultBox;         // Guessed hot end parameter - see above
	private JTextField bDefaultBox;         // Guessed hot end parameter - see above
	private JTextField lagDefaultBox;       // Guessed time lag between power change and effect (secs)

//	PID variables

	private JTextField targetBox;           // Target hot end temperature
	private JTextField heatingTimeBox;      // Time to get to temp ideally
	private JTextField KpDefaultBox;        // Guessed PID parameter
	private JTextField KiDefaultBox;        // Guessed PID parameter
	private JTextField KdDefaultBox;        // Guessed PID parameter
	private JTextField clampDefaultBox;     // Guessed PID parameter (The integral clamps to +/- this value)
	//private JTextField experimentalDataBox; // File name
	private JPanel radioPanel;
	
	/**
	 * Set my parameters - called by the hot end model class to copy its values into
	 * the GUI text boxes.
	 * 
	 * @param dtp
	 * @param endTimep
	 * @param pTestp
	 * @param rangep
	 * @param aDefaultp
	 * @param bDefaultp
	 * @param lagDefaultp
	 * @param targetp
	 * @param heatingTimep
	 * @param KpDefaultp
	 * @param KiDefaultp
	 * @param KdDefaultp
	 * @param clampDefaultp
	 * @param experimentalDatap
	 */
	public void setParameters(double dtp, double endTimep, double pTestp, double rangep, double aDefaultp, 
			double bDefaultp, double lagDefaultp, double targetp, double heatingTimep, double KpDefaultp, 
			double KiDefaultp, double KdDefaultp, double clampDefaultp, String experimentalDatap)
	{
		dt = dtp; 
		endTime = endTimep;
		pTest = pTestp;
		range = rangep;
		aDefault = aDefaultp; 
		bDefault = bDefaultp;
		lagDefault = lagDefaultp; 
		target = targetp;
		heatingTime = heatingTimep; 
		KpDefault = KpDefaultp;
		KiDefault = KiDefaultp;
		KdDefault = KdDefaultp;
		clampDefault = clampDefaultp; 
		//experimentalData = experimentalDatap;
	}
	
	/**
	 * Setup a text box with a double in it
	 * @param name
	 * @param val
	 * @return
	 */
	private JTextField setUp(String name, double val)
	{
		JLabel jLabel = new JLabel();
		radioPanel.add(jLabel);
		jLabel.setText(name);
		jLabel.setHorizontalAlignment(SwingConstants.LEFT);
		JTextField jtf = new JTextField(Double.toString(val));
		jtf.setSize(20, 10);
		jtf.setHorizontalAlignment(SwingConstants.LEFT);
		radioPanel.add(jtf);
		return jtf;
	}
	
	/**
	 * Set up and display the data panel
	 * @param h
	 * @param oderss
	 * @param pidrss
	 */
	public void run(HotEndModel h, double oderss, double pidrss)
	{
		JFrame f = new JFrame();
		dialog = new JDialog(f, "Set Parameters");
		//dialog.setLocation(500, 400);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

		hem = h;

		radioPanel = new JPanel(new GridLayout(0, 1));
		radioPanel.setSize(300,200);

		JLabel jLabel1 = new JLabel();
		radioPanel.add(jLabel1);
		jLabel1.setText("RepRapPro Ltd Hot End Modeller");
		jLabel1.setHorizontalAlignment(SwingConstants.CENTER);

		dtBox = setUp("time increment (s)", dt);
		endTimeBox = setUp("end time (s)", endTime);
		pTestBox = setUp("test power [0, 1]", pTest);
		rangeBox = setUp("D->A range", range);
		aDefaultBox = setUp("ODE a value", aDefault);
		bDefaultBox = setUp("ODE b value", bDefault);
		lagDefaultBox = setUp("ODE lag value (s)", lagDefault);
		targetBox = setUp("target temperature (C)", target);
		heatingTimeBox = setUp("heating time (s)", heatingTime);
		KpDefaultBox = setUp("PID Kp", KpDefault);
		KiDefaultBox = setUp("PID Ki", KiDefault);
		KdDefaultBox = setUp("PID Kd", KdDefault);
		clampDefaultBox = setUp("PID integral clamp", clampDefault);
		
//		JLabel jLabel2 = new JLabel();
//		radioPanel.add(jLabel2);
//		jLabel2.setText("Experiment file name");
//		jLabel2.setHorizontalAlignment(SwingConstants.LEFT);
//		experimentalDataBox = new JTextField(experimentalData);
//		experimentalDataBox.setSize(20, 10);
//		experimentalDataBox.setHorizontalAlignment(SwingConstants.RIGHT);
//		radioPanel.add(experimentalDataBox);
		
		JLabel jLabel3 = new JLabel();
		radioPanel.add(jLabel3);
		jLabel3.setText("ODE RMS error: " + oderss);
		jLabel3.setHorizontalAlignment(SwingConstants.LEFT);
		
		JLabel jLabel4 = new JLabel();
		radioPanel.add(jLabel4);
		jLabel4.setText("PID RMS error: " + pidrss);
		jLabel4.setHorizontalAlignment(SwingConstants.LEFT);
		
		if(experimentalData == null)
		{
			final JFileChooser fc = new JFileChooser();
			fc.setDialogTitle("Experimental heating file");
			fc.showOpenDialog(this);
			experimentalData = fc.getSelectedFile().toString();
		}
		
		
		try
		{

			JButton runButton = new JButton();
			radioPanel.add(runButton);
			runButton.setText("Run");
			runButton.addActionListener(new ActionListener() 
			{
				public void actionPerformed(ActionEvent evt) 
				{
					runHandler();
				}
			}
			);

			add(radioPanel, BorderLayout.LINE_START);
			setBorder(BorderFactory.createEmptyBorder(20,20,20,20));
			
			JButton exitButton = new JButton();
			radioPanel.add(exitButton);
			exitButton.setText("Exit");
			exitButton.addActionListener(new ActionListener() 
			{
				public void actionPerformed(ActionEvent evt) 
				{
					exitHandler();
				}
			}
			);

			add(radioPanel, BorderLayout.LINE_START);
			setBorder(BorderFactory.createEmptyBorder(20,20,20,20));

		} catch (Exception ex)
		{
			ex.printStackTrace();
		}	
		
		setOpaque(true); //content panes must be opaque
		dialog.setContentPane(this);

		//Display the window.
		dialog.pack();
		dialog.setModalityType(JDialog.DEFAULT_MODALITY_TYPE);
		dialog.setVisible(true);
	}
	
	/**
	 * Constructor does very little
	 * @param h
	 * @param oderss
	 * @param pidrss
	 */
	public UserInteraction()
	{
		super(new BorderLayout());
	} 

	public void runHandler()
	{
		hem.setParameters(
		Double.parseDouble(dtBox.getText()), 
		Double.parseDouble(endTimeBox.getText()),
		Double.parseDouble(pTestBox.getText()),
		Double.parseDouble(rangeBox.getText()),
		Double.parseDouble(aDefaultBox.getText()), 
		Double.parseDouble(bDefaultBox.getText()),
		Double.parseDouble(lagDefaultBox.getText()), 
		Double.parseDouble(targetBox.getText()),
		Double.parseDouble(heatingTimeBox.getText()), 
		Double.parseDouble(KpDefaultBox.getText()),
		Double.parseDouble(KiDefaultBox.getText()),
		Double.parseDouble(KdDefaultBox.getText()),
		Double.parseDouble(clampDefaultBox.getText()), 
		experimentalData
		);
		dialog.dispose();
	}
	
	public void exitHandler()
	{
		dialog.dispose();
		System.exit(0);
	}
}


public class HotEndModel {

	private double timePenalty = 100.0; // If this is 0, no penalty. Otherwise later errors weigh more than early ones
	private double dt = 1.0;      // Time increment for Euler integration and experiment
	private double endTime = 240; // How long to run the simulation (secs)
	private double T0 = 25;       // Room temperature
	private double gradEst = 0.01;// Used in numerical estimation of derivatives
	private double pTest = 0.25;   // Power used in the experimental test
	private double range = 255;   // Actual power values to multiply p by

	// Hot end parameters

	private double aDefault = 0.05;   // Guessed hot end parameter - see above
	private double bDefault = 0.02;   // Guessed hot end parameter - see above
	private double lagDefault = 5;    // Guessed time lag between power change and effect (secs)

	// PID variables

	private double target = 205;      // Target hot end temperature
	private double heatingTime = 40;  // Time to get to temp ideally
	private double KpDefault = 0.0015;// Guessed PID parameter
	private double KiDefault = 0.005; // Guessed PID parameter
	private double KdDefault = 0.02;  // Guessed PID parameter
	private double clampDefault = 70; // Guessed PID parameter (The integral clamps to +/- this value)

	private static String experimentalData = "t04"; // The name of the file with the experimental data

	// Hot end parameters

	static final int a = 0;   // Array locations of hot end parameters
	static final int b = 1;
	static final int lag = 2;

	// PID variables

	static final int Kp = 0;  // Array locations of PID parameters
	static final int Ki = 1;
	static final int Kd = 2;
	static final int clamp = 3;
	double e = 0;             // Error
	double eLast = 0;         // Last error
	double eIntegral = 0;     // Integral of error

	double experiment[];      // Array for the experimental data
	double ideal[];           // Array for the ideal heating profile
	
	double oderss = 0;        // Residual error from ODE fitting
	double pidrss = 0;        // Residual error from ODE fitting

	JDialog dialog;			  // Interaction with the user
	
	/**
	 * Ring buffer to hold power values to achieve time lag
	 * 
	 */
	class Ring {
		int l;
		double buffer[];
		int in;
		int out;

		Ring(double v, double lagTime) {
			l = (int) (1 + lagTime / dt);
			buffer = new double[l];
			in = 0;
			out = l;
			for (int i = 0; i < l; i++)
				buffer[i] = v;
		}

		void add(double v) {
			buffer[in] = v;
			in++;
			if (in >= l)
				in = 0;
		}

		double get() {
			out++;
			if (out >= l)
				out = 0;
			return buffer[out];
		}
	}
	
	/**
	 * Set parameter values.  Called by the GUI panel to load the data typed by the user
	 * into the local parameters.
	 * 
	 * @param dtp
	 * @param endTimep
	 * @param pTestp
	 * @param rangep
	 * @param aDefaultp
	 * @param bDefaultp
	 * @param lagDefaultp
	 * @param targetp
	 * @param heatingTimep
	 * @param KpDefaultp
	 * @param KiDefaultp
	 * @param KdDefaultp
	 * @param clampDefaultp
	 * @param experimentalDatap
	 */
	public void setParameters(double dtp, double endTimep, double pTestp, double rangep, double aDefaultp, 
			double bDefaultp, double lagDefaultp, double targetp, double heatingTimep, double KpDefaultp, 
			double KiDefaultp, double KdDefaultp, double clampDefaultp, String experimentalDatap)
	{
		dt = dtp; 
		endTime = endTimep;
		pTest = pTestp;
		range = rangep;
		aDefault = aDefaultp; 
		bDefault = bDefaultp;
		lagDefault = lagDefaultp; 
		target = targetp;
		heatingTime = heatingTimep; 
		KpDefault = KpDefaultp;
		KiDefault = KiDefaultp;
		KdDefault = KdDefaultp;
		clampDefault = clampDefaultp; 
		experimentalData = experimentalDatap;
	}
	
	/**
	 * Set everything up as at the start, but with whatever current parameter values
	 * there are.
	 */
	public void initiate()
	{
		UserInteraction ui = new UserInteraction();
		ui.setParameters(dt, endTime, pTest, range, aDefault, 
				bDefault, lagDefault, target, heatingTime, KpDefault, 
				KiDefault, KdDefault, clampDefault, experimentalData);
		ui.run(this, oderss, pidrss);
		
		experiment = new double[5 + (int) (endTime / dt)];
		ideal = new double[5 + (int) (endTime / dt)];
		int i = 0;
		try {
			BufferedReader in = new BufferedReader(new FileReader(
					experimentalData));
			String str;
			while ((str = in.readLine()) != null) {
				if (!str.isEmpty()) // Data files sometimes have extra newlines
				{
					experiment[i] = Double.parseDouble(str);
					i++;
				}
			}
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		T0 = experiment[0];		
	}

	
	/**
	 * Constructor reads in the experimental data and initialises the guesses
	 * from which to iterate.
	 */
	public HotEndModel() 
	{
		initiate();
	}

	/**
	 * Set up the ideal heating profile
	 * 
	 * @param pr
	 */
	public void setIdeal() {
		double t = 0;
		int i = 0;
		while (t < endTime) 
		{
			//if (t < heatingTime) 
			if 	(t < lagDefault)
			{
				ideal[i] = T0; // + (target - T0) * t / heatingTime;
			} else
				ideal[i] = target;
			i++;
			t += dt;
		}
	}

	/**
	 * Plot a graph of the results
	 * @param dataA
	 * @param dataB
	 * @param dataC
	 * @param dataD
	 * @throws Exception
	 */
	public void buildChart(XYSeries dataA, XYSeries dataB, XYSeries dataC,
			XYSeries dataD) throws Exception {
		XYSeriesCollection dataset = new XYSeriesCollection();
		if (dataA != null)
			dataset.addSeries(dataA);
		if (dataB != null)
			dataset.addSeries(dataB);
		if (dataC != null)
			dataset.addSeries(dataC);
		if (dataD != null)
			dataset.addSeries(dataD);
		JFreeChart chart = ChartFactory.createXYLineChart("Hot-end graph",
				"Time (secs)", "Temp (C)", dataset, PlotOrientation.VERTICAL,
				true, false, false);
		XYPlot plot = (XYPlot) chart.getPlot();
		// plot.setBackgroundPaint(Color.lightGray);
		// plot.setDomainGridlinePaint(Color.white);
		// plot.setRangeGridlinePaint(Color.white);
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
	 * Take one Euler step in the ODE TODO: Change to Runge-Kutta?
	 * 
	 * @param T
	 * @param p
	 * @param hotEnd
	 * @param r
	 * @return
	 */
	public double nextT(double T, double p, double[] hotEnd, Ring r) {
		double result = T
				+ (hotEnd[a] * r.get() * range - hotEnd[b] * (T - T0)) * dt;
		r.add(p);
		return result;
	}

	/**
	 * Return the power demanded by the PID
	 */
	public double PID(double T, double[] pid) 
	{
		eLast = e;
		e = target - T;
		double result = pid[Kp] * e + pid[Ki] * eIntegral + pid[Kd]
				* (e - eLast) / dt;
		eIntegral += e * dt;
		if (eIntegral < -pid[clamp]) eIntegral = -pid[clamp];
		if (eIntegral > pid[clamp])	eIntegral = pid[clamp];
		if (result < 0) result = 0;
		if (result > 1)	result = 1;
		return result;
	}

	/**
	 * Calculate the PID response given the hot end parameters
	 * and a guess at the PID parameters.
	 * 
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
		e = 0;
		eLast = 0;
		eIntegral = 0;
		while (t < endTime) {
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
	 * 
	 * @param pid
	 * @return
	 */
	public XYSeries PIDPower(double[] pid, double[] hotEnd) {
		XYSeries result = new XYSeries("PID");
		double r[] = pidPredict(pid, hotEnd);
		int i = 0;
		double t = 0;
		while (t < endTime) {
			result.add(t, r[i]);
			t += dt;
			i++;
		}
		return result;
	}

	/**
	 * Return the analytical ODE solution at time t for constant power pTest
	 * 
	 * @param t
	 * @param hotEnd
	 * @return
	 */
	public double odeSolution(double t, double[] hotEnd) {
		if (t <= hotEnd[lag])
			return T0;
		return T0 + hotEnd[a] * pTest * range
				* (1 - Math.exp(hotEnd[b] * (hotEnd[lag] - t))) / hotEnd[b];
	}

	/**
	 * Plot curve for the analytical solution at constant power pTest
	 * 
	 * @param hotEnd
	 * @return
	 */
	public XYSeries analyticalSolution(double[] hotEnd) {
		XYSeries result = new XYSeries("Theory");
		double t = 0;
		while (t < endTime) {
			result.add(t, odeSolution(t, hotEnd));
			t += dt;
		}
		return result;
	}

	/**
	 * Compute the residual sum of squares between the analytical solution for
	 * hotEnd and the experiment.
	 * 
	 * @param hotEnd
	 * @return
	 */
	public double heRSS(double[] hotEnd) {
		double r = 0;
		double t = 0;
		int i = 0;
		while (t < endTime) {
			double d = odeSolution(t, hotEnd) - experiment[i];
			r += d * d;
			t += dt;
			i++;
		}
		return r;
	}

	/**
	 * Compute the residual sum of squares between the PID and the ideal.
	 * 
	 * @param pid
	 * @return
	 */
	public double pidRSS(double[] pid, double[] hotEnd) {
		double v[] = pidPredict(pid, hotEnd);
		double r = 0;
		double t = 0;
		int i = 0;
		while (t < endTime) {
			double d = Math.abs(v[i] - ideal[i]);
			d += timePenalty*d*t/endTime;
			r += d * d;
			t += dt;
			i++;
		}
		return r;
	}

	/**
	 * Compute the partial derivative of the experiment RSS w.r.t. hotEnd[i]
	 * 
	 * @param r
	 * @param hotEnd
	 * @param i
	 * @return
	 */
	public double heDrByDp(double r, double[] hotEnd, int i) {
		double delta = Math.abs(hotEnd[i] * gradEst);
		double[] temp = new double[hotEnd.length];
		for (int j = 0; j < hotEnd.length; j++)
			temp[j] = hotEnd[j];
		temp[i] = hotEnd[i] + delta;
		double newR = heRSS(temp);
		return (newR - r) / delta;
	}

	/**
	 * Compute the partial derivative of the PID RSS w.r.t. pid[i]
	 * 
	 * @param r
	 * @param pid
	 * @param i
	 * @return
	 */
	public double pidDrByDp(double r, double[] pid, double[] hotEnd, int i) {
		double delta = Math.abs(pid[i] * gradEst);
		double[] temp = new double[pid.length];
		for (int j = 0; j < pid.length; j++)
			temp[j] = pid[j];
		temp[i] = pid[i] + delta;
		double newR = pidRSS(temp, hotEnd);
		return (newR - r) / delta;
	}

	/**
	 * Dumb descend-the-slope optimiser to fit the parameters to the experiment.
	 * 
	 * @return
	 */
	public double[] heFit() {
		double[] hotEnd = { aDefault, bDefault, lagDefault };
		double r = heRSS(hotEnd);
		double[] newHotEnd;
		double newR;
		int count = 0;
		double g;
		while (count < 200 && Math.sqrt(r) * dt / endTime > 0.1) {
			for (int i = 0; i < hotEnd.length; i++) {
				g = heDrByDp(r, hotEnd, i);
				newHotEnd = new double[hotEnd.length];
				for (int j = 0; j < hotEnd.length; j++)
					newHotEnd[j] = hotEnd[j];
				if (g > 0)
					newHotEnd[i] = hotEnd[i] - Math.abs(hotEnd[i] * gradEst);
				else
					newHotEnd[i] = hotEnd[i] + Math.abs(hotEnd[i] * gradEst);
				newR = heRSS(newHotEnd);
				if (newR < r) {
					r = newR;
					hotEnd = newHotEnd;
				}
			}
			count++;
		}
//		System.out.println(" Hot end RMS error: " + Math.sqrt(r) * dt / endTime
//				+ ", iterations: " + count);
//		System.out.println(" Hot end parameters: a = " + hotEnd[a] + ", b = "
//				+ hotEnd[b] + ", lag = " + hotEnd[lag]);
		aDefault = hotEnd[a];
		bDefault = hotEnd[b];
		lagDefault = hotEnd[lag];
		oderss = Math.sqrt(r) * dt / endTime;
		return hotEnd;
	}

	/**
	 * Dumb descend-the-slope optimiser to fit the pid to the ideal.
	 * 
	 * @return
	 */
	public double[] pidFit(double[] hotEnd) {
		double[] pid = { KpDefault, KiDefault, KdDefault, clampDefault };
		double r = pidRSS(pid, hotEnd);
		double[] newPid;
		double newR;
		int count = 0;
		double g;
		while (count < 200 && Math.sqrt(r) * dt / endTime > 0.2) {
			for (int i = 0; i < pid.length; i++) {
				g = pidDrByDp(r, pid, hotEnd, i);
				newPid = new double[pid.length];
				for (int j = 0; j < pid.length; j++)
					newPid[j] = pid[j];
				if (g > 0)
					newPid[i] = pid[i] - Math.abs(pid[i] * gradEst);
				else
					newPid[i] = pid[i] + Math.abs(pid[i] * gradEst);
				newR = pidRSS(newPid, hotEnd);
				if (newR < r) {
					r = newR;
					pid = newPid;
				}
			}
			count++;
		}
//		System.out.println(" PID RMS error: " + Math.sqrt(r) * dt / endTime
//				+ ", iterations: " + count);
//		System.out.println(" PID parameters: Kp = " + pid[Kp] + ", Ki = "
//				+ pid[Ki] + ", Kd = " + pid[Kd] + ", clamp = " + pid[clamp]);
		KpDefault = pid[Kp];
		KiDefault = pid[Ki];
		KdDefault = pid[Kd];
		clampDefault = pid[clamp];
		pidrss = Math.sqrt(r) * dt / endTime;
		return pid;
	}

	/**
	 * Plot curve for the experimental results
	 * 
	 * @return
	 */
	public XYSeries experimentalObservations() {
		XYSeries result = new XYSeries("Experiment");
		double t = 0;
		int i = 0;
		while (t < endTime) {
			result.add(t, experiment[i]);
			t += dt;
			i++;
		}
		return result;
	}

	/**
	 * Plot curve for the ideal profile
	 * 
	 * @return
	 */
	public XYSeries idealProfile() {
		XYSeries result = new XYSeries("Ideal");
		double t = 0;
		int i = 0;
		while (t < endTime) {
			result.add(t, ideal[i]);
			t += dt;
			i++;
		}
		return result;
	}

	/**
	 * Main program just sets everything up, computes the Hot End parameters
	 * from the experiment then uses those to optimise the PID.  Finally it plots graphs.
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		HotEndModel hem = new HotEndModel();
		for(;;)
		{
			double[] hotEnd = hem.heFit();
			hem.setIdeal();
			double[] pid = hem.pidFit(hotEnd);
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
			hem.initiate();
		}
	}
}
