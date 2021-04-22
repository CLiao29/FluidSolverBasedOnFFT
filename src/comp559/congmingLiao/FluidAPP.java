package comp559.congmingLiao;

import java.util.Random;
import java.util.concurrent.TimeUnit;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;

/**
 * 
 * @author Congming Liao
 * Base class for building the simple GUI and initialize force.
 *
 */
public class FluidAPP {
	public static final int n = 20; 
	public static Complex[] v0 = new Complex[n*(n+2)];
	public static Complex[] u0 = new Complex[n*(n+2)];
	public static Complex[] u  = new Complex[n*n];
	public static Complex[] v  = new Complex[n*n];
	public static final double visc = 0.001;
	public static final double dt   = 1;
	public static double t    = 0;
	public static final JFrame jf = new JFrame();

	//Labels to show on the frame;
	public static final JLabel label1 = new JLabel("t = "+t);
	public static final JLabel label2 = new JLabel("dt = "+dt);
	public static final JLabel label3 = new JLabel("visc = "+visc);
	
	public static void main(String[] args) {
		init();
		initFrame();
		drawGrid(jf.getGraphics());
		
	}
	
	public static void initFrame()
	{
		jf.setTitle("A Simple Fluid Solver based on the FFT");
		jf.setSize(800,800);
		jf.setResizable(false);
		jf.setLocationRelativeTo(null);
		jf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		jf.setVisible(true);
		Mouse mouse = new Mouse(jf.getGraphics());
		jf.addMouseListener(mouse);
	}
	
	public static void init() {
		Random r = new Random();
		for(int i = 0; i<n*n; i++)
		{
			v[i] = new Complex(0,0);
			u[i] = new Complex(0,0);
		}
		
		for(int j = 0; j<n*(n+2); j++)
		{	if(r.nextDouble()<0.1)
			{
			double aSign = r.nextDouble();
			double a = r.nextDouble()*30+100 ;
			a = aSign>0.5? a:-a;
			double bSign = r.nextDouble();
			double b = r.nextDouble()*30+100 ;
			b = bSign>0.5? b:-b;
			
			u0[j] = new Complex(a,0);
			v0[j] = new Complex(b,0);
			}
		else {
			u0[j] = new Complex(0,0);
			v0[j] = new Complex(0,0);
			
		}
		
		}


		
	}
	
	public static void play() {
			Fluid.stable_solve(u, v, u0, v0, visc, dt);
	}
	
	public static void drawGrid(Graphics g)
	{			
		JPanel panel = (JPanel) jf.getContentPane();
		panel.setLayout(null);
    
		panel.add(label1);
		panel.add(label2);
		panel.add(label3);
		Dimension size1  = label1.getPreferredSize();
		Dimension size2  = label2.getPreferredSize();
		Dimension size3  = label3.getPreferredSize();
		label1.setBounds(10, 10, size1.width, size1.height);
		label2.setBounds(70, 10, size2.width, size2.height);
		label3.setBounds(130, 10, size3.width, size3.height);
	
		Graphics2D g2 = (Graphics2D) g;
		int offset = 100;

		for(int i =0; i<n+1;i++)
		{
			g2.draw(new Line2D.Float(offset,offset+i*30,offset+600,offset+i*30));
			g2.draw(new Line2D.Float(offset+i*30,offset,offset+i*30,offset+600));
		}
	}
	
	public static class Mouse implements MouseListener{
		private Graphics aGraphics;
		
		public Mouse(Graphics pGraphics)
		{	
			aGraphics = pGraphics;
			
		}

		@Override
		public void mouseClicked(MouseEvent e) {
			int time =0;
			
			aGraphics.clearRect(70, 70, 800, 800);
			t += dt;
			label1.setText("t = "+t);
			drawGrid(jf.getGraphics());
			play();
			int offset =115;
			for(int i =0; i<n; i++)
			{
				for(int j =0; j<n; j++)
				{	
					Graphics2D g2 = (Graphics2D) aGraphics;
					int x1 = offset+i*30;
					int y1 = offset+j*30;
					double x = (u[i+n*j].re()*1e20);
					while(Math.abs(x)>40) {x = x/10;}
					int x2 = (int)x;
					
					double y =  (v[i+n*j].re()*1e20);
					while(Math.abs(y)>40) {y = y/10;}
					int y2  = (int)y;
					Line2D.Double line = new Line2D.Double(x1,y1,x1+x2,y1+y2);
					g2.draw(line);
					
			    
				}
			}
			
		}

		@Override
		public void mousePressed(MouseEvent e) {
	
			
		}

		@Override
		public void mouseReleased(MouseEvent e) {
			
		}

		@Override
		public void mouseEntered(MouseEvent e) {
			
		}

		@Override
		public void mouseExited(MouseEvent e) {

			
		}
	}
}
