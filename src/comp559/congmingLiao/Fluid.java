package comp559.congmingLiao;


import org.jtransforms.fft.DoubleFFT_1D;

public class Fluid {
	public static final int n = 20; 
	
	/**
	 * 
	 * @param x
	 * @return index to map in the grids
	 */
	public static int convert(double x) {
		if(x >=0.0)
		{
			return (int)x;
		}
		else {
			return -(int)(1-(x));
		}
	}
	
	public static void FFT(Complex[] a) {
		DoubleFFT_1D fftDo = new DoubleFFT_1D(a.length);
		//convert to double[] because now we are actually dealing with the real value.
		double[] input = new double[a.length];
		for(int i = 0; i < a.length; i++)
		{
			input[i] = a[i].re();
		}
		double[] output = new double[a.length*2];
		System.arraycopy(input, 0, output, 0, input.length);
		
		//FFT
		fftDo.realForwardFull(output);
		
		//convert back to Complex and put the value into a.
		for(int i = 0; i<a.length; i++)
		{
			a[i] = new Complex(output[2*i],output[2*i+1]);
		}
	}
	
	public static void IFFT(Complex[] a) {
		DoubleFFT_1D ifftDo = new DoubleFFT_1D(a.length);
		double[] ifft = new double[a.length*2];
		for(int i =0; i<a.length;i++)
		{
			ifft[2*i] = a[i].re();
			ifft[2*i+1] = a[i].im();
		}
		//IFFT
		ifftDo.complexInverse(ifft,true);
        for (int i =0; i<a.length; i++) {
        	a[i] = new Complex(ifft[i*2],0);
        }
	}
	
	
	
	public static void stable_solve( Complex[] u, Complex[] v, Complex[] u0, Complex[] v0, double visc, double dt)
	{	
		//add force
		for(int i =0; i< n*n; i++)
		{	
			u[i]= u[i].plus(u0[i].scale(dt));
			u0[i] = new Complex(u[i].re(),0);
			v[i]= v[i].plus(v0[i].scale(dt));
			v0[i] = new Complex(v[i].re(),0);
		}
		
		//advect
		double x = 0.5/n;//half grid
		double y = 0.5/n;
		for(int i =0 ; i<n; i++)
		{	
			for(int j =0; j<n; j++)
			{	
				double x0 = n*(x-u0[i+n*j].scale(dt).re())-0.5; //i-dt*u0[i+n*j].re()*n;
				double y0 = n*(y-v0[i+n*j].scale(dt).re())-0.5; //j-dt*v0[i+n*j].re()*n;
				int i0 = convert(x0);
				double s = x0-i0; 
				i0 = (n+(i0%n))%n;
				int i1 = (i0+1)%n;
				int j0 = convert(y0);
				double t = y0-j0; 
				j0 = (n+(j0%n))%n;
				int j1 = (j0+1)%n;
				u[i+n*j] = new Complex((1-s)*((1-t)*u0[i0+n*j0].re()+t*u0[i0+n*j1].re())+
                        				 s  *((1-t)*u0[i1+n*j0].re()+t*u0[i1+n*j1].re()),0);
	
				v[i+n*j] = new Complex((1-s)*((1-t)*v0[i0+n*j0].re()+t*v0[i0+n*j1].re())+
                                         s  *((1-t)*v0[i1+n*j0].re()+t*v0[i1+n*j1].re()),0);
				
				y+= 1/n;
			}
			x+=1.0/n;
		}
		
		 for (int i=0 ; i<n ; i++ )
		        for (int j=0 ; j<n ; j++ )
		            { u0[i+(n+2)*j] = new Complex(u[i+n*j].re(),0); 
		              v0[i+(n+2)*j] = new Complex(v[i+n*j].re(),0); 
		            }
		
		//transform our velocity field in the Fourier domain
		// now all the values I'm dealing with are "complex"
		FFT(u0);
		FFT(v0);
		
		
		//filter out the higher spatial frequencies by multiplying the Fourier transform by a filter whose decay depends on the magnitude of wavenumber, the viscosity and the time step
		Complex[] Ustore = new Complex[2];
		Complex[] Vstore = new Complex[2];
		for (int i=0 ; i<=n ; i+=2 ) {
		        	x = 0.5*i;
		        for (int j=0 ; j <n ; j++ ) {
		            	y = j<=n/2 ? j : j-n;
		            	double r = x*x+y*y;
		            	if ( r==0.0 ) continue;
		            	double f = Math.exp(-r*dt*visc);
		            	//projecting each velocity vector onto the line (or plane in three dimensions) perpendicular
		            	//to the wavenumber
		            	Ustore[0] = new Complex(u0[i  +(n+2)*j].re(),u0[i  +(n+2)*j].im());
		            	Vstore[0] = new Complex(v0[i  +(n+2)*j].re(),v0[i  +(n+2)*j].im());
		            	Ustore[1] = new Complex(u0[i+1+(n+2)*j].re(),u0[i+1+(n+2)*j].im());
		            	Vstore[1] = new Complex(v0[i+1+(n+2)*j].re(),v0[i+1+(n+2)*j].im());
		            	u0[i  +(n+2)*j] =  Ustore[0].scale(1-x*x/r).plus(Vstore[0].scale(-x*y/r)).scale(f);
		            	u0[i+1+(n+2)*j] =  Ustore[1].scale(1-x*x/r).plus(Vstore[1].scale(-x*y/r)).scale(f);
		            	v0[i+  (n+2)*j] =  Ustore[0].scale(-y*x/r).plus(Vstore[0].scale(1-y*y/r)).scale(f);
		            	v0[i+1+(n+2)*j] =  Ustore[1].scale(-y*x/r).plus(Vstore[1].scale(1-y*y/r)).scale(f);
		        }
		    }
	    //transform the velocity field back into the spatial domain
	    IFFT(u0);
	    IFFT(v0);
		
	    //normalize it
	    double f = 1.0/(n*n);
	    for (int i=0 ; i<n ; i++ )
	        for (int j=0 ; j<n ; j++ )
	            { u[i+n*j] = u0[i+(n+2)*j].scale(f);
	              v[i+n*j] = v0[i+(n+2)*j].scale(f); }
	    

	}
}

