﻿/*
 * Created by SharpDevelop.
 * User: jeffA
 * Date: 12/7/2015
 * Time: 11:44 AM
 * 
 * Filename: Daubechies.cs
 * This is part of a wavelet library developed for SDR# for use by SDR# plugins.
 * This is simply a C# port of a wavelet library developed by Ian Kaplan.  His original
 * source code has been included at the branch of the ifdef directive.
 * This code does wavelet transforms using the Daubechies D4 basis function.
 */
 
 
  
 
#define CSHARPCODE

using System;
using System.Diagnostics;
using SDRSharp.WaveletLib;

#if CSHARPCODE 

namespace SDRSharp.WaveletLib
{
    public unsafe class Daubechies
    {
    	protected void predict( double[] vec, int N, LiftBase.transDirection direction )  //edit: T& repby double[]
  		{
  			Debug.Assert( false );
  		} // predict
  
    	protected virtual void update( double[] vec, int N, LiftBase.transDirection direction )  //edit: T& repby double[]
  		{
    		Debug.Assert( false );
  		} // update
  
  		/** forward transform scaling coefficients */
  		private double h0, h1, h2, h3;
  
  
 	 	/** forward transform wave coefficients */
  		private double g0, g1, g2, g3;
  
  		private double Ih0, Ih1, Ih2, Ih3;
  		private double Ig0, Ig1, Ig2, Ig3;
        
        /* Forward Daubechies D4 transform step */
        public void forwardStep( double[] a, int n )  //edit: T& repby double[], const int n repby int n
  		{
    		if (n >= 4) {
      			int i, j;
      			int half = n >> 1;	//const int repby int
      
      			double[] tmp = new double[n];
      
      			for (i = 0, j = 0; j < n-3; j += 2, i++) {
					tmp[i]      = a[j]*h0 + a[j+1]*h1 + a[j+2]*h2 + a[j+3]*h3;
					tmp[i+half] = a[j]*g0 + a[j+1]*g1 + a[j+2]*g2 + a[j+3]*g3;
      			}
      
      			tmp[i]      = a[n-2]*h0 + a[n-1]*h1 + a[0]*h2 + a[1]*h3;
      			tmp[i+half] = a[n-2]*g0 + a[n-1]*g1 + a[0]*g2 + a[1]*g3;
      
      			for (i = 0; i < n; i++) {
					a[i] = tmp[i];
      			}
      			//delete [] tmp;
    		}
  		}
        
		/**
			Forward Daubechies D4 transform step, where the locations
    		for the high and low pass filters are reversed.
    	*/
    	void forwardStepRev( double[] a, int n )  //edit: T& repby double[], const int n repby int n
  		{
    		if (n >= 4) {
      			int i, j;
      			int half = n >> 1;
      
      			double[] tmp = new double[n];
      
      			for (i = 0, j = 0; j < n-3; j += 2, i++) {
					tmp[i+half] = a[j]*h0 + a[j+1]*h1 + a[j+2]*h2 + a[j+3]*h3;
					tmp[i]      = a[j]*g0 + a[j+1]*g1 + a[j+2]*g2 + a[j+3]*g3;
      			}
      
      			tmp[i+half] = a[n-2]*h0 + a[n-1]*h1 + a[0]*h2 + a[1]*h3;
      			tmp[i]      = a[n-2]*g0 + a[n-1]*g1 + a[0]*g2 + a[1]*g3;
      
      			for (i = 0; i < n; i++) {
					a[i] = tmp[i];
      			}
      			//delete [] tmp;
  			}
    	}
  
  		/* Inverse Daubechies D4 transform */
  		public void inverseStep( double[] a, int n )  //edit: T& repby double[], const int n repby int n
  		{
   			if (n >= 4) {
    			int i, j;
      			int half = n >> 1;
      			int halfPls1 = half + 1;
      
      			double[] tmp = new double[n];
      
      			//      last smooth val  last coef.  first smooth  first coef
      			tmp[0] = a[half-1]*Ih0 + a[n-1]*Ih1 + a[0]*Ih2 + a[half]*Ih3;
      			tmp[1] = a[half-1]*Ig0 + a[n-1]*Ig1 + a[0]*Ig2 + a[half]*Ig3;
      			for (i = 0, j = 2; i < half-1; i++) {
				//     smooth val     coef. val       smooth val    coef. val
					tmp[j++] = a[i]*Ih0 + a[i+half]*Ih1 + a[i+1]*Ih2 + a[i+halfPls1]*Ih3;
					tmp[j++] = a[i]*Ig0 + a[i+half]*Ig1 + a[i+1]*Ig2 + a[i+halfPls1]*Ig3;
      			}
      			for (i = 0; i < n; i++) {
					a[i] = tmp[i];
      			}
      			//delete [] tmp;
  			}
  		} // inverseStep
  

  		/**
    		Initialize the filter constants used in the Daubechies 
    		D4 transform.
   		*/
  		public Daubechies() 
  		{
  			double sqrt_3 = Math.Sqrt( 3 );
    		double denom = 4 * Math.Sqrt( 2 );
    
    		//
    		// forward transform scaling (smoothing) coefficients
    		//
    		h0 = (1 + sqrt_3)/denom;
    		h1 = (3 + sqrt_3)/denom;
    		h2 = (3 - sqrt_3)/denom;
    		h3 = (1 - sqrt_3)/denom;
    		//
    		// forward transform wavelet coefficients (a.k.a. high
    		// pass filter coefficients)
    		//
    		g0 =  h3;
    		g1 = -h2;
    		g2 =  h1;
    		g3 = -h0;
    
    		Ih0 = h2;
    		Ih1 = g2;  // h1
    		Ih2 = h0;
    		Ih3 = g0;  // h3
    
    		Ig0 = h3;
    		Ig1 = g3;  // -h0
    		Ig2 = h1;
    		Ig3 = g1;  // -h2
  		}

  		/**
    		Forward Daubechies D4 transform
   		*/
   		public void forwardTrans( double[] ts, int N )  //edit: T& repby double[]
  		{
    		int n;
    		for (n = N; n >= 4; n >>= 1) {
    			forwardStep( ts, n );
    		}
  		} // forwardTrans
  
  
  	/**
    	Inverse Daubechies D4 transform
   	*/
   	public void inverseTrans( double[] coef, int N )  //edit: T& repby double[]
   	{
    	int n;
    	for (n = 4; n <= N; n <<= 1) {
      		inverseStep( coef, n );
    	}
  	} // inverseTrans
    
}
#else
  //#include <math.h>

/**
  Daubechies D4 wavelet transform (D4 denotes four coefficients)

  I have to confess up front that the comment here does not even come
  close to describing wavelet algorithms and the Daubechies D4
  algorithm in particular.  I don't think that it can be described in
  anything less than a journal article or perhaps a book.  I even have
  to apologize for the notation I use to describe the algorithm, which
  is barely adequate.  But explaining the correct notation would take
  a fair amount of space as well.  This comment really represents some
  notes that I wrote up as I implemented the code.  If you are
  unfamiliar with wavelets I suggest that you look at the bearcave.com
  web pages and at the wavelet literature.  I have yet to see a really
  good reference on wavelets for the software developer.  The best
  book I can recommend is <i>Ripples in Mathematics</i> by Jensen and
  Cour-Harbo.

  All wavelet algorithms have two components, a wavelet function and a
  scaling function.  These are sometime also referred to as high pass
  and low pass filters respectively.

  The wavelet function is passed two or more samples
  and calculates a wavelet coefficient.  In the case of
  the Haar wavelet this is 

  <pre>
  coef<sub>i</sub> = odd<sub>i</sub> - even<sub>i</sub>
  or 
  coef<sub>i</sub> = 0.5 * (odd<sub>i</sub> - even<sub>i</sub>)
  </pre>

  depending on the version of the Haar algorithm used.

  The scaling function produces a smoother version of the
  original data.  In the case of the Haar wavelet algorithm
  this is an average of two adjacent elements.

  The Daubechies D4 wavelet algorithm also has a wavelet
  and a scaling function.  The coefficients for the 
  scaling function are denoted as h<sub>i</sub> and the
  wavelet coefficients are g<sub>i</sub>.


  Mathematicians like to talk about wavelets in terms of
  a wavelet algorithm applied to an infinite data set.
  In this case one step of the forward transform can be expressed
  as the infinite matrix of wavelet coefficients
  represented below multiplied by the infinite signal
  vector.

  <pre>
     a<sub>i</sub> = ...h0,h1,h2,h3, 0, 0, 0, 0, 0, 0, 0, ...   s<sub>i</sub>
     c<sub>i</sub> = ...g0,g1,g2,g3, 0, 0, 0, 0, 0, 0, 0, ...   s<sub>i+1</sub>
   a<sub>i+1</sub> = ...0, 0, h0,h1,h2,h3, 0, 0, 0, 0, 0, ...   s<sub>i+2</sub>
   c<sub>i+1</sub> = ...0, 0, g0,g1,g2,g3, 0, 0, 0, 0, 0, ...   s<sub>i+3</sub>
   a<sub>i+2</sub> = ...0, 0, 0, 0, h0,h1,h2,h3, 0, 0, 0, ...   s<sub>i+4</sub>
   c<sub>i+2</sub> = ...0, 0, 0, 0, g0,g1,g2,g3, 0, 0, 0, ...   s<sub>i+5</sub>
   a<sub>i+3</sub> = ...0, 0, 0, 0, 0, 0, h0,h1,h2,h3, 0, ...   s<sub>i+6</sub>
   c<sub>i+3</sub> = ...0, 0, 0, 0, 0, 0, g0,g1,g2,g3, 0, ...   s<sub>i+7</sub>
  </pre>

  The dot product (inner product) of the infinite vector and
  a row of the matrix produces either a smoother version of the
  signal (a<sub>i</sub>) or a wavelet coefficient (c<sub>i</sub>).

  In an ordered wavelet transform, the smoothed (a<sub>i</sub>) are
  stored in the first half of an <i>n</i> element array region.  The
  wavelet coefficients (c<sub>i</sub>) are stored in the second half
  the <i>n</i> element region.  The algorithm is recursive.  The
  smoothed values become the input to the next step.

  The transpose of the forward transform matrix above is used
  to calculate an inverse transform step.  Here the dot product is
  formed from the result of the forward transform and an inverse
  transform matrix row.

  <pre>
      s<sub>i</sub> = ...h2,g2,h0,g0, 0, 0, 0, 0, 0, 0, 0, ...  a<sub>i</sub>
    s<sub>i+1</sub> = ...h3,g3,h1,g1, 0, 0, 0, 0, 0, 0, 0, ...  c<sub>i</sub>
    s<sub>i+2</sub> = ...0, 0, h2,g2,h0,g0, 0, 0, 0, 0, 0, ...  a<sub>i+1</sub>
    s<sub>i+3</sub> = ...0, 0, h3,g3,h1,g1, 0, 0, 0, 0, 0, ...  c<sub>i+1</sub>
    s<sub>i+4</sub> = ...0, 0, 0, 0, h2,g2,h0,g0, 0, 0, 0, ...  a<sub>i+2</sub>
    s<sub>i+5</sub> = ...0, 0, 0, 0, h3,g3,h1,g1, 0, 0, 0, ...  c<sub>i+2</sub>
    s<sub>i+6</sub> = ...0, 0, 0, 0, 0, 0, h2,g2,h0,g0, 0, ...  a<sub>i+3</sub>
    s<sub>i+7</sub> = ...0, 0, 0, 0, 0, 0, h3,g3,h1,g1, 0, ...  c<sub>i+3</sub>
  </pre>

  Using a standard dot product is grossly inefficient since most
  of the operands are zero.  In practice the wavelet coefficient 
  values are moved along the signal vector and a four element 
  dot product is calculated.  Expressed in terms of arrays, for
  the forward transform this would be:

  <pre>
  a<sub>i</sub> = s[i]*h0 + s[i+1]*h1 + s[i+2]*h2 + s[i+3]*h3
  c<sub>i</sub> = s[i]*g0 + s[i+1]*g1 + s[i+2]*g2 + s[i+3]*g3
  </pre>

  This works fine if we have an infinite data set, since we don't
  have to worry about shifting the coefficients "off the end" of
  the signal.

  I sometimes joke that I left my infinite data set in my other bear
  suit.  The only problem with the algorithm described so far is that
  we don't have an infinite signal.  The signal is finite.  In fact
  not only must the signal be finite, but it must have a power of two
  number of elements.

  If i=N-1, the i+2 and i+3 elements will be beyond the end of 
  the array.  There are a number of methods for handling the 
  wavelet edge problem.  This version of the algorithm acts 
  like the data is periodic, where the data at the start of 
  the signal wraps around to the end.

  This algorithm uses a temporary array.  A Lifting Scheme version of
  the Daubechies D4 algorithm does not require a temporary.  The
  matrix discussion above is based on material from <i>Ripples in
  Mathematics</i>, by Jensen and Cour-Harbo.  Any error are mine.

  <b>Author</b>: Ian Kaplan<br>
  <b>Use</b>: You may use this software for any purpose as long
  as I cannot be held liable for the result.  Please credit me
  with authorship if use use this source code.

  This comment is formatted for the doxygen documentation generator

 */
template <class T>
class Daubechies : public liftbase<T, double> {
  
protected:
  void predict( T& vec, int N, transDirection direction )
  {
    assert( false );
  } // predict
  
  virtual void update( T& vec, int N, transDirection direction )
  {
    assert( false );
  } // update
  
private:
  /** forward transform scaling coefficients */
  double h0, h1, h2, h3;
  /** forward transform wave coefficients */
  double g0, g1, g2, g3;
  
  double Ih0, Ih1, Ih2, Ih3;
  double Ig0, Ig1, Ig2, Ig3;

public:
  
  /**
    Forward Daubechies D4 transform step
    
    */
  void forwardStep( T& a, const int n )
  {
    if (n >= 4) {
      int i, j;
      const int half = n >> 1;
      
      double* tmp = new double[n];
      
      for (i = 0, j = 0; j < n-3; j += 2, i++) {
	tmp[i]      = a[j]*h0 + a[j+1]*h1 + a[j+2]*h2 + a[j+3]*h3;
	tmp[i+half] = a[j]*g0 + a[j+1]*g1 + a[j+2]*g2 + a[j+3]*g3;
      }
      
      tmp[i]      = a[n-2]*h0 + a[n-1]*h1 + a[0]*h2 + a[1]*h3;
      tmp[i+half] = a[n-2]*g0 + a[n-1]*g1 + a[0]*g2 + a[1]*g3;
      
      for (i = 0; i < n; i++) {
	a[i] = tmp[i];
      }
      delete [] tmp;
    }
  }
  
  /**
    Forward Daubechies D4 transform step, where the locations
    for the high and low pass filters are reversed.
    
    */
  void forwardStepRev( T& a, const int n )
  {
    if (n >= 4) {
      int i, j;
      const int half = n >> 1;
      
      double* tmp = new double[n];
      
      for (i = 0, j = 0; j < n-3; j += 2, i++) {
	tmp[i+half] = a[j]*h0 + a[j+1]*h1 + a[j+2]*h2 + a[j+3]*h3;
	tmp[i]      = a[j]*g0 + a[j+1]*g1 + a[j+2]*g2 + a[j+3]*g3;
      }
      
      tmp[i+half] = a[n-2]*h0 + a[n-1]*h1 + a[0]*h2 + a[1]*h3;
      tmp[i]      = a[n-2]*g0 + a[n-1]*g1 + a[0]*g2 + a[1]*g3;
      
      for (i = 0; i < n; i++) {
	a[i] = tmp[i];
      }
      delete [] tmp;
    }
  }
  
  /**
    Inverse Daubechies D4 transform
    */
  void inverseStep( T& a, const int n )
  {
    if (n >= 4) {
      int i, j;
      const int half = n >> 1;
      const int halfPls1 = half + 1;
      
      double* tmp = new double[n];
      
      //      last smooth val  last coef.  first smooth  first coef
      tmp[0] = a[half-1]*Ih0 + a[n-1]*Ih1 + a[0]*Ih2 + a[half]*Ih3;
      tmp[1] = a[half-1]*Ig0 + a[n-1]*Ig1 + a[0]*Ig2 + a[half]*Ig3;
      for (i = 0, j = 2; i < half-1; i++) {
	//     smooth val     coef. val       smooth val    coef. val
	tmp[j++] = a[i]*Ih0 + a[i+half]*Ih1 + a[i+1]*Ih2 + a[i+halfPls1]*Ih3;
	tmp[j++] = a[i]*Ig0 + a[i+half]*Ig1 + a[i+1]*Ig2 + a[i+halfPls1]*Ig3;
      }
      for (i = 0; i < n; i++) {
	a[i] = tmp[i];
      }
      delete [] tmp;
    }
  } // inverseStep
  

  /**
    Initialize the filter constants used in the Daubechies 
    D4 transform.
   */
  Daubechies() 
  {
    const double sqrt_3 = sqrt( 3 );
    const double denom = 4 * sqrt( 2 );
    
    //
    // forward transform scaling (smoothing) coefficients
    //
    h0 = (1 + sqrt_3)/denom;
    h1 = (3 + sqrt_3)/denom;
    h2 = (3 - sqrt_3)/denom;
    h3 = (1 - sqrt_3)/denom;
    //
    // forward transform wavelet coefficients (a.k.a. high
    // pass filter coefficients)
    //
    g0 =  h3;
    g1 = -h2;
    g2 =  h1;
    g3 = -h0;
    
    Ih0 = h2;
    Ih1 = g2;  // h1
    Ih2 = h0;
    Ih3 = g0;  // h3
    
    Ig0 = h3;
    Ig1 = g3;  // -h0
    Ig2 = h1;
    Ig3 = g1;  // -h2
  }

  /**
    Forward Daubechies D4 transform
   */
  void forwardTrans( T& ts, int N )
  {
    int n;
    for (n = N; n >= 4; n >>= 1) {
      forwardStep( ts, n );
    }
  } // forwardTrans
  
  
  /**
    Inverse Daubechies D4 transform
   */
  void inverseTrans( T& coef, int N )
  {
    int n;
    for (n = 4; n <= N; n <<= 1) {
      inverseStep( coef, n );
    }
  } // inverseTrans
  
}; // Daubechies

#endif
}  //namespace