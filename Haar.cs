/*
 * Created by SharpDevelop.
 * User: jeffA
 * Date: 12/7/2015
 * Time: 11:41 AM
 * 
 * Filename: Haar.cs
 * This is part of a wavelet library developed for SDR# for use by SDR# plugins.
 * This is simply a C# port of a wavelet library developed by Ian Kaplan.  His original
 * source code has been included at the branch of the ifdef directive.
 * This code does wavelet transforms using the Haar basis function.
 */

#define CSHARPCODE
 
using System;
using SDRSharp.WaveletLib;

#if CSHARPCODE 

namespace SDRSharp.WaveletLib
	
{
    public class Haar
    {

  		private enum bogus{ numPts = 4 };
  		private LinTerp fourPtInterp;

  /**
    Copy four points or <i>N</i> (which ever is less) data points from
    <i>vec</i> into <i>d</i> These points are the "known" points used
    in the polynomial interpolation.

    @param vec[] the input data set on which the wavelet is calculated
    @param d[] an array into which <i>N</i> data points, starting at
               <i>start</i> are copied.
    @param N the number of polynomial interpolation points
    @param start the index in <i>vec</i> from which copying starts
   */
  		private void fill( double[] vec, double[] d, int N, int start ) //edit: T& repby double[]
  		{
  			int n = 4;	//was assigned to enumeration bogus;  changed to const 4
    		if (n > N)
      			n = N;
    		int end = start + n;
    		int j = 0;

    		for (int i = start; i < end; i++) {
      			d[j] = vec[i];
      			j++;
    		}
  		} // fill

   /**
      The update step of the wavelet Lifting Scheme

      <b>transform</b>

      In the Lifting Scheme transform the update step follows
      the predict step.  After the predict step, the differences
      have been calculated in-place writing over the even (b)
      values.  The update step calculates the Haar average using
      an even/odd pair.  The average is placed in the even member
      of the pair.

    */
   		protected void update( double[] vec, int N, LiftBase.transDirection direction)  //edit: T& repby double[]
   		{
      		int i;
      		int half = N >> 1;
      		int j = half;
      		for (i = 0; i < half; i++) {
         		if (direction == LiftBase.transDirection.forward) {  // forward transform stage
            		vec[i] = vec[i] + (vec[j]/2.0);
         		}
         		else if (direction == LiftBase.transDirection.inverse) { // inverse transform step
            		vec[i] = vec[i] - (vec[j]/2.0);
         		}
	 			else {
	   				Console.Write("update: bad direction value\n");
	   				break;
	 			}
         		j++;
      		} // for
   		} // update


  /**
    <p>
    Predict an odd point from the even points, using 4-point
    polynomial interpolation.
    </p>
    <p>
    The four points used in the polynomial interpolation are
    the even points.  We pretend that these four points
    are located at the x-coordinates 0,1,2,3.  The first 
    odd point interpolated will be located between the first
    and second even point, at 0.5.  The next N-3 points are
    located at 1.5 (in the middle of the four points).
    The last two points are located at 2.5 and 3.5.
    For complete documentation see
    </p>
    <pre>
  <a href="http://www.bearcave.com/misl/misl_tech/wavelets/lifting/index.html">
  http://www.bearcave.com/misl/misl_tech/wavelets/lifting/index.html</a>
    </pre>

    <p>
    The difference between the predicted (interpolated) value
    and the actual odd value replaces the odd value in the
    forward transform.
    </p>

    <p>
    As the recursive steps proceed, N will eventually be 4 and
    then 2.  When N = 4, linear interpolation is used.
    When N = 2, Haar interpolation is used (the prediction for
    the odd value is that it is equal to the even value).
    </p>

    @param vec the input data on which the forward or inverse
               transform is calculated.
    @param N the area of vec over which the transform is calculated
    @param direction forward or inverse transform

   */
  		protected void interp( double[] vec, int N, LiftBase.transDirection direction )  //edit: T& repby double[]
  		{
    		int half = N >> 1;
    		double[] d = new Double[4];

    		for (int i = 0; i < half; i++) {
      			double predictVal;

      			if (i == 0) {
					if (half == 1) {
	  					// e.g., N == 2, and we use Haar interpolation
	  					predictVal = vec[0];
					}
					else {
	  					fill( vec, d, N, 0 );
	  					predictVal = fourPtInterp.interpPoint( 0.5, half, d );
					}
      			}
      			else if (i == 1) {
					predictVal = fourPtInterp.interpPoint( 1.5, half, d );
      			}
      			else if (i == half-2) {
					predictVal = fourPtInterp.interpPoint( 2.5, half, d );
      			}
      			else if (i == half-1) {
					predictVal = fourPtInterp.interpPoint( 3.5, half, d );
      			}
      			else {
					fill( vec, d, N, i-1);
					predictVal = fourPtInterp.interpPoint( 1.5, half, d );
      			}

      			int j = i + half;
      			if (direction == LiftBase.transDirection.forward) {
					vec[j] = vec[j] - predictVal;
      			}
      			else if (direction == LiftBase.transDirection.inverse) {
					vec[j] = vec[j] + predictVal;
      			}
      			else {
					Console.Write("interp: bad direction value\n");
					break;
      			}
    		} // for
  		} // interp

   /**
      Predict step of the wavelet Lifting Scheme.

      The term "predict" comes from <i>Building Your Own Wavelets
      at Home</i>.

      <b>transform</b>

      In the wavelet transform, calculate the difference between even
      and odd element pairs.  This is a slightly modified version of
      the classic haar difference.  If the even elements are labeled
      as <i>a</i> and the odd elements are labeled as <i>b</i> (where
      we start counting at zero) the difference is simply

<pre>
       d<sub>i</sub> = b<sub>i</sub> - a<sub>i</sub>
</pre>

      Since an "in-place" algorithm is used, where we update the
      odd elements with the difference, we get
<pre>
       b<sub>i</sub> = b<sub>i</sub> - a<sub>i</sub>
</pre>

      <b>inverse transform</b>

      Reverse the transform predict step by adding the average
      (even element) to the difference (odd element).

    */
   		protected void predict( double[] vec, int N, LiftBase.transDirection direction)  //edit: T& repby double[]
   		{
      		int i;
      		int half = N >> 1;
      		int j = 0;
      		for (i = half; i < N; i++) {
         		if (direction == LiftBase.transDirection.forward) { // forward transform stage
	   				vec[i] = vec[i] - vec[j];
         		}
         		else if (direction == LiftBase.transDirection.inverse ) { // inverse transform stage
	   				vec[i] = vec[i] + vec[j];
         		}
	 			else {
	   				Console.Write("predict: bad direction\n");
	   				break;
	 			}
         		j++;
      		}
   		} // predict

   		public void forwardStep( double[] vec, int n )  //edit: T& repby double[], const int n repby int n
  		{
   			LiftBase tmp = new LiftBase();
    		tmp.split( vec, n );
    		predict( vec, n, LiftBase.transDirection.forward );
    		update( vec, n, LiftBase.transDirection.forward );
    		interp( vec, n, LiftBase.transDirection.forward );
  		} // forwardStep

   		public virtual void inverseStep( double[] vec, int n )  //edit: T& repby double[], const int n repby int n
  		{
    		interp( vec, n, LiftBase.transDirection.inverse );
    		update( vec, n, LiftBase.transDirection.inverse );
    		predict( vec, n, LiftBase.transDirection.inverse );
    		LiftBase tmp = new LiftBase();
    		tmp.merge( vec, n );
    		
  		}
    }

#else
/*
<h4>
   Copyright and Use
</h4>

<p>
   You may use this source code without limitation and without
   fee as long as you include:
</p>
<blockquote>
     This software was written and is copyrighted by Ian Kaplan, Bear
     Products International, www.bearcave.com, 2002.
</blockquote>
<p>
   This software is provided "as is", without any warranty or
   claim as to its usefulness.  Anyone who uses this source code
   uses it at their own risk.  Nor is any support provided by
   Ian Kaplan and Bear Products International.
<p>
   Please send any bug fixes or suggested source changes to:
<pre>
     iank@bearcave.com
</pre>

  @author Ian Kaplan


 */

//#include <vector>

//#include "liftbase.h"
//#include "polyinterp.h"


/**
   This class implements the wavelet Lifting Scheme with polynomial
   interpolation.

   This version of wavelets is based on Wim Sweldens' tutorial paper
   <i>Building Your Own Wavelets at Home</i>.  This tutorial was
   presented at SIGGraph.  The tutorial is available on the 
   Web.

   One of the attractive features of Lifting Scheme wavelets is that
   the transform and the inverse transform are exact mirrors of each other.
   To arrive at the inverse transform the order of the operations is
   reversed and subtraction and addition operations are exchanged.
   This allows a generic Lifting Scheme base class to be constructed,
   where the ordering is defined by the sub-classes.

 */
template<class T>
class polyHaar : public liftbase<T, double> {

public:
  /**
    polyHaar class constructor
   */
  polyHaar() {}
  ~polyHaar() {}
  // declare, but don't define to disallow copy constructor
  polyHaar(const polyHaar &rhs );

private:
  typedef enum { numPts = 4 } bogus;
  polyinterp fourPtInterp;

  /**
    Copy four points or <i>N</i> (which ever is less) data points from
    <i>vec</i> into <i>d</i> These points are the "known" points used
    in the polynomial interpolation.

    @param vec[] the input data set on which the wavelet is calculated
    @param d[] an array into which <i>N</i> data points, starting at
               <i>start</i> are copied.
    @param N the number of polynomial interpolation points
    @param start the index in <i>vec</i> from which copying starts
   */
  void fill( T &vec, double d[], int N, int start )
  {
    int n = numPts;
    if (n > N)
      n = N;
    int end = start + n;
    int j = 0;

    for (int i = start; i < end; i++) {
      d[j] = vec[i];
      j++;
    }
  } // fill

protected:


   /**
      The update step of the wavelet Lifting Scheme

      <b>transform</b>

      In the Lifting Scheme transform the update step follows
      the predict step.  After the predict step, the differences
      have been calculated in-place writing over the even (b)
      values.  The update step calculates the Haar average using
      an even/odd pair.  The average is placed in the even member
      of the pair.

    */
   void update( T& vec, int N, transDirection direction)
   {
      int i;
      int half = N >> 1;
      int j = half;
      for (i = 0; i < half; i++) {
         if (direction == forward) {  // forward transform stage
            vec[i] = vec[i] + (vec[j]/2.0);
         }
         else if (direction == inverse) { // inverse transform step
            vec[i] = vec[i] - (vec[j]/2.0);
         }
	 else {
	   Cionsole.Write("update: bad direction value\n");
	   break;
	 }
         j++;
      } // for
   } // update


  /**
    <p>
    Predict an odd point from the even points, using 4-point
    polynomial interpolation.
    </p>
    <p>
    The four points used in the polynomial interpolation are
    the even points.  We pretend that these four points
    are located at the x-coordinates 0,1,2,3.  The first 
    odd point interpolated will be located between the first
    and second even point, at 0.5.  The next N-3 points are
    located at 1.5 (in the middle of the four points).
    The last two points are located at 2.5 and 3.5.
    For complete documentation see
    </p>
    <pre>
  <a href="http://www.bearcave.com/misl/misl_tech/wavelets/lifting/index.html">
  http://www.bearcave.com/misl/misl_tech/wavelets/lifting/index.html</a>
    </pre>

    <p>
    The difference between the predicted (interpolated) value
    and the actual odd value replaces the odd value in the
    forward transform.
    </p>

    <p>
    As the recursive steps proceed, N will eventually be 4 and
    then 2.  When N = 4, linear interpolation is used.
    When N = 2, Haar interpolation is used (the prediction for
    the odd value is that it is equal to the even value).
    </p>

    @param vec the input data on which the forward or inverse
               transform is calculated.
    @param N the area of vec over which the transform is calculated
    @param direction forward or inverse transform

   */
  void interp( T &vec, int N, transDirection direction )
  {
    int half = N >> 1;
    double d[4];

    for (int i = 0; i < half; i++) {
      double predictVal;

      if (i == 0) {
	if (half == 1) {
	  // e.g., N == 2, and we use Haar interpolation
	  predictVal = vec[0];
	}
	else {
	  fill( vec, d, N, 0 );
	  predictVal = fourPtInterp.interpPoint( 0.5, half, d );
	}
      }
      else if (i == 1) {
	predictVal = fourPtInterp.interpPoint( 1.5, half, d );
      }
      else if (i == half-2) {
	predictVal = fourPtInterp.interpPoint( 2.5, half, d );
      }
      else if (i == half-1) {
	predictVal = fourPtInterp.interpPoint( 3.5, half, d );
      }
      else {
	fill( vec, d, N, i-1);
	predictVal = fourPtInterp.interpPoint( 1.5, half, d );
      }

      int j = i + half;
      if (direction == forward) {
	vec[j] = vec[j] - predictVal;
      }
      else if (direction == inverse) {
	vec[j] = vec[j] + predictVal;
      }
      else {
	Console.Write("interp: bad direction value\n");
	break;
      }
    } // for
  } // interp

   /**
      Predict step of the wavelet Lifting Scheme.

      The term "predict" comes from <i>Building Your Own Wavelets
      at Home</i>.

      <b>transform</b>

      In the wavelet transform, calculate the difference between even
      and odd element pairs.  This is a slightly modified version of
      the classic haar difference.  If the even elements are labeled
      as <i>a</i> and the odd elements are labeled as <i>b</i> (where
      we start counting at zero) the difference is simply

<pre>
       d<sub>i</sub> = b<sub>i</sub> - a<sub>i</sub>
</pre>

      Since an "in-place" algorithm is used, where we update the
      odd elements with the difference, we get
<pre>
       b<sub>i</sub> = b<sub>i</sub> - a<sub>i</sub>
</pre>

      <b>inverse transform</b>

      Reverse the transform predict step by adding the average
      (even element) to the difference (odd element).

    */
   void predict( T& vec, int N, transDirection direction)
   {
      int i;
      int half = N >> 1;
      int j = 0;
      for (i = half; i < N; i++) {
         if (direction == forward) { // forward transform stage
	   vec[i] = vec[i] - vec[j];
         }
         else if (direction == inverse ) { // inverse transform stage
	   vec[i] = vec[i] + vec[j];
         }
	 else {
	   Console.Write("predict: bad direction\n");
	   break;
	 }
         j++;
      }
   } // predict


public:

  void forwardStep( T& vec, const int n )
  {
    split( vec, n );
    predict( vec, n, forward );
    update( vec, n, forward );
    interp( vec, n, forward );
  } // forwardStep

  virtual void inverseStep( T& vec, const int n )
  {
    interp( vec, n, inverse );
    update( vec, n, inverse );
    predict( vec, n, inverse );
    merge( vec, n );
  }

}; // polyHaar

#endif

}//namespace