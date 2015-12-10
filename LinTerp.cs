/*
 * Created by SharpDevelop.
 * User: jeffA
 * Date: 12/7/2015
 * Time: 11:43 AM
 * 
 * Filename: LinTerp.cs
 * This is part of a wavelet library developed for SDR# for use by SDR# plugins.
 * This is simply a C# port of a wavelet library developed by Ian Kaplan.  His original
 * source code has been included at the branch of the ifdef directive.
 * This code does wavelet transforms using the Linear Interpolation basis function.
 */
 
#define CSHARPCODE

using System;
using SDRSharp.WaveletLib;

#if CSHARPCODE 

namespace SDRSharp.WaveletLib
{
    public class LinTerp
    {
    	private double[,] fourPointTable = new double[4,4];
    	private double[,] twoPointTable = new double[4,4];

   /**
     The polynomial interpolation algorithm assumes that the known
     points are located at x-coordinates 0, 1,.. N-1.  An interpolated
     point is calculated at <b><i>x</i></b>, using N coefficients.  The
     polynomial coefficients for the point <b><i>x</i></b> can be 
     calculated staticly, using the Lagrange method.

     @param x the x-coordinate of the interpolated point
     @param N the number of polynomial points.
     @param c[] an array for returning the coefficients

     */
    	private void lagrange( double x, int N, double[] c )
  		{
    		double num, denom;

    		for (int i = 0; i < N; i++) {
      			num = 1;
      			denom = 1;
      			for (int k = 0; k < N; k++) {
					if (i != k) {
	  					num = num * (x - k);
	  					denom = denom * (i - k);
					}
      			} // for k
      			c[i] = num / denom;
    		} // for i
  		} // lagrange


  /**
    For a given N-point polynomial interpolation, fill the coefficient
    table, for points 0.5 ... (N-0.5).

   */
  		private void fillTable( int N, double[,] table )//edit: const int n repby int n
  		{
    		double x;
    		double n = N;
    		int i = 0;
    		double[,] tmptable = new double[4,4];
    		tmptable = table;
    		double[] tmp = new double[4];

    		for (x = 0.5; x < n; x = x + 1.0) {
    			tmp[0] = tmptable[i,0];
    			tmp[1] = tmptable[i,1];
    			tmp[2] = tmptable[i,2];
    			tmp[3] = tmptable[i,3];
      			lagrange( x, N, tmp );
      			i++;
    		}
  		} // fillTable


  /**
    Print an N x N table polynomial coefficient table
   */
  		private void printTable( double[,] table, int N)
  		{
    		Console.Write("%d-point interpolation table:\n", N);
    		double x = 0.5;
    		for (int i = 0; i < N; i++) {
      			Console.Write("%4.2f: ", x);
      			for (int j = 0; j < N; j++) {
					Console.Write("%6.4f", table[i,j] );
					if (j < N-1)
	  					Console.Write(", ");
      			}
      			Console.Write("\n");
      			x = x + 1.0;
    		}    
  		}


  /**
    Print the 4-point and 2-point polynomial coefficient
    tables.
   */
  		private void printTables()
  		{
    		printTable( fourPointTable, 4 );
    		printTable( twoPointTable, 2 );
    		Console.Write("\n");
  		} // printTables


    /**
    For the polynomial interpolation point x-coordinate 
    <b><i>x</i></b>, return the associated polynomial
    interpolation coefficients.

    @param x the x-coordinate for the interpolated pont
    @param n the number of polynomial interpolation points
    @param c[] an array to return the polynomial coefficients
   */
  		private void getCoef( double x, int n, double[] c)
  		{
    		int j = (int)x;

    		if (j < 0 || j >= n) {
      			Console.Write("poly::getCoef: n = %d, bad x value = %f\n", n, x);
    		}

    		if (n == 2) {
      			c[2] = 0.0;
      			c[3] = 0.0;
    		}
    		else if (n != 4) {
      			Console.Write("poly::getCoef: bad value for N\n");
      			return;
    		}

    		for (int i = 0; i < n; i++) {
      			if (n == 4) {
					c[i] = fourPointTable[j,i];
      			}
      			else { // n == 2
					c[i] = twoPointTable[j,i];
      			}
    		}

  		} // getCoef


   /**
      The polyinterp (polynomial interpolation) class constructor
      calculates the coefficients for a four point polynomial
      interpolated at -0.5, 0.5, 1.5, 2.5 and 3.5.
   */
   		public LinTerp() 
   		{
    		// Fill in the 4-point polynomial interplation table
    		// for the points 0.5, 1.5, 2.5, 3.5
    		fillTable( 4, fourPointTable );
    
    		// Fill in the 2-point polynomial interpolation table
    		// for 0.5 and 1.5
    		fillTable( 2, twoPointTable );
   		} // polyinterp class constructor


  /**
    Given four points at the x,y coordinates {0,d<sub>0</sub>},
    {1,d<sub>1</sub>}, {2,d<sub>2</sub>}, {3,d<sub>3</sub>}
    return the y-coordinate value for the polynomial interpolated
    point at <b><i>x</i></b>.
    
    @param x the x-coordinate for the point to be interpolated
    @param N the number of interpolation points
    @param d[] an array containing the y-coordinate values for
           the known points (which are located at x-coordinates
	   0..N-1).
   */
  public double interpPoint( double x, int N, double[] d)
  		{
  	double[] c = new Double[4];
    		double point = 0;
    
    		int n = 4;
    		if (N < 4)
      			n = N;

    		getCoef( x, n, c );

    		if (n == 4) {
      			point = c[0]*d[0] + c[1]*d[1] + c[2]*d[2] + c[3]*d[3]; 
    		}
    		else if (n == 2) {
      			point = c[0]*d[0] + c[1]*d[1];
    		}

    		return point;
  		} // interpPoint

    }
}

#else
/**
   \file

   This file contains code to implement Lifting Scheme wavelets.
   Lifting Scheme wavelets are described in Wim Sweldens' tutorial
   paper <i>Building Your Own Wavelets at Home</i> which is available
   on the Web.

   Lifting Scheme wavelets are a conceptual extension of Haar wavelets.
   One of the disadvantages of Haar wavelets is the high frequency
   (largest) coefficient spectrum can miss detail (even to odd
   transitions, for example).  Lifting Scheme wavelets properly
   represent change in all coefficient spectrum.  This makes lifting
   scheme wavelets a better choice for some algorithms which do
   filtering based on the absolute standard deviation calculated on the
   high frequency coefficients.

 */


/**
   Support for four point polynomial interpolation, using the Lagrange
   formula.

 */
class polyinterp {

   private:
   double fourPointTable[4][4];
   double twoPointTable[4][4];

   /**
     The polynomial interpolation algorithm assumes that the known
     points are located at x-coordinates 0, 1,.. N-1.  An interpolated
     point is calculated at <b><i>x</i></b>, using N coefficients.  The
     polynomial coefficients for the point <b><i>x</i></b> can be 
     calculated staticly, using the Lagrange method.

     @param x the x-coordinate of the interpolated point
     @param N the number of polynomial points.
     @param c[] an array for returning the coefficients

     */
  void lagrange( double x, int N, double c[] )
  {
    double num, denom;

    for (int i = 0; i < N; i++) {
      num = 1;
      denom = 1;
      for (int k = 0; k < N; k++) {
	if (i != k) {
	  num = num * (x - k);
	  denom = denom * (i - k);
	}
      } // for k
      c[i] = num / denom;
    } // for i
  } // lagrange


  /**
    For a given N-point polynomial interpolation, fill the coefficient
    table, for points 0.5 ... (N-0.5).

   */
  void fillTable( const int N, double table[4][4] )
  {
    double x;
    double n = N;
    int i = 0;

    for (x = 0.5; x < n; x = x + 1.0) {
      lagrange( x, N, table[i] );
      i++;
    }
  } // fillTable


  /**
    Print an N x N table polynomial coefficient table
   */
  void printTable( double table[4][4], int N)
  {
    Console.Write("%d-point interpolation table:\n", N);
    double x = 0.5;
    for (int i = 0; i < N; i++) {
      Console.Write("%4.2f: ", x);
      for (int j = 0; j < N; j++) {
	Console.Write("%6.4f", table[i][j] );
	if (j < N-1)
	  Console.Write(", ");
      }
      Coneole.Write("\n");
      x = x + 1.0;
    }    
  }


  /**
    Print the 4-point and 2-point polynomial coefficient
    tables.
   */
  void printTables()
  {
    printTable( fourPointTable, 4 );
    printTable( twoPointTable, 2 );
    Console.Write("\n");
  } // printTables


    /**
    For the polynomial interpolation point x-coordinate 
    <b><i>x</i></b>, return the associated polynomial
    interpolation coefficients.

    @param x the x-coordinate for the interpolated pont
    @param n the number of polynomial interpolation points
    @param c[] an array to return the polynomial coefficients
   */
  void getCoef( double x, int n, double c[])
  {
    int j = (int)x;

    if (j < 0 || j >= n) {
      Console.Write("poly::getCoef: n = %d, bad x value = %f\n", n, x);
    }

    if (n == 2) {
      c[2] = 0.0;
      c[3] = 0.0;
    }
    else if (n != 4) {
      Console.Write("poly::getCoef: bad value for N\n");
      return;
    }

    for (int i = 0; i < n; i++) {
      if (n == 4) {
	c[i] = fourPointTable[j][i];
      }
      else { // n == 2
	c[i] = twoPointTable[j][i];
      }
    }

  } // getCoef


   public:

   /**
      The polyinterp (polynomial interpolation) class constructor
      calculates the coefficients for a four point polynomial
      interpolated at -0.5, 0.5, 1.5, 2.5 and 3.5.
   */
   polyinterp() 
   {
    // Fill in the 4-point polynomial interplation table
    // for the points 0.5, 1.5, 2.5, 3.5
    fillTable( 4, fourPointTable );
    
    // Fill in the 2-point polynomial interpolation table
    // for 0.5 and 1.5
    fillTable( 2, twoPointTable );
   } // polyinterp class constructor


  /**
    Given four points at the x,y coordinates {0,d<sub>0</sub>},
    {1,d<sub>1</sub>}, {2,d<sub>2</sub>}, {3,d<sub>3</sub>}
    return the y-coordinate value for the polynomial interpolated
    point at <b><i>x</i></b>.
    
    @param x the x-coordinate for the point to be interpolated
    @param N the number of interpolation points
    @param d[] an array containing the y-coordinate values for
           the known points (which are located at x-coordinates
	   0..N-1).
   */
  double interpPoint( double x, int N, double d[])
  {
    double c[4];
    double point = 0;
    
    int n = 4;
    if (N < 4)
      n = N;

    getCoef( x, n, c );

    if (n == 4) {
      point = c[0]*d[0] + c[1]*d[1] + c[2]*d[2] + c[3]*d[3]; 
    }
    else if (n == 2) {
      point = c[0]*d[0] + c[1]*d[1];
    }

    return point;
  } // interpPoint

};  // class interpolation


#endif