/*
 * Created by SharpDevelop.
 * User: jeffA
 * Date: 12/7/2015
 * Time: 4:17 PM
 * 
 * Filename: LiftBase.cs
 * This is part of a wavelet library developed for SDR# for use by SDR# plugins.
 * This is simply a C# port of a wavelet library developed by Ian Kaplan.  His original
 * source code has been included at the branch of the ifdef directive.
 * This code supports higher level wavelet transform functions.
 */
 
#define CSHARPCODE

using System;
using System.Diagnostics;

#if CSHARPCODE
namespace SDRSharp.WaveletLib
{
	/// <summary>
	/// Description of LiftBase.
	/// </summary>
	public class LiftBase
	{
		//public LiftBase()
		//{
		//}
		
		public enum transDirection{
			/** "enumeration" for forward wavelet transform */
    		forward = 1,
  			/** "enumeration" for inverse wavelet transform */
    		inverse = 2
		};

  /**
    Split the <i>vec</i> into even and odd elements,
    where the even elements are in the first half
    of the vector and the odd elements are in the
    second half.
   */
  		public void split( double[] vec, int N )
  		{
    
    		int start = 1;
    		int end = N - 1;

    		while (start < end) {
      			for (int i = start; i < end; i = i + 2) {
					double tmp = vec[i];
					vec[i] = vec[i+1];
					vec[i+1] = tmp;
      			}
      			start = start + 1;
      			end = end - 1;
    		}
  		}

  /**
    Merge the odd elements from the second half of the N element
    region in the array with the even elements in the first
    half of the N element region.  The result will be the
    combination of the odd and even elements in a region
    of length N.
    
   */
  		public void merge( double[] vec, int N )
  		{
    		int half = N >> 1;
    		int start = half-1;
    		int end = half;
    
    		while (start > 0) {
      			for (int i = start; i < end; i = i + 2) {
					double tmp = vec[i];
					vec[i] = vec[i+1];
					vec[i+1] = tmp;
      			}
      			start = start - 1;
      			end = end + 1;
    		}
  		}


  /** 
    Predict step, to be defined by the subclass

    @param vec input array
    @param N size of region to act on (from 0..N-1)
    @param direction forward or inverse transform

   */
  		protected extern virtual void predict( double[] vec, int N, transDirection direction ); //edit: ; repby =0;

  /**
    Reverse predict step.  

    The predict step applied the high pass filter to the data
    set and places the result in the upper half of the array.
    The reverse predict step applies the high pass filter and
    places the result in the lower half of the array.

    This reverse predict step is only used by wavelet packet
    frequency analysis algorithms.  The default version
    of this algorihtm does nothing.
   */
  		protected extern virtual void predictRev( double[] vec, int N, transDirection direction );  //edit: ; repby {};
  

  /** 
    Update step, to be defined by the subclass 

    @param vec input array
    @param N size of region to act on (from 0..N-1)
    @param direction forward or inverse transform

  */
 		protected extern virtual void update( double[] vec, int N, transDirection direction );  //edit: ; repby =0;


  /**
    Reverse update step
   */
  		protected void updateRev( double[] vec, int N, transDirection direction ) {}



  /**
    One step in the forward wavelet transform
   */
  		public virtual void forwardStep( double[] vec, int n )  //edit: const int n repby int n
  		{
  			split(vec,n);
  			predict(vec,n,transDirection.forward);
  			update(vec,n,transDirection.forward);	//edit: forwardStep repby forward
  		} // forwardStep

  /**
    Reverse forward transform step.  The result of the high
    pass filter is stored in the lower half of the array
    and the result of the low pass filter is stored in the
    upper half.

    This function should be defined by any subclass that
    is used for wavelet frequency analysis.
   */
  		public virtual void forwardStepRev( double[] vec, int N )  //edit: const int n repby int n
  		{
  			Debug.Assert(false);
  		}

  /**
    Simple wavelet Lifting Scheme forward transform

     forwardTrans is passed an indexable object.  The object must
    contain a power of two number of data elements.  Lifting Scheme
    wavelet transforms are calculated in-place and the result is
    returned in the argument array.

    The result of forwardTrans is a set of wavelet coefficients
    ordered by increasing frequency and an approximate average
    of the input data set in vec[0].  The coefficient bands
    follow this element in powers of two (e.g., 1, 2, 4, 8...).

  */
 		public virtual void forwardTrans( double[] vec, int N )  //edit: const int n repby int n
  		{

    		for (int n = N; n > 1; n = n >> 1) {
      			forwardStep( vec, n );
    		}
  		} // forwardTrans


  /**
    One inverse wavelet transform step
   */
  		public virtual void inverseStep( double[] vec, int n )  //edit: const int n repby int n
  		{
    		update( vec, n, transDirection.inverse );
    		predict( vec, n, transDirection.inverse );
    		merge( vec, n );
  		}

  /** 
    Reverse inverse transform step.  Calculate the inverse transform
    from a high pass filter result stored in the lower half of the
    array and a low pass filter result stored in the upper half.

    This function should be defined by any subclass that
    is used for wavelet frequency analysis.
   */
  		public virtual void inverseStepRev( double[] vec, int n )  //edit: const int n repby int n
  		{
  			Debug.Assert(false);
  		}


  /**
    Default two step Lifting Scheme inverse wavelet transform

    inverseTrans is passed the result of an ordered wavelet 
    transform, consisting of an average and a set of wavelet
    coefficients.  The inverse transform is calculated
    in-place and the result is returned in the argument array.

   */
  		public virtual void inverseTrans( double[] vec, int N )  //edit: const int n repby int n
  		{

    		for (int n = 2; n <= N; n = n << 1) {
      			inverseStep( vec, n );
    		}
  		} // inverseTrans

	}
}
#else
protected:

  typedef enum { 
  /** "enumeration" for forward wavelet transform */
    forward = 1,
  /** "enumeration" for inverse wavelet transform */
    inverse = 2 
  } transDirection;


  /**
    Split the <i>vec</i> into even and odd elements,
    where the even elements are in the first half
    of the vector and the odd elements are in the
    second half.
   */
  void split( T& vec, int N )
  {
    
    int start = 1;
    int end = N - 1;

    while (start < end) {
      for (int i = start; i < end; i = i + 2) {
	T_elem tmp = vec[i];
	vec[i] = vec[i+1];
	vec[i+1] = tmp;
      }
      start = start + 1;
      end = end - 1;
    }
  }

  /**
    Merge the odd elements from the second half of the N element
    region in the array with the even elements in the first
    half of the N element region.  The result will be the
    combination of the odd and even elements in a region
    of length N.
    
   */
  void merge( T& vec, int N )
  {
    int half = N >> 1;
    int start = half-1;
    int end = half;
    
    while (start > 0) {
      for (int i = start; i < end; i = i + 2) {
	T_elem tmp = vec[i];
	vec[i] = vec[i+1];
	vec[i+1] = tmp;
      }
      start = start - 1;
      end = end + 1;
    }
  }


  /** 
    Predict step, to be defined by the subclass

    @param vec input array
    @param N size of region to act on (from 0..N-1)
    @param direction forward or inverse transform

   */
  virtual void predict( T& vec, int N, transDirection direction ) = 0;

  /**
    Reverse predict step.  

    The predict step applied the high pass filter to the data
    set and places the result in the upper half of the array.
    The reverse predict step applies the high pass filter and
    places the result in the lower half of the array.

    This reverse predict step is only used by wavelet packet
    frequency analysis algorithms.  The default version
    of this algorihtm does nothing.
   */
  virtual void predictRev( T& vec, int N, transDirection direction ) {};
  

  /** 
    Update step, to be defined by the subclass 

    @param vec input array
    @param N size of region to act on (from 0..N-1)
    @param direction forward or inverse transform

  */
  virtual void update( T& vec, int N, transDirection direction ) = 0;


  /**
    Reverse update step
   */
  virtual void updateRev( T& vec, int N, transDirection direction ) {}

public:

  /**
    One step in the forward wavelet transform
   */
  virtual void forwardStep( T& vec, const int n )
  {
    split( vec, n );
    predict( vec, n, forward );
    update( vec, n, forward );
  } // forwardStep

  /**
    Reverse forward transform step.  The result of the high
    pass filter is stored in the lower half of the array
    and the result of the low pass filter is stored in the
    upper half.

    This function should be defined by any subclass that
    is used for wavelet frequency analysis.
   */
  virtual void forwardStepRev( T& vec, const int N )
  {
    assert(false);
  }

  /**
    Simple wavelet Lifting Scheme forward transform

     forwardTrans is passed an indexable object.  The object must
    contain a power of two number of data elements.  Lifting Scheme
    wavelet transforms are calculated in-place and the result is
    returned in the argument array.

    The result of forwardTrans is a set of wavelet coefficients
    ordered by increasing frequency and an approximate average
    of the input data set in vec[0].  The coefficient bands
    follow this element in powers of two (e.g., 1, 2, 4, 8...).

  */
  virtual void forwardTrans( T& vec, const int N )
  {

    for (int n = N; n > 1; n = n >> 1) {
      forwardStep( vec, n );
    }
  } // forwardTrans


  /**
    One inverse wavelet transform step
   */
  virtual void inverseStep( T& vec, const int n )
  {
    update( vec, n, inverse );
    predict( vec, n, inverse );
    merge( vec, n );
  }

  /** 
    Reverse inverse transform step.  Calculate the inverse transform
    from a high pass filter result stored in the lower half of the
    array and a low pass filter result stored in the upper half.

    This function should be defined by any subclass that
    is used for wavelet frequency analysis.
   */
  virtual void inverseStepRev( T& vec, const int n )
  {
    assert( false );
  }


  /**
    Default two step Lifting Scheme inverse wavelet transform

    inverseTrans is passed the result of an ordered wavelet 
    transform, consisting of an average and a set of wavelet
    coefficients.  The inverse transform is calculated
    in-place and the result is returned in the argument array.

   */
  virtual void inverseTrans( T& vec, const int N )
  {

    for (int n = 2; n <= N; n = n << 1) {
      inverseStep( vec, n );
    }
  } // inverseTrans

#endif
