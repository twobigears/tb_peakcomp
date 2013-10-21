/*
The MIT License (MIT)

Copyright (c) 2013 Two Big Ears Ltd.
http://twobigears.com

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/* 
Stereo L+R peak compressor with variable knee-smooth, attack, release and makeup gain.
Varun Nair. varun@twobigears.com

Inspired by:
http://www.eecs.qmul.ac.uk/~josh/documents/GiannoulisMassbergReiss-dynamicrangecompression-JAES2012.pdf
https://ccrma.stanford.edu/~jos/filters/Nonlinear_Filter_Example_Dynamic.htm
*/

#include "m_pd.h" 
#include <math.h>
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

 static t_class *peakcomp_tilde_class; 

 typedef struct _peakcomp_tilde{

	t_object	x_obj;				// Pointer to object 
	double		f_float2Sig;		// Dummy variable for conversion of floats in inlet to signals
	float		fFs;				// Sample rate
	
	// Arrays for delay lines 
	float 		kneeRecursive[2];
	float 		attackRecursive[2];
	float 		releaseRecursive[2];	

	// Values
	float 		attack;
	float 		release;
	float 		ratio;
	float 		threshold;
	float 		knee;
	float 		kneeCoeffs;
	float 		kneeCoeffsMinus;
	float 		attackCoeffs;
	float 		attackCoeffsMinus;
	float 		releaseCoeff;
	float 		releaseCoeffMinus;
	float 		onebyfFS;
	float 		makeupGain;

	
} t_peakcomp_tilde; 

// Forward declarations

// DSP processing
t_int *peakcomp_tilde_perform(t_int *w);
// DSP setup 
void peakcomp_tilde_dsp(t_peakcomp_tilde *x, t_signal **sp);
// Update comp controls
void peakcomp_tilde_update(t_peakcomp_tilde *x);
// Constructor
void *peakcomp_tilde_new(t_floatarg th, t_floatarg ra, t_floatarg at, t_floatarg re, t_floatarg mg);
// Class init
void tb_peakcomp_tilde_setup(void);
// Threshold message in
void peakcomp_tilde_threshold(t_peakcomp_tilde *x, t_floatarg thresh);
// Ratio message in
void peakcomp_tilde_ratio(t_peakcomp_tilde *x, t_floatarg ratio);
// Attack message in
void peakcomp_tilde_attack(t_peakcomp_tilde *x, t_floatarg attack);
// Release message in
void peakcomp_tilde_release(t_peakcomp_tilde *x, t_floatarg release);
// Makeup gain message in
void peakcomp_tilde_makeup(t_peakcomp_tilde *x, t_floatarg makeup);
// Knee smooth message in
void peakcomp_tilde_knee(t_peakcomp_tilde *x, t_floatarg knee);
// Gain reduction outout toggle message in
void peakcomp_tilde_reductionOut(t_peakcomp_tilde *x, t_floatarg knee);

t_int *peakcomp_tilde_perform(t_int *w) {

	t_peakcomp_tilde *x = (t_peakcomp_tilde *)(w[1]);	// Object ref
	t_sample *in1     = (t_sample *)(w[2]);				// Input samples (left)
	t_sample *in2     = (t_sample *)(w[3]);				// Input samples (right)
	t_sample *out1    = (t_sample *)(w[4]);				// Output samples (left)
	t_sample *out2	  = (t_sample *)(w[5]);		 		// Output samples (right)
	t_sample *outGR	  = (t_sample *)(w[6]);		 		// Output samples (gain reduction)
	int		 n		  = (int)		(w[7]);      		// Vector size
	
	// DSP loop
	for (int i=0; i<n; i++) {

		float input1 = (float)in1[i];
		float input2 = (float)in2[i];

		// L+R sum
		float peakEnv = fabsf(input1) + fabsf(input2);
		//release recursive
		x->releaseRecursive[0] = (x->releaseCoeffMinus * peakEnv) + (x->releaseCoeff * MAX(peakEnv, x->releaseRecursive[1]));
		//attack recursive
		x->attackRecursive[0] = ((x->attackCoeffsMinus * x->releaseRecursive[0]) + (x->attackCoeffs * x->attackRecursive[1]));	
		//knee smoothening and gain reduction
		x->kneeRecursive[0] = (x->kneeCoeffsMinus * MAX(MIN(((x->threshold + (x->ratio * (x->attackRecursive[0] - x->threshold))) / x->attackRecursive[0]), 1.f), 0.f)) + (x->kneeCoeffs * x->kneeRecursive[1]);

		out1[i] = input1 * x->kneeRecursive[0] * x->makeupGain;
		out2[i] = input2 * x->kneeRecursive[0] * x->makeupGain;
		outGR[i] = x->kneeRecursive[0];

		x->releaseRecursive[1] = x->releaseRecursive[0];
		x->attackRecursive[1] = x->attackRecursive[0];
		x->kneeRecursive[1] = x->kneeRecursive[0];

	}

	return (w + 8); 
}

void peakcomp_tilde_dsp(t_peakcomp_tilde *x, t_signal **sp) 
{
	// Initialise paramters

	for (int i = 0; i < 2; i++)
	{
		x->kneeRecursive[i] = 0.f;
		x->attackRecursive[i] = 0.f;
		x->releaseRecursive[i] = 0.f;
	}

	// Get sample rate
	float fs = sp[0]->s_sr; 

	// Sample rate change
	if(x->fFs != fs) 
	{
		x->fFs = sp[0]->s_sr;
		x->onebyfFS = 1.f / x->fFs;
	}

	// Add perform() function to dsp graph. 
	dsp_add(peakcomp_tilde_perform, 
		7,					  
		x,					  
			sp[0]->s_vec,         // Inlet1 vector
			sp[1]->s_vec,		  // Inlet2 vector
			sp[2]->s_vec,		  // Outlet1 vector
			sp[3]->s_vec,		  // Outlet2 vector
			sp[4]->s_vec,		  // Outlet3 vector
			sp[0]->s_n			  // Vector size
			) ;

}

void peakcomp_tilde_update(t_peakcomp_tilde *x) {

	x->kneeCoeffs = expf(0.f - (x->onebyfFS / x->knee));
	x->kneeCoeffsMinus = 1.f - x->kneeCoeffs;
	x->attackCoeffs = expf(0.f - (x->onebyfFS / x->attack));
	x->attackCoeffsMinus = 1.f - x->attackCoeffs;
	x->releaseCoeff = expf(0.f - (x->onebyfFS / x->release));
	x->releaseCoeffMinus = 1.f - x->releaseCoeff;
}

void peakcomp_tilde_threshold(t_peakcomp_tilde *x, t_floatarg thresh) {
	if (thresh <= 0) {
		// dB to linear (could use the function from m_pd.h)
		x->threshold = pow(10, (thresh*0.05f));	
	}
	else x->threshold = 0;
}

void peakcomp_tilde_ratio(t_peakcomp_tilde *x, t_floatarg ratio) {
	if (ratio >= 1) {
		x->ratio = 1.f/ratio;
	}
	else x->ratio = 1;
}

void peakcomp_tilde_attack(t_peakcomp_tilde *x, t_floatarg attack) {
	if (attack >= 0.001) {
		x->attack = attack * 0.001;
	}
	else x->attack = 0.000001;
	peakcomp_tilde_update(x);
}

void peakcomp_tilde_release(t_peakcomp_tilde *x, t_floatarg release) {
	if (release >= 0.001) {
		x->release = release * 0.001;
	}
	else x->release = 0.000001;
	peakcomp_tilde_update(x);
}

void peakcomp_tilde_makeup(t_peakcomp_tilde *x, t_floatarg makeup) {
	// dB to linear (could use the function from m_pd.h)
	x->makeupGain = pow(10, (makeup * 0.05));	
}

void peakcomp_tilde_knee(t_peakcomp_tilde *x, t_floatarg knee) {
	if (knee >= 0 && knee <= 1)
	{	
		// knee value (0 to 1) is scaled from 0 (hard) to 0.02 (smooth). Could be scaled to a larger number.
		x->knee = knee * 0.02;
		peakcomp_tilde_update(x);
	}	
}

// Constructor
void *peakcomp_tilde_new(t_floatarg th, t_floatarg ra, t_floatarg at, t_floatarg re, t_floatarg mg)
{	

	post("tb_peakcomp~ v1.0 2013 Two Big Ears");

	//Make a new instance of peakcomp~
	t_peakcomp_tilde *x = (t_peakcomp_tilde *) pd_new(peakcomp_tilde_class); 

	//Values from object arguments
	peakcomp_tilde_threshold(x, th);
	peakcomp_tilde_ratio(x, ra);
	peakcomp_tilde_attack(x, at);
	peakcomp_tilde_release(x, re);
	peakcomp_tilde_makeup(x, mg);

	//Make new signal inlets and outlets
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
	outlet_new(&x->x_obj, &s_signal);
	outlet_new(&x->x_obj, &s_signal);
	outlet_new(&x->x_obj, &s_signal);

	x->fFs = sys_getsr(); 

	x->onebyfFS = 1.f / x->fFs;

	peakcomp_tilde_update(x);

	//cast back to void* to return pointer to Pd
	return (void *) x; 
}

void tb_peakcomp_tilde_setup(void) {

	peakcomp_tilde_class = class_new(gensym("tb_peakcomp~"),				
							(t_newmethod)peakcomp_tilde_new,				
							0,												
							sizeof(t_peakcomp_tilde),					
							CLASS_DEFAULT,									
							A_DEFFLOAT, A_DEFFLOAT, A_DEFFLOAT, A_DEFFLOAT, A_DEFFLOAT,	// Threshold (dB), ratio, attack(ms), release (ms), makeup gain (dB)
							0);												

	class_addmethod(peakcomp_tilde_class, 
		(t_method)peakcomp_tilde_dsp, gensym("dsp"), 0); 
	
	// Messages as input
	class_addmethod(peakcomp_tilde_class,
		(t_method)peakcomp_tilde_threshold, gensym("threshold"),
		A_DEFFLOAT, 0);

	class_addmethod(peakcomp_tilde_class,
		(t_method)peakcomp_tilde_ratio, gensym("ratio"),
		A_DEFFLOAT, 0);

	class_addmethod(peakcomp_tilde_class,
		(t_method)peakcomp_tilde_attack, gensym("attack"),
		A_DEFFLOAT, 0);

	class_addmethod(peakcomp_tilde_class,
		(t_method)peakcomp_tilde_release, gensym("release"),
		A_DEFFLOAT, 0);

	class_addmethod(peakcomp_tilde_class,
		(t_method)peakcomp_tilde_makeup, gensym("makeup"),
		A_DEFFLOAT, 0);

	class_addmethod(peakcomp_tilde_class,
		(t_method)peakcomp_tilde_knee, gensym("knee"),
		A_DEFFLOAT, 0);
	
	// Enables 'floats' as signals on first inlet
	CLASS_MAINSIGNALIN(peakcomp_tilde_class, t_peakcomp_tilde, f_float2Sig);
}
