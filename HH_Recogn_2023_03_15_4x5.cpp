#define _CRT_SECURE_NO_WARNINGS

#define _USE_MATH_DEFINES
#include "math.h"
#include <stdlib.h>
#include <stdio.h>
#include <locale.h>
#include <time.h>
#include <stdbool.h>
#include <list>
#include <omp.h>
#include <vector>

using namespace std;

#define Neuron_width 4
#define Neuron_height 5

#define Neuron_layer (int)(Neuron_width * Neuron_height)
#define Neuron_count (int)(1 + 3 * (Neuron_width * Neuron_height))

#define Equations_per_neuron 4

#define Neuron_equations_count Equations_per_neuron * Neuron_count
#define Equations_count Neuron_equations_count

double *f;
double *f_diff;

double **k;
double *phi_k1;
double *phi_k2;
double *phi_k3;

bool write_file = false;
bool omp_enable = false;

double I_syn_buffer[Neuron_count];

double g_syn_RN_param;
double g_syn_SN_param;
double g_syn_IN_param;
double g_syn_CN_param;

double I_app_RN_param;
double I_app_SN_param;
double I_app_IN_param;
double I_app_CN_param;

// https://neuronaldynamics.epfl.ch/online/Ch2.S2.html
double C_m = 1;	  // muF/cm^2
double g_K = 35;  // mS/cm^2
double g_Na = 40; // mS/cm^2
double g_L = 0.3; // mS/cm^2
double E_K = -77; // mV
double E_Na = 55; // mV
double E_L = -65; // mV

double p_rewir;
double p_inhib;

double I_app_min;
double I_app_max;

double *g_syn;
double k_syn = 0.2;
double *E_syn;

double *I_app;

double **N_A;
double **N_B;
double *N_C;

list<double> *V_spikes;
list<double> *V_spikes_Freq;

FILE *fp_I_syn;

double *Meander_start_from_zero;
double *Meander_width;
double *Meander_height;
double *Meander_interval;
int *Meander_count;

const double Offset = 0.00;
const double Duration = 1000.0;
const double Interval = 100.0;
const double Magnitude = 0.0;
const int Count = 1000;

const double t_start = 0;
const double t_max = 10.0; // 100 msec = 0.1 sec
const double dt = 0.00005; // 0.01 msec = 0.00001 sec; 0.1 msec = 0.0001 sec; 1 msec = 0.001 sec

double **tau;
double** E_syn_matrix;

#define tau_lag 20 // ms

#define ms_to_step (int)(0.001 / dt) // (0.001 / dt) !!! Don't forget !!!

#define Max_delay (int)(tau_lag * ms_to_step)
double **V_old_array;

int sign(int value)
{
	if (value >= 0)
		return 1;
	else if (value < 0)
		return -1;

	//return 0;
}

// bool thread_count_printed = false;

double I_stim(int i, double t)
{
	if (t < Meander_start_from_zero[i])
		return 0;

	int count = t / (Meander_width[i] + Meander_interval[i]);

	if (count > Meander_count[i])
		return 0;

	t -= Meander_start_from_zero[i];
	t = fmod(t, Meander_width[i] + Meander_interval[i]);

	return t < Meander_width[i] ? Meander_height[i] : 0;
}

double V(int i)
{
	return f[i * Equations_per_neuron];
}

void SetV(int i, double value)
{
	f[i * Equations_per_neuron] = value;
}

double m(int i)
{
	return f[i * Equations_per_neuron + 1];
}

void Setm(int i, double value)
{
	f[i * Equations_per_neuron + 1] = value;
}

double n(int i)
{
	return f[i * Equations_per_neuron + 2];
}

void Setn(int i, double value)
{
	f[i * Equations_per_neuron + 2] = value;
}

double h(int i)
{
	return f[i * Equations_per_neuron + 3];
}

void Seth(int i, double value)
{
	f[i * Equations_per_neuron + 3] = value;
}

///

double V_old(int i, int delay)
{
	if (delay == 0)
		return V(i);

	return V_old_array[i][Max_delay - delay];
}

int RandomI(int min, int max)
{
	return ((double)rand() / (RAND_MAX - 1)) * (max - min) + min;
}

double RandomD(double min, double max)
{
	return ((double)rand() / RAND_MAX) * (max - min) + min;
}

double alpha_m(double *f, int i)
{
	return 0.182 * (V(i) + 35) / (1 - exp(-(V(i) + 35) / 9));
}

double beta_m(double *f, int i)
{
	return -0.124 * (V(i) + 35) / (1 - exp((V(i) + 35) / 9));
}

double alpha_n(double *f, int i)
{
	return 0.02 * (V(i) - 25) / (1 - exp(-(V(i) - 25) / 9));
}

double beta_n(double *f, int i)
{
	return -0.002 * (V(i) - 25) / (1 - exp((V(i) - 25) / 9));
}

double alpha_h(double *f, int i)
{
	return 0.25 * exp(-(V(i) + 90) / 12);
}

double beta_h(double *f, int i)
{
	return 0.25 * exp((V(i) + 62) / 6) / exp((V(i) + 90) / 12);
}

//

double HodgkinHuxley(int i, double *f, double t)
{
	int in = i / Equations_per_neuron;
	int il = i % Equations_per_neuron;

	switch (il)
	{
		case 0: // V
		{
			double I_syn = 0;

			for (int j = 0; j < N_C[in]; j++)
			{
				int input = N_B[in][j];
				I_syn += g_syn[input] * (E_syn_matrix[input][in] - V(in)) / (1 + exp(-(V_old(input, tau[input][in]) / k_syn)));
			}

			I_syn_buffer[in] = I_syn;

			return 1000 * ((g_Na * pow(m(in), 3) * h(in) * (E_Na - V(in)) + g_K * n(in) * (E_K - V(in)) + g_L * (E_L - V(in)) + I_app[in] + I_syn + I_stim(in, t)) / C_m); // V
		}

		case 1: // m
		{
			return 1000 * (alpha_m(f, in) * (1 - m(in)) - beta_m(f, in) * m(in)); // m
		}

		case 2: // n
		{
			return 1000 * (alpha_n(f, in) * (1 - n(in)) - beta_n(f, in) * n(in)); // n
		}

		case 3: // h
		{
			return 1000 * (alpha_h(f, in) * (1 - h(in)) - beta_h(f, in) * h(in)); // h
		}
	}

	return 0;
}

void Euler(double t, double dt, double *f, double *f_next)
{
	#pragma omp parallel for
	for (int i = 0; i < Equations_count; i++)
	{
		f_next[i] = f[i] + HodgkinHuxley(i, f, t) * dt;
	}
}

void RungeKutta(double t, double dt, double *f, double *f_next)
{
	#pragma omp parallel for
	for (int i = 0; i < Equations_count; i++)
	{
		// if (!thread_count_printed)
		//{
		//	thread_count_printed = true;
		//	printf("Threads = %d\n", omp_get_num_threads());
		// }

		k[i][0] = HodgkinHuxley(i, f, t) * dt;
		phi_k1[i] = f[i] + k[i][0] / 2;
		k[i][1] = HodgkinHuxley(i, phi_k1, t) * dt;
		phi_k2[i] = f[i] + k[i][1] / 2;
		k[i][2] = HodgkinHuxley(i, phi_k2, t) * dt;
		phi_k3[i] = f[i] + k[i][2] / 2;
		k[i][3] = HodgkinHuxley(i, phi_k3, t) * dt;

		f_next[i] = f[i] + (k[i][0] + 2 * k[i][1] + 2 * k[i][2] + k[i][3]) / 6;
	}
}

void CopyArray(double *source, double *target, int N)
{
	for (int i = 0; i < N; i++)
		target[i] = source[i];
}

bool Approximately(double a, double b)
{
	if (a < 0)
		a = -a;

	if (b < 0)
		b = -b;

	return a - b <= 0.000001;
}

// http://preshing.com/20111007/how-to-generate-random-timings-for-a-poisson-process/
double nextTime(double rateParameter)
{
	return -log(1.0 - (double)rand() / (RAND_MAX)) / rateParameter;
}

void GenerateMeander()
{
	for (int i = 0; i < Neuron_count; i++)
	{
		if (i == 0)
		{
			Meander_start_from_zero[i] = Offset;
			Meander_interval[i] = Interval;
			Meander_width[i] = Duration;
			Meander_height[i] = Magnitude;
			Meander_count[i] = Count;
		}
		else
		{
			Meander_start_from_zero[i] = 0;
			Meander_interval[i] = 0;
			Meander_width[i] = 0;
			Meander_height[i] = 0;
			Meander_count[i] = 0;
		}
	}
}

void FillMatrixZero()
{
	for (int i = 0; i < Neuron_count; i++)
	{
		for (int j = 0; j < Neuron_count; j++)
		{
			N_A[i][j] = 0;
			E_syn_matrix[i][j] = -90;
		}
	}
}

void Fill_N_BCMatrix()
{
	for (int i = 0; i < Neuron_count; i++)
	{
		int bIndex = 0;
		N_C[i] = 0;

		for (int j = 0; j < Neuron_count; j++)
		{
			if (N_A[j][i] == 1)
			{
				N_B[i][bIndex] = j;
				bIndex++;
				N_C[i]++;
			}
		}
	}
}

bool CheckNeuronSameLine(int i, int j)
{
	return i / Neuron_width == j / Neuron_width;
}

bool IsNeuronWireNeighbors(int i, int j)
{
	if (CheckNeuronSameLine(i, j) && (i == j - 1 || i == j + 1))
		return true;

	if (i == j - Neuron_width || i == j + Neuron_width)
		return true;

	return false;
}

void FillNeuronMatrix()
{
	// 0
	for (int i = 1; i < 2 * Neuron_layer + 1; i++)
		N_A[0][i] = 1;

	int first_start = 1;
	int first_end = Neuron_layer + 1;

	for (int i = first_start; i < first_end; i++)
		N_A[i][i + Neuron_layer] = 1;

	int second_start = first_end;
	int second_end = second_start + Neuron_layer;

	for (int i = second_start; i < second_end; i++)
	{
		for (int j = Neuron_layer; j < 2 * Neuron_layer; j++)
		{
			int a = second_start + j;
			N_A[i][second_start + j] = 1;
		}
	}
}

void FillgsynArray()
{
	// 0
	g_syn[0] = g_syn_RN_param;

	int first_start = 1;
	int first_end = Neuron_layer + 1;

	for (int i = first_start; i < first_end; i++)
		g_syn[i] = g_syn_SN_param;

	int second_start = first_end;
	int second_end = second_start + Neuron_layer;

	for (int i = second_start; i < second_end; i++)
	{
		g_syn[i] = g_syn_IN_param;
	}

	int third_start = second_end;
	int third_end = third_start + Neuron_layer;

	for (int i = third_start; i < third_end; i++)
	{
		g_syn[i] = g_syn_CN_param;
	}
}

void FillIappArray()
{
	// 0
	I_app[0] = I_app_RN_param;

	int first_start = 1;
	int first_end = Neuron_layer + 1;

	for (int i = first_start; i < first_end; i++)
		I_app[i] = I_app_SN_param;

	int second_start = first_end;
	int second_end = second_start + Neuron_layer;

	for (int i = second_start; i < second_end; i++)
	{
		I_app[i] = I_app_IN_param;
	}

	int third_start = second_end;
	int third_end = third_start + Neuron_layer;

	for (int i = third_start; i < third_end; i++)
	{
		I_app[i] = I_app_CN_param;
	}
}

void FillEsynMatrix()
{
	// 0
	for (int i = 1; i < 2 * Neuron_layer + 1; i++)
		E_syn_matrix[0][i] = 0;

	// 1
	FILE* fpSN;
	fpSN = fopen("SN_input_array.txt", "r");
	int sn[Neuron_layer];

	for (int i = 0; i < Neuron_layer; i++)
 	{
  		fscanf(fpSN, "%d", &sn[i]);
 	}
 	fclose(fpSN);

	int first_start = 1;
	int first_end = Neuron_layer + 1;
	int sn_index = 0;

	for (int i = first_start; i < first_end; i++)
	{
		E_syn_matrix[i][i + Neuron_layer] = sn[sn_index] == 1 ? 0 : -90;
		sn_index++;
	}

	// 2
	fpSN = fopen("Patterns.txt", "r");
	vector<double*> patterns;
	double* current_pattern;
	int index = 0;

	while (!feof(fpSN))
	{
		if (index == 0)
			current_pattern = new double[Neuron_layer];

		int value = 0;
		fscanf(fpSN, "%d", &value);

		current_pattern[index] = value;
		index++;

		if (index == Neuron_layer)
		{
			index = 0;
			patterns.push_back(current_pattern);
		}
	}

 	fclose(fpSN);

	int second_start = first_end;
	int second_end = second_start + Neuron_layer;

	// Direct
	double a = 45;
	double b = -45;
	// Inversion
	//double a = -45;
	//double b = -45;
	

	for (int i = second_start; i < second_end; i++)
	{
		for (int j = Neuron_layer; j < 2 * Neuron_layer; j++)
		{
			int sum = 0;

			for (int k = 0; k < patterns.size(); k++)
				sum += patterns[k][i - Neuron_layer - 1] * patterns[k][second_start + j - 2 * Neuron_layer - 1];

			E_syn_matrix[i][second_start + j] = a * sign(sum) + b;
		}
	}
}

void FillVOldFromCurrent()
{
	if (Max_delay == 0)
		return;

	for (int i = 0; i < Neuron_count; i++)
		for (int j = 0; j < Max_delay; j++)
			// V_old_array[i][j] = V(i);
			V_old_array[i][j] = -60;
}

void UpdateVOld()
{
	if (Max_delay == 0)
		return;

	for (int i = 0; i < Neuron_count; i++)
	{
		for (int j = 1; j < Max_delay; j++)
			V_old_array[i][j - 1] = V_old_array[i][j];

		V_old_array[i][Max_delay - 1] = V(i);
	}
}
//
void FillTauMatrix()
{
	for (int i = 0; i < Neuron_count; i++)
	{
		for (int j = 0; j < Neuron_count; j++)
		{
			if (i == j || N_A[i][j] == 0)
			{
				tau[i][j] = 0;
				continue;
			}

			tau[i][j] = tau_lag * ms_to_step;
		}
	}
}

int main(int argc, char *argv[])
{
    int write = 0;
	sscanf(argv[1], "%i", &write);
	write_file = write;
	printf("write_file = %i\n", write_file);

	int omp = 0;
	sscanf(argv[2], "%i", &omp);
	omp_enable = omp;
	printf("omp_enable = %i\n", omp_enable);	

	sscanf(argv[3], "%lf", &g_syn_RN_param);
	printf("g_syn_RN = %f\n", g_syn_RN_param);
    
	sscanf(argv[4], "%lf", &g_syn_SN_param);
	printf("g_syn_SN = %f\n", g_syn_SN_param);
	
	sscanf(argv[5], "%lf", &g_syn_IN_param);
	printf("g_syn_IN = %f\n", g_syn_IN_param);
		
	sscanf(argv[6], "%lf", &g_syn_CN_param);
	printf("g_syn_CN = %f\n", g_syn_CN_param);
	
	////

	sscanf(argv[7], "%lf", &I_app_RN_param);
	printf("I_app_RN = %f\n", I_app_RN_param);
    
	sscanf(argv[8], "%lf", &I_app_SN_param);
	printf("I_app_SN = %f\n", I_app_SN_param);
	
	sscanf(argv[9], "%lf", &I_app_IN_param);
	printf("I_app_IN = %f\n", I_app_IN_param);
		
	sscanf(argv[10], "%lf", &I_app_CN_param);
	printf("I_app_CN = %f\n", I_app_CN_param);
	
	//
	
	double E_syn_RN_param;
	sscanf(argv[11], "%lf", &E_syn_RN_param);
	printf("E_syn_RN = %f\n", E_syn_RN_param);
    
    double E_syn_SN_param;
	sscanf(argv[12], "%lf", &E_syn_SN_param);
	printf("E_syn_SN = %f\n", E_syn_SN_param);

	double E_syn_IN_param;
	sscanf(argv[13], "%lf", &E_syn_IN_param);
	printf("E_syn_IN = %f\n", E_syn_IN_param);
	
	double E_syn_CN_param;
	sscanf(argv[14], "%lf", &E_syn_CN_param);
	printf("E_syn_CN = %f\n", E_syn_CN_param);
	
	//

	f = new double[Equations_count];
	f_diff = new double[Equations_count];

	int a = Neuron_count;

	E_syn = new double[Neuron_count];
	g_syn = new double[Neuron_count];
	I_app = new double[Neuron_count];
	V_spikes = new list<double>[Neuron_count];

	Meander_start_from_zero = new double[Neuron_count];
	Meander_width = new double[Neuron_count];
	Meander_height = new double[Neuron_count];
	Meander_interval = new double[Neuron_count];
	Meander_count = new int[Neuron_count];

	tau = new double *[Neuron_count];
	for (int i = 0; i < Neuron_count; i++)
		tau[i] = new double[Neuron_count];

	V_old_array = new double *[Neuron_count];
	for (int i = 0; i < Neuron_count; i++)
		V_old_array[i] = new double[Max_delay];

	E_syn_matrix = new double* [Neuron_count];
	for (int i = 0; i < Neuron_count; i++)
		E_syn_matrix[i] = new double[Neuron_count];

	FILE *fp0;
	FILE *fp_I_stim = 0;
	FILE *fp_V = 0;
	FILE *fp_V_spikes = 0;
	FILE *fp_Esyn;

	// FILE* fp_res;
	srand(time(NULL));

	// for (int i = 0; i < Node_count; i++)
	//	V_old_length[i] = 0;

	N_A = new double *[Neuron_count];
	for (int i = 0; i < Neuron_count; i++)
		N_A[i] = new double[Neuron_count];

	N_B = new double *[Neuron_count];
	for (int i = 0; i < Neuron_count; i++)
		N_B[i] = new double[Neuron_count];

	N_C = new double[Neuron_count];

	FillMatrixZero();
	FillNeuronMatrix();
	Fill_N_BCMatrix();
	// FillFullTauMatrix();
	FillTauMatrix();
	FillgsynArray();
	FillIappArray();
	FillEsynMatrix();

	fp0 = fopen("N_A.txt", "w+");
	for (int i = 0; i < Neuron_count; i++)
	{
		for (int j = 0; j < Neuron_count; j++)
		{
			fprintf(fp0, "%d\t", (int)N_A[i][j]);
		}
		fprintf(fp0, "\n");
	}
	fclose(fp0);

	fp0 = fopen("tau.txt", "w+");
	for (int i = 0; i < Neuron_count; i++)
	{
		for (int j = 0; j < Neuron_count; j++)
		{
			fprintf(fp0, "%f\t", tau[i][j] / ms_to_step);
		}
		fprintf(fp0, "\n");
	}
	fclose(fp0);

	fp0 = fopen("E_syn_matrix.txt", "w+");
	for (int i = 0; i < Neuron_count; i++)
	{
		for (int j = 0; j < Neuron_count; j++)
		{
			fprintf(fp0, "%f\t", E_syn_matrix[i][j]);
		}
		fprintf(fp0, "\n");
	}
	fclose(fp0);

	/*
	// Write to file number of links for each neuron
	fp0 = fopen("links.txt", "w+");
	for (int i = 0; i < Node_count; i++)
	{
		int links_count = 0;
		for (int j = 0; j < Node_count; j++)
		{
			if (N_A[i][j] == 1)
			{
				links_count++;
			}
		}
		fprintf(fp0, "%d\n", (int)links_count);
	}
	fclose(fp0);*/

	// setlocale(LC_NUMERIC, "French_Canada.1252");

	fp0 = fopen("N_B.txt", "w+");
	for (int i = 0; i < Neuron_count; i++)
	{
		for (int j = 0; j < N_C[i]; j++)
		{
			fprintf(fp0, "%d\t", (int)N_B[i][j]);
		}
		fprintf(fp0, "\n");
	}
	fclose(fp0);

	fp0 = fopen("N_C.txt", "w+");
	for (int i = 0; i < Neuron_count; i++)
	{
		fprintf(fp0, "%d\n", (int)N_C[i]);
	}
	fclose(fp0);

	fp0 = fopen("Matrix_sizes.txt", "w+");
	fprintf(fp0, "%d\t%d\n", (int)Neuron_width, (int)Neuron_height);
	fclose(fp0);

	fp0 = fopen("E_syn_matrix.txt", "w+");
	for (int i = 0; i < Neuron_count; i++)
	{
		for (int j = 0; j < Neuron_count; j++)
		{
			fprintf(fp0, "%f\t", E_syn_matrix[i][j]);
		}
		fprintf(fp0, "\n");
	}
	fclose(fp0);

	// Initial values
	/*for (int i = 0; i < Equations_count; i++)
	{
	  f[i] = 0;
    }*/

	// I_app_min = 1.1;
	// I_app_max = 1.5;

	FILE *fp_I_app;
	fp_I_app = fopen("I_app.txt", "w");
	for (int i = 0; i < Neuron_count; i++)
	{
		//I_app[i] = 0.8; // RandomD(I_app_min, I_app_max);
		fprintf(fp_I_app, "%f\n", I_app[i]);
	}
	fclose(fp_I_app);

	double percent_stable_state = 0.50;
	double eps_persent = 0.20;

	double V0 = -58.7085;
	double m0 = 0.0953;
	double n0 = 0.000913;
	double h0 = 0.3662;

	double V1 = 14.8409;
	double m1 = 0.9174;
	double n1 = 0.0140;
	double h1 = 0.0539;

	/*for (int i = 0; i < Node_count; i++) // init only for neurons
	{
	  double random = RandomD(0, 1);

	  SetV(i, random < percent_stable_state ? V0 + RandomD(-V0 * eps_persent, V0 * eps_persent) : V1 + RandomD(-V1 * eps_persent, V1 * eps_persent)); // V
	  Setm(i, random < percent_stable_state ? m0 + RandomD(-m0 * eps_persent, m0 * eps_persent) : m1 + RandomD(-m1 * eps_persent, m1 * eps_persent)); // m
	  Setn(i, random < percent_stable_state ? n0 + RandomD(-n0 * eps_persent, n0 * eps_persent) : n1 + RandomD(-n1 * eps_persent, n1 * eps_persent)); // n
	  Seth(i, random < percent_stable_state ? h0 + RandomD(-h0 * eps_persent, h0 * eps_persent) : h1 + RandomD(-h1 * eps_persent, h1 * eps_persent)); // h
	}*/

	for (int i = 0; i < Neuron_count; i++) // init only for neurons
	{
		// double random = RandomD(0, 1);

		/*SetV(i, random < percent_stable_state ? V0 : V1); // V
		Setm(i, random < percent_stable_state ? m0 : m1); // m
		Setn(i, random < percent_stable_state ? n0 : n1); // n
		Seth(i, random < percent_stable_state ? h0 : h1); // h*/

		// SetV(i, RandomD(-80, 20)); // V
		// Setm(i, RandomD(0, 1)); // m
		// Setn(i, RandomD(0, 1)); // n
		// Seth(i, RandomD(0, 1)); // h

		/*SetV(i, V0); // V
		Setm(i, m0); // m
		Setn(i, n0); // n
		Seth(i, h0); // h*/

		/*SetV_P(i, random < percent_stable_state ? V0 : V1); // V
		Setm_P(i, random < percent_stable_state ? m0 : m1); // m
		Setn_P(i, random < percent_stable_state ? n0 : n1); // n
		Seth_P(i, random < percent_stable_state ? h0 : h1); // h*/

		double random = RandomD(0, 1);

		SetV(i, random < percent_stable_state ? V0 + RandomD(-V0 * eps_persent, V0 * eps_persent) : V1 + RandomD(-V1 * eps_persent, V1 * eps_persent)); // V
		Setm(i, random < percent_stable_state ? m0 + RandomD(-m0 * eps_persent, m0 * eps_persent) : m1 + RandomD(-m1 * eps_persent, m1 * eps_persent)); // m
		Setn(i, random < percent_stable_state ? n0 + RandomD(-n0 * eps_persent, n0 * eps_persent) : n1 + RandomD(-n1 * eps_persent, n1 * eps_persent)); // n
		Seth(i, random < percent_stable_state ? h0 + RandomD(-h0 * eps_persent, h0 * eps_persent) : h1 + RandomD(-h1 * eps_persent, h1 * eps_persent)); // h
	}

	// double E_syn0 = 0;	 // Excitatory neuron
	// double E_syn1 = -90; // Inhibitory neuron
	// double p_inhib = 0.0;

	//fp_Esyn = fopen("E_syn.txt", "w+");
	
	//for (int i = 0; i < Neuron_count; i++)
	//{
		/*
		double x = RandomD(0, 1);

		if (x < p_inhib)
			E_syn[i] = E_syn1;
		else
			E_syn[i] = E_syn0;*/

	//	fprintf(fp_Esyn, "%f\n", E_syn[i]);
	//}
	//fclose(fp_Esyn);

	fp0 = fopen("g_syn.txt", "w+");
	
	for (int i = 0; i < Neuron_count; i++)
	{
		//g_syn[i] = g_syn_param;

		fprintf(fp0, "%f\n", g_syn[i]);
	}
	fclose(fp0);

	GenerateMeander();

	double t = t_start;

	double start_rk4, end_rk4;
	clock_t omp_start_rk4, omp_end_rk4;

	// if (omp_enable)
	//	omp_start_rk4 = omp_get_wtime();
	// else
	start_rk4 = clock();

	int lastPercent = -1;

	FillVOldFromCurrent();

	k = new double *[Equations_count];
	for (int i = 0; i < Equations_count; i++)
		k[i] = new double[4];

	phi_k1 = new double[Equations_count];
	phi_k2 = new double[Equations_count];
	phi_k3 = new double[Equations_count];

	if (write_file)
	{
		fp_I_stim = fopen("results_I_stim.txt", "w+");
		fp_I_syn = fopen("results_I_syn.txt", "w+");
		fp_V = fopen("results_V.txt", "w+");
		fp_V_spikes = fopen("results_V_spikes.txt", "w+");
	}

	double *f_next = new double[Equations_count];

	while (t < t_max || Approximately(t, t_max))
	{
		if (write_file)
		{
			fprintf(fp_I_stim, "%f\t", t);
			fprintf(fp_V, "%f\t", t);
			fprintf(fp_V_spikes, "%f\t", t);
		}

		for (int i = 0; i < Neuron_count; i++)
		{
			/*
			if (t > last_meander_end[i])
			{
				GenerateRandomMeander(i, t);
				last_meander_end[i] = Meander_start_from_zero[i] + Duration;
			}*/

			if (write_file)
				fprintf(fp_I_stim, "%f\t", I_stim(i, t));
		}

		if (write_file)
		{
			fprintf(fp_I_stim, "\n");

			for (int i = 0; i < Neuron_count; i++)
			{
				if (isnan(V(i)))
					return 1;

				fprintf(fp_V, i == Neuron_count - 1 ? "%f" : "%f\t", V(i)); // V
			}

			// for (int i = 5; i < Equations_count; i += Equations_count)
			//	fprintf(fp_m, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // m

			// for (int i = 6; i < Equations_count; i += Equations_count)
			//	fprintf(fp_n, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // n

			// for (int i = 7; i < Equations_count; i += Equations_count)
			//	fprintf(fp_h, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // h

			fprintf(fp_V, "\n");
			// fprintf(fp_m, "\n");
			// fprintf(fp_n, "\n");
			// fprintf(fp_h, "\n");
		}

		RungeKutta(t, dt, f, f_next);
		// Euler(t, dt, f, f_next);

		#pragma omp parallel for
		for (int i = 0; i < Neuron_count; i++)
		{
			int index = Equations_per_neuron * i;
			double diff = f_next[index] - f[index];

			fprintf(fp_V_spikes,
					i == Neuron_equations_count - 1 ? "%d" : "%d\t",
					diff < 0 && f_diff[index] > 0 && f[index] > -10 && (V_spikes[i].size() == 0 || t - V_spikes[i].back() > 0.001) ? 1 : 0);

			if (diff < 0 && f_diff[index] > 0 && f[index] > -10 && (V_spikes[i].size() == 0 || t - V_spikes[i].back() > 0.001))
			{
				V_spikes[i].push_back(t);
			}

			f_diff[index] = diff;
		}

		fprintf(fp_V_spikes, "\n");

		CopyArray(f_next, f, Equations_count);

		UpdateVOld();

		fprintf(fp_I_syn, "%f\t", t);

		for (int i = 0; i < Neuron_count; i++)
			fprintf(fp_I_syn, i == Neuron_count - 1 ? "%f" : "%f\t", I_syn_buffer[i]);

		fprintf(fp_I_syn, "\n");

		t += dt;

		int percent = (int)(100 * (t - t_start) / (t_max - t_start));
		if (percent != lastPercent)
		{
			printf("Progress: %d%%\n", percent);
			lastPercent = percent;
		}
	}

	bool anySpikes = false;

	for (int i = 0; i < Neuron_count; i++)
	{
		if (V_spikes[i].size() > 0)
		{
			anySpikes = true;
			break;
		}
	}

	if (anySpikes)
	{
		// FREQS
		{
			double *V_mean_freqs = new double[Neuron_count];
			double *V_STD_freqs = new double[Neuron_count];

			// #pragma omp parallel for
			for (int i = Neuron_count - 1; i < Neuron_count; i++)
			{
				list<double>::iterator it_V_spikes = V_spikes[i].begin();

				while (it_V_spikes != V_spikes[i].end() && *it_V_spikes < 0.25 * t_max)
				{
					V_spikes[i].pop_front();
					it_V_spikes = V_spikes[i].begin();
				}

				list<double> V_freqs;

				it_V_spikes = V_spikes[i].begin();

				V_mean_freqs[i] = 0;
				V_STD_freqs[i] = 0;

				if (V_spikes[i].size() <= 1)
				{
					continue;
				}
				else
				{
					for (int j = 1; j < V_spikes[i].size(); j++)
					{
						double first = *it_V_spikes;
						advance(it_V_spikes, 1);
						double next = *it_V_spikes;

						double T = next - first;
						V_freqs.push_back(1 / T);
					}
				}

				list<double>::iterator it_Freq = V_freqs.begin();

				for (int j = 0; j < V_freqs.size(); j++)
				{
					V_mean_freqs[i] += *it_Freq;
					advance(it_Freq, 1);
				}

				V_mean_freqs[i] /= V_freqs.size();

				it_Freq = V_freqs.begin();

				for (int j = 0; j < V_freqs.size(); j++)
				{
					V_STD_freqs[i] += pow(*it_Freq - V_mean_freqs[i], 2);
					advance(it_Freq, 1);
				}

				V_STD_freqs[i] /= V_freqs.size();
				V_STD_freqs[i] = sqrt(V_STD_freqs[i]);
			}

			fp0 = fopen("V_mean_freqs.txt", "w+");
			for (int i = Neuron_count - 1; i < Neuron_count; i++)
			{
				fprintf(fp0, "%f\t", V_mean_freqs[i]);
			}
			fclose(fp0);

			fp0 = fopen("V_STD_freqs.txt", "w+");
			for (int i = Neuron_count - 1; i < Neuron_count; i++)
			{
				fprintf(fp0, "%f\t", V_STD_freqs[i]);
			}
			fclose(fp0);

			delete[] V_mean_freqs;
			delete[] V_STD_freqs;
		}

		// SCF
		/*
		{
			fp0 = fopen("SCF.txt", "w+");

			double tau_scf_start = 0;
			double tau_scf_end = 0.500;
			double tau_scf_step = 0.0001;

			int first = 0;

			for (double tau_scf = tau_scf_start; tau_scf <= tau_scf_end; tau_scf += tau_scf_step)
			{
				fprintf(fp0, "%f\t", tau_scf);

				for (int second = 0; second < Neuron_count; second++)
				{
					double delta_scf = 0.001;
					int match_count = 0;

					list<double>::iterator second_spikes = V_spikes[second].begin();

					double correlation_scf = 0;

					if (V_spikes[second].size() > 0)
					{
						for (int i = 0; i < V_spikes[second].size(); i++)
						{
							double second_spike_time = *second_spikes;
							double target_first_spike_time = second_spike_time - tau_scf;

							list<double>::iterator first_spikes = V_spikes[first].begin();
							bool match = false;

							for (int j = 0; j < V_spikes[first].size(); j++)
							{
								double first_spike = *first_spikes;

								if (first_spike >= target_first_spike_time - delta_scf && first_spike <= target_first_spike_time + delta_scf)
								{
									match = true;
									break;
								}

								advance(first_spikes, 1);
							}

							if (match)
								match_count++;

							advance(second_spikes, 1);
						}

						correlation_scf = (double)match_count / V_spikes[second].size();
					}

					fprintf(fp0, "%f\t", correlation_scf);
				}

				fprintf(fp0, "\n");
			}

			fclose(fp0);
		}
		*/
	}
	else
		printf("Spikes are empty");

	// delete[] f_next;

	if (write_file)
	{
		// fclose(fp_res);
		fclose(fp_I_stim);
		fclose(fp_I_syn);
		fclose(fp_V);
		// fclose(fp_m);
		// fclose(fp_n);
		// fclose(fp_h);
		fclose(fp_V_spikes);
	}

	double extime_rk4 = 0;

	if (omp_enable)
	{
		omp_end_rk4 = omp_get_wtime();
		extime_rk4 = (double)(omp_end_rk4 - omp_start_rk4);// / CLOCKS_PER_SEC;
	}
	else
	{
		end_rk4 = clock();
		extime_rk4 = (double)(end_rk4 - start_rk4); // / CLOCKS_PER_SEC;
	}

	int minutes = (int)extime_rk4 / 60;
	int seconds = (int)extime_rk4 % 60;
	printf("\nExecution time was: %d minutes %d seconds\n ", minutes, seconds);

	/*int nth;
	#pragma omp parallel
	{
		#pragma omp master
		nth = omp_get_num_threads();
	}*/

	fp0 = fopen("time_exec.txt", "a");
	// fprintf(fp0, "%d %lf\n", nth, extime_rk4);
	fprintf(fp0, "%lf\n", extime_rk4);
	fclose(fp0);

	delete[] f;

	for (int i = 0; i < Neuron_count; i++)
		delete[] N_A[i];

	delete[] N_A;

	for (int i = 0; i < Neuron_count; i++)
		delete[] N_B[i];

	delete[] N_B;

	delete[] N_C;

	delete[] E_syn;
	delete[] I_app;
	delete[] V_spikes;
	// delete[] V_spikes_Freq;

	delete[] Meander_start_from_zero;
	delete[] Meander_width;
	delete[] Meander_height;
	delete[] Meander_interval;

	for (int i = 0; i < Neuron_count; i++)
		delete[] E_syn_matrix[i];

	delete[] E_syn_matrix;

	for (int i = 0; i < Neuron_count; i++)
		delete[] tau[i];

	delete[] tau;

	for (int i = 0; i < Neuron_count; i++)
		delete[] V_old_array[i];

	delete[] V_old_array;

	for (int i = 0; i < Equations_count; i++)
		delete[] k[i];

	delete[] k;

	delete[] phi_k1;
	delete[] phi_k2;
	delete[] phi_k3;
}