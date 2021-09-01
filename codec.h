#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include<string.h>
#include "arithmetic_codec.h"

#define N 256   // buffer_size
#define N_C 8	// maximum condition number 

const char path_source_file[] = "..\\Markov_source_09.txt";
const char path_side_information_file[] = "..\\Markov_source_1_si_1.txt";
const char path_code_word[] = "..\\1";
const char path_output_file[] = "..\\deoutput.txt";

int Length = 0;

char* input;		//source
char* SI;			//side information 
char* output;		//output

double Probability[N_C];
int Count_0[N_C];
int Count_s[N_C];

int condition_num = 1; //condition number	
int prob_c = 0;
int max_length = pow(2.0, condition_num) - 1;

const int p_n = 2;	//purging rate = 1/p_n;

double p_c = 0.01;	//BSC cross rate

double py0 = 0;		//P(Y=0)
double py1 = 0;		//P(Y=1)

int _pos = 1;

double mu = 0.01;		//forbidden length
typedef struct BiTNode
{
	Arithmetic_Codec   Codec;
	Static_Data_Model   static_model;
	int in_num;
	char X;
	long double lamda;
	long double Lambda;
	char* output;
	int prob_c;
	double Probability[N_C];
	int Count_0[N_C];
	int Count_s[N_C];
}BiTNode, * BiTree;

BiTree buffer[N];		//decoding tree buffer

unsigned Encoder();
unsigned SPDCAC_encode(Static_Data_Model& model, Arithmetic_Codec& encoder);
unsigned SPDCAC_encode_2(Static_Data_Model& model, Arithmetic_Codec& encoder);	// without forbidden symbol
void statistics();
void Decoder();
BiTree SPDCAC_decode(BiTree T, int n);
BiTree SPDCAC_decode_2(BiTree T, int n);
BiTree SPDCAC_assume(BiTree T, int n);
void posterior_probaility_decode(BiTree T, char X, char Y);
void DeleteNode();
void compare();
void posterior_probaility_assume(BiTree T, char X, char Y, int X_n_);
double adaptive_prob(BiTree T);
void sort(int n);
int quicksort(int len);
int quicksort_r(int start, int end);