
#include<time.h>
#include<malloc.h>
#include<math.h>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include "arithmetic_codec.h"

#define N 256   // buffer_size
#define N_C 4	// maximum conditional number 
#define Length 1024		//source length

const char path_source_file[] = "..\\Markov_source.txt";
const char path_side_information_file[] = "..\\Markov_source_si.txt";
const char path_code_word[] = "..\\code_word";
const char path_output_file[] = "..\\output.txt";

char* input;		//source
char* SI;			//side information 
char* output;		//output

double Probability_0[N_C];	//conditional probability
int Count_0[N_C];
int Count_s[N_C];

int condition_num = 1; //conditional number of context model	
int prob_c = 0;		//index of context
int max_length = pow(2.0, condition_num) - 1;

const int p_n = 2;		//purging rate = 1/p_n;

double p_c = 0.0137;		//BSC cross probability

int _pos = 1;

double mu = 0.01;		//forbidden length
typedef struct BiTNode
{
	Arithmetic_Codec   Codec;
	Static_Data_Model   static_model;
	int in_num;
	char X;
	long double lamda;	//posterior probability 
	long double Lambda;	//accumulative posterior probability
	char output[Length];	//decoded results of this branch
	int prob_c;
	double Probability_0[N_C];
	int Count_0[N_C];
	int Count_s[N_C];
}BiTNode, * BiTree;

BiTree buffer[N];		//decoding tree buffer

unsigned Encoder();
unsigned SPDCAC_encoder(Static_Data_Model& model, Arithmetic_Codec& encoder);
unsigned SPDCAC_encoder_2(Static_Data_Model& model, Arithmetic_Codec& encoder);	// without forbidden symbol
void Decoder();
BiTree SPDCAC_decoder(BiTree T, int n);
BiTree SPDCAC_decoder_2(BiTree T, int n);	// without forbidden symbol
BiTree SPDCAC_assume(BiTree T, int n);
void posterior_probaility_decode(BiTree T, char X, char Y);
void DeleteNode();
void compare();
void posterior_probaility_assume(BiTree T, char X, char Y, int X_n_);
double update_prob(BiTree T);
int quicksort(int len);
int quicksort_r(int start, int end);

int main()
{
	clock_t start, finish;
	start = clock();

	FILE* fi = fopen(path_source_file, "r");
	input = new char[Length];
	fread(input, 1, Length, fi);
	fclose(fi);

	SI = new char[Length];
	FILE* fs = fopen(path_side_information_file, "r");
	fread(SI, 1, Length, fs);
	fclose(fs);

	//int ccc = 0;
	//for (int i = 0; i < Length; i++)
	//{
	//	if (input[i] != SI[i])
	//		ccc++;
	//}
	//p_c = (double)ccc / Length;

	int code_length = Encoder();
	printf("code_word_length= %d bytes\n", code_length);
	printf("code rate= %f bit/symbol\n", (double)code_length / (Length / 8));	//in the *.txt file, each binary symbol is a bit
	Decoder();

	finish = clock();
	printf("running time= %fs\n", (double)(finish - start) / CLOCKS_PER_SEC);

	system("pause");
	return 0;
}

unsigned Encoder()
{
	Arithmetic_Codec   Codec;
	Static_Data_Model   static_model;
	int code_word_length = 0;

	for (int i = 0; i < N_C; i++)
	{
		Probability_0[i] = 0.5;
		Count_0[i] = 1;
		Count_s[i] = 2;
	}

	FILE* fid = fopen(path_code_word, "wb");

	Codec.set_buffer(Length, 0);
	SPDCAC_encoder(static_model, Codec);
	//SPDCAC_encoder_2(static_model, Codec);	// without forbidden symbol
	code_word_length = Codec.write_to_file(fid);

	fclose(fid);
	return code_word_length;
}

unsigned SPDCAC_encoder(Static_Data_Model& model, Arithmetic_Codec& encoder)
{
	encoder.start_encoder();
	prob_c = 0;
	unsigned bit = 0;

	for (int i = 0; i < Length; i++)
	{
		prob_c &= max_length;

		if ((i + 1) % p_n != 0)
		{
			double probability[3];

			if (Probability_0[prob_c] > (1 - Probability_0[prob_c]))
			{
				probability[0] = Probability_0[prob_c] * mu;
				probability[1] = Probability_0[prob_c] * (1 - mu);
				probability[2] = 1 - Probability_0[prob_c];
			}
			else
			{
				probability[0] = Probability_0[prob_c];
				probability[1] = (1 - Probability_0[prob_c]) * (1 - mu);
				probability[2] = (1 - Probability_0[prob_c]) * mu;
			}
			for (int k = 0; k < 3; k++)
			{
				if (probability[k] < 0.0001)
					probability[k] = 0.0001;
				else if (probability[k] > 0.9999)
					probability[k] = 0.9999;
			}

			model.set_distribution(3, probability);

			if (Probability_0[prob_c] > (1 - Probability_0[prob_c]))
			{
				switch (input[i] - 48)
				{
				case 0:
					bit = 1;
					break;
				case 1:
					bit = 2;
					break;
				}
				encoder.encode(bit, model);
			}
			else
				encoder.encode(unsigned(input[i] - 48), model);
		}

		prob_c &= max_length;
		if ((input[i] - 48) == 0)
			Count_0[prob_c]++;
		Count_s[prob_c]++;
		Probability_0[prob_c] = (double)Count_0[prob_c] / Count_s[prob_c];

		prob_c = (prob_c << 1) | (input[i] - 48);
	}
	return 0;
}

unsigned SPDCAC_encoder_2(Static_Data_Model& model, Arithmetic_Codec& encoder)	// without forbidden symbol
{
	encoder.start_encoder();
	prob_c = 0;

	for (int i = 0; i < Length; i++)
	{
		prob_c &= max_length;

		if ((i + 1) % p_n != 0)
		{
			double probability[2];
			probability[0] = Probability_0[prob_c];
			probability[1] = 1 - Probability_0[prob_c];
			model.set_distribution(2, probability);
			encoder.encode(unsigned(input[i] - 48), model);
		}

		prob_c &= max_length;
		if ((input[i] - 48) == 0)
			Count_0[prob_c]++;
		Count_s[prob_c]++;
		Probability_0[prob_c] = (double)Count_0[prob_c] / Count_s[prob_c];

		prob_c = (prob_c << 1) | (input[i] - 48);
	}
	return 0;
}

void Decoder()
{
	buffer[0] = (BiTree)malloc(sizeof(BiTNode));

	for (int i = 0; i < N_C; i++)
	{
		buffer[0]->Probability_0[i] = 0.5;
		buffer[0]->Count_0[i] = 1;
		buffer[0]->Count_s[i] = 2;
	}

	buffer[0]->Codec.set_buffer(Length, 0);
	FILE* fp = fopen(path_code_word, "rb");
	buffer[0]->Codec.read_from_file(fp);
	fclose(fp);

	buffer[0]->in_num = 0;
	buffer[0]->Lambda = 0;
	buffer[0]->lamda = 0;
	buffer[0]->prob_c = 0;

	while (1)
	{
		compare();
		if (buffer[0]->in_num >= Length)
		{
			break;
		}
		else
		{
			int pos_i = _pos;
			for (int i = 0; i < pos_i; i++)
			{
				if (buffer[i] != NULL)
				{
					buffer[i]->prob_c &= max_length;
					if ((buffer[i]->in_num + 1) % p_n != 0)
						SPDCAC_decoder(buffer[i], i);
						//SPDCAC_decoder_2(buffer[i], i);	// without forbidden symbol
					else
						SPDCAC_assume(buffer[i], i);
				}
			}
		}
	}

	int BER = 0;
	for (int i = 0; i < Length; i++)
	{
		if (buffer[0]->output[i] != input[i])
			BER++;
	}
	printf("BER= %f\n", (double)BER / Length);

	FILE* fo = fopen(path_output_file, "w");
	fwrite(buffer[0]->output, 1, Length, fo);
	fclose(fo);

	for (int i = 0; i < N; i++)
	{
		if (buffer[i])
		{
			free(buffer[i]);
		}
	}
}

BiTree SPDCAC_decoder(BiTree T, int n)
{
	double probability[3];
	if (T->Probability_0[T->prob_c] > (1 - T->Probability_0[T->prob_c]))
	{
		probability[0] = T->Probability_0[T->prob_c] * mu;
		probability[1] = T->Probability_0[T->prob_c] * (1 - mu);
		probability[2] = 1 - T->Probability_0[T->prob_c];
	}
	else
	{
		probability[0] = T->Probability_0[T->prob_c];
		probability[1] = (1 - T->Probability_0[T->prob_c]) * (1 - mu);
		probability[2] = (1 - T->Probability_0[T->prob_c]) * mu;
	}
	for (int i = 0; i < 3; i++)
	{
		if (probability[i] < 0.0001)
			probability[i] = 0.0001;
		else if (probability[i] > 0.9999)
			probability[i] = 0.9999;
	}
	T->static_model.set_distribution(3, probability);
	T->X = T->Codec.decode(T->static_model) + 48;
	if (T->Probability_0[T->prob_c] > (1 - T->Probability_0[T->prob_c]))
	{
		switch (T->X)
		{
		case '1':
			T->X = '0';
			break;
		case '2':
			T->X = '1';
			break;
		default:
			T->X = '2';
			break;
		}
	}
	if (T->X == '2')
	{
		free(buffer[n]);
		buffer[n] = NULL;
	}
	else
	{
		posterior_probaility_decode(T, T->X, SI[T->in_num]);
		update_prob(T);
		T->prob_c = (T->prob_c << 1) | (T->X - 48);
		T->prob_c &= max_length;
		T->output[T->in_num] = T->X;
		T->in_num++;
		T->Lambda = T->Lambda + T->lamda;
	}
	return T;
}

BiTree SPDCAC_decoder_2(BiTree T, int n)	// without forbidden symbol
{
	double probability[2];
	probability[0] = T->Probability_0[T->prob_c];
	probability[1] = 1 - T->Probability_0[T->prob_c];
	T->static_model.set_distribution(2, probability);
	T->X = T->Codec.decode(T->static_model) + 48;
	posterior_probaility_decode(T, T->X, SI[T->in_num]);
	update_prob(T);
	T->output[T->in_num] = T->X;
	T->in_num++;
	T->Lambda = T->Lambda + T->lamda;
	T->prob_c = (T->prob_c << 1) | (T->X - 48);
	T->prob_c &= max_length;
	return T;
}

BiTree SPDCAC_assume(BiTree T, int n)
{
	char che = SI[T->in_num];

	buffer[_pos] = (BiTree)malloc(sizeof(BiTNode));
	memcpy(buffer[_pos], buffer[n], sizeof(BiTNode));

	//the branch of symbol 0
	buffer[n]->X = '0';
	buffer[n]->output[buffer[n]->in_num] = buffer[n]->X;
	posterior_probaility_assume(buffer[n], '0', che, buffer[n]->prob_c);
	buffer[n]->Lambda = buffer[n]->Lambda + buffer[n]->lamda;
	buffer[n]->in_num = buffer[n]->in_num + 1;
	update_prob(buffer[n]);
	buffer[n]->prob_c = (buffer[n]->prob_c << 1) | 0;
	buffer[n]->prob_c &= max_length;

	//the branch of symbol 1
	buffer[_pos]->X = '1';
	buffer[_pos]->output[buffer[_pos]->in_num] = buffer[_pos]->X;
	posterior_probaility_assume(buffer[_pos], '1', che, buffer[_pos]->prob_c);
	buffer[_pos]->Lambda = buffer[_pos]->Lambda + buffer[_pos]->lamda;
	buffer[_pos]->in_num = buffer[_pos]->in_num + 1;
	update_prob(buffer[_pos]);
	buffer[_pos]->prob_c = (buffer[_pos]->prob_c << 1) | 1;
	buffer[_pos]->prob_c &= max_length;

	_pos++;
	return T;
}

void posterior_probaility_decode(BiTree T, char X, char Y)
{
	double _p = 1 - p_c;
	long double pxy = 0;
	if (X == Y)
	{
		pxy = _p;
	}
	else
	{
		pxy = p_c;
	}
	if (pxy == 0) 		
		pxy = 0.000000001;
	T->lamda = log(pxy) / log(2.0);
}

void posterior_probaility_assume(BiTree T, char X, char Y, int X_n_)
{
	double _p = 1 - p_c;
	long double pxy = 0;
	if (X == Y)
	{
		if (X == '0')
			pxy = _p * T->Probability_0[X_n_];
		else if (X == '1')
			pxy = _p * (1 - T->Probability_0[X_n_]);
	}
	else
	{
		if (X == '0')
			pxy = p_c * T->Probability_0[X_n_];
		else if (X == '1')
			pxy = p_c * (1 - T->Probability_0[X_n_]);
	}
	if (pxy == 0) 		
		pxy = 0.000000001;
	T->lamda = log(pxy) / log(2.0);
}

void compare()
{
	int i = 0, j;
	while (i <= N - 1)
	{
		if (buffer[i] == NULL)
		{
			for (j = i + 1; j < N; j++)
			{
				if (buffer[j] != NULL)
				{
					buffer[i] = buffer[j];
					buffer[j] = NULL;
					break;
				}
			}
		}
		i++;
	}

	int nu = N;
	for (int i = 0; i < N; i++)
	{
		if (buffer[i] == NULL)
		{
			nu = i;
			break;
		}
	}
	quicksort(nu);

	DeleteNode();

	if(_pos>1)
		for (int i = 2; i < N; i++)
		{
			if (buffer[i] == NULL)
			{
				_pos = i;
				break;
			}
		}
}

void DeleteNode()
{
	for (int i = N / 2; i <= N - 1; i++)
	{
		if (buffer[i])
		{
			free(buffer[i]);
			buffer[i] = NULL;
		}
	}
}

int quicksort(int len)
{
	quicksort_r(0, len - 1);
	return 0;
}
int quicksort_r(int start, int end)
{
	if (start >= end) {
		return 0;
	}
	BiTree pivot;
	pivot = buffer[end];
	BiTree swp;
	//set a pointer to divide array into two parts
	//one part is smaller than pivot and another larger
	int pointer = start;
	int i;
	for (i = start; i < end; i++) {
		if (buffer[i]->Lambda > pivot->Lambda) {
			if (pointer != i) {
				//swap a[i] with a[pointer]
				//a[pointer] behind larger than pivot
				swp = buffer[i];
				buffer[i] = buffer[pointer];
				buffer[pointer] = swp;
			}
			pointer++;
		}
	}
	//swap back pivot to proper position
	swp = buffer[end];
	buffer[end] = buffer[pointer];
	buffer[pointer] = swp;

	quicksort_r(start, pointer - 1);
	quicksort_r(pointer + 1, end);
	return 0;
}

double update_prob(BiTree T)
{
	if (T->X == '0')
		T->Count_0[T->prob_c]++;
	T->Count_s[T->prob_c]++;
	T->Probability_0[T->prob_c] = (double)T->Count_0[T->prob_c] / T->Count_s[T->prob_c];

	return T->Probability_0[T->prob_c];
}


