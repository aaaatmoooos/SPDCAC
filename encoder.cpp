#include "codec.h"

unsigned Encoder()
{
	Arithmetic_Codec   Codec;
	Static_Data_Model   static_model;
	int code_word_length = 0;

	for (int i = 0; i < N_C; i++)
	{
		Probability[i] = 0;
		Count_0[i] = 1;
		Count_s[i] = 2;
	}

	FILE* fid = fopen(path_code_word, "wb");

	Codec.set_buffer(Length, 0);
	SPDCAC_encode(static_model, Codec);
	//SPDCAC_encode_2(static_model, Codec);	// without forbidden symbol
	code_word_length = Codec.write_to_file(fid);

	fclose(fid);
	return code_word_length;
}

unsigned SPDCAC_encode(Static_Data_Model& model, Arithmetic_Codec& encoder)
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

			if (Probability[prob_c] > (1 - Probability[prob_c]))
			{
				probability[0] = Probability[prob_c] * mu;
				probability[1] = Probability[prob_c] * (1 - mu);
				probability[2] = 1 - Probability[prob_c];
			}
			else
			{
				probability[0] = Probability[prob_c];
				probability[1] = (1 - Probability[prob_c]) * (1 - mu);
				probability[2] = (1 - Probability[prob_c]) * mu;
			}
			for (int k = 0; k < 3; k++)
			{
				if (probability[k] < 0.0001)
					probability[k] = 0.0001;
				else if (probability[k] > 0.9999)
					probability[k] = 0.9999;
			}

			model.set_distribution(3, probability);

			if (Probability[prob_c] > (1 - Probability[prob_c]))
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
		Probability[prob_c] = (double)Count_0[prob_c] / Count_s[prob_c];

		prob_c = (prob_c << 1) | (input[i] - 48);
	}
	return 0;
}

unsigned SPDCAC_encode_2(Static_Data_Model& model, Arithmetic_Codec& encoder)
{
	encoder.start_encoder();
	prob_c = 0;

	for (int i = 0; i < Length; i++)
	{
		prob_c &= max_length;

		if ((i + 1) % p_n != 0)
		{
			double probability[2];
			probability[0] = Probability[prob_c];
			probability[1] = 1 - Probability[prob_c];
			model.set_distribution(2, probability);
			encoder.encode(unsigned(input[i] - 48), model);
		}

		prob_c &= max_length;
		if ((input[i] - 48) == 0)
			Count_0[prob_c]++;
		Count_s[prob_c]++;
		Probability[prob_c] = (double)Count_0[prob_c] / Count_s[prob_c];

		prob_c = (prob_c << 1) | (input[i] - 48);
	}
	return 0;
}
