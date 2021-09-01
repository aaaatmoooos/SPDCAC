
#ifndef ARITHMETIC_CODEC
#define ARITHMETIC_CODEC

#include <stdio.h>


class Static_Data_Model                       
{
public:

  Static_Data_Model(void);
 ~Static_Data_Model(void);

  unsigned model_symbols(void) { return data_symbols; }

  void set_distribution(unsigned number_of_symbols,
                        const double probability[] = 0);    

//private:  //  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
  unsigned * distribution, * decoder_table;
  unsigned data_symbols, last_symbol, table_size, table_shift;
  friend class Arithmetic_Codec;
};

class Arithmetic_Codec
{
public:

  Arithmetic_Codec(void);
 ~Arithmetic_Codec(void);
  Arithmetic_Codec(unsigned max_code_bytes,
                   unsigned char * user_buffer = 0);         

  unsigned char * buffer(void) { return code_buffer; }

  void set_buffer(unsigned max_code_bytes,
                  unsigned char * user_buffer = 0);          

  void     start_encoder(void);
  void     start_decoder(void);
  int     read_from_file(FILE * code_file);  

  unsigned stop_encoder(void);                 
  unsigned write_to_file(FILE * code_file);   
  void     stop_decoder(void);

  void     put_bit(unsigned bit);
  unsigned get_bit(void);

  void     put_bits(unsigned data, unsigned number_of_bits);
  unsigned get_bits(unsigned number_of_bits);

  void     encode(unsigned data,
                  Static_Data_Model &);
  unsigned decode(Static_Data_Model &);


//private:  //  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
  void propagate_carry(void);
  void renorm_enc_interval(void);
  void renorm_dec_interval(void);
  unsigned char * code_buffer, * new_buffer, * ac_pointer;
  unsigned base, value, length;                     
  unsigned buffer_size, mode;     
};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#endif
