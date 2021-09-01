
#include <stdlib.h>
#include<string.h>
#include "arithmetic_codec.h"

const unsigned AC__MinLength = 0x01000000U;   // threshold for renormalization
const unsigned AC__MaxLength = 0xFFFFFFFFU;      // maximum AC interval length

                                           // Maximum values for binary models
const unsigned BM__LengthShift = 13;     // length bits discarded before mult.
const unsigned BM__MaxCount    = 1 << BM__LengthShift;  // for adaptive models

                                          // Maximum values for general models
const unsigned DM__LengthShift = 15;     // length bits discarded before mult.
const unsigned DM__MaxCount    = 1 << DM__LengthShift;  // for adaptive models

static void AC_Error(const char * msg)
{
  fprintf(stderr, "\n\n -> Arithmetic coding error: ");
  fputs(msg, stderr);
  fputs("\n Execution terminated!\n", stderr);
  exit(1);
}


inline void Arithmetic_Codec::propagate_carry(void)
{
  unsigned char * p;            // carry propagation on compressed data buffer
  for (p = ac_pointer - 1; *p == 0xFFU; p--) *p = 0;
  ++*p;
}

inline void Arithmetic_Codec::renorm_enc_interval(void)
{
  do {                                          // output and discard top byte
    *ac_pointer++ = (unsigned char)(base >> 24);
    base <<= 8;
  } while ((length <<= 8) < AC__MinLength);        // length multiplied by 256
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline void Arithmetic_Codec::renorm_dec_interval(void)
{
  do {                                          // read least-significant byte
    value = (value << 8) | unsigned(*++ac_pointer);
  } while ((length <<= 8) < AC__MinLength);        // length multiplied by 256
}

void Arithmetic_Codec::encode(unsigned data,
                              Static_Data_Model & M)
{
#ifdef _DEBUG
  if (mode != 1) AC_Error("encoder not initialized");
  if (data >= M.data_symbols) AC_Error("invalid data symbol");
#endif

  unsigned x, init_base = base;
                                                           // compute products
  if (data == M.last_symbol) {
    x = M.distribution[data] * (length >> DM__LengthShift);
    base   += x;                                            // update interval
    length -= x;                                          // no product needed
  }
  else {
    x = M.distribution[data] * (length >>= DM__LengthShift);
    base   += x;                                            // update interval
    length  = M.distribution[data+1] * length - x;
  }
             
  if (init_base > base) propagate_carry();                 // overflow = carry

  if (length < AC__MinLength) renorm_enc_interval();        // renormalization
}

unsigned Arithmetic_Codec::decode(Static_Data_Model & M)
{
#ifdef _DEBUG
  if (mode != 2) AC_Error("decoder not initialized");
#endif

  unsigned n, s, x, y = length;

  if (M.decoder_table) {              // use table look-up for faster decoding

    unsigned dv = value / (length >>= DM__LengthShift);
    unsigned t = dv >> M.table_shift;

    s = M.decoder_table[t];         // initial decision based on table look-up
    n = M.decoder_table[t+1] + 1;

    while (n > s + 1) {                        // finish with bisection search
      unsigned m = (s + n) >> 1;
      if (M.distribution[m] > dv) n = m; else s = m;
    }
                                                           // compute products
    x = M.distribution[s] * length;
    if (s != M.last_symbol) y = M.distribution[s+1] * length;
  }

  else {                                  // decode using only multiplications

    x = s = 0;
    length >>= DM__LengthShift;
    unsigned m = (n = M.data_symbols) >> 1;
                                                // decode via bisection search
    do {
      unsigned z = length * M.distribution[m];
      if (z > value) {
        n = m;
        y = z;                                             // value is smaller
      }
      else {
        s = m;
        x = z;                                     // value is larger or equal
      }
    } while ((m = (s + n) >> 1) != s);
  }

  value -= x;                                               // update interval
  length = y - x;

  if (length < AC__MinLength) renorm_dec_interval();        // renormalization

  return s;
}

Arithmetic_Codec::Arithmetic_Codec(void)
{
  mode = buffer_size = 0;
  new_buffer = code_buffer = 0;
}

Arithmetic_Codec::Arithmetic_Codec(unsigned max_code_bytes,
                                   unsigned char * user_buffer)
{
  mode = buffer_size = 0;
  new_buffer = code_buffer = 0;
  set_buffer(max_code_bytes, user_buffer);
}

Arithmetic_Codec::~Arithmetic_Codec(void)
{
  delete [] new_buffer;
}


void Arithmetic_Codec::set_buffer(unsigned max_code_bytes,
                                  unsigned char * user_buffer)
{
                                                  // test for reasonable sizes
  if ((max_code_bytes < 16) || (max_code_bytes > 0x1000000U))
    AC_Error("invalid codec buffer size");
  //if (mode != 0) AC_Error("cannot set buffer while encoding or decoding");

  if (user_buffer != 0) {                       // user provides memory buffer
    buffer_size = max_code_bytes;
    code_buffer = user_buffer;               // set buffer for compressed data
    delete [] new_buffer;                 // free anything previously assigned
    new_buffer = 0;
    return;
  }

  buffer_size = 0;
  new_buffer = 0;
  if (max_code_bytes <= buffer_size) return;               // enough available

  buffer_size = max_code_bytes;                           // assign new memory
  delete [] new_buffer;                   // free anything previously assigned
  if ((new_buffer = new unsigned char[buffer_size+16]) == 0) // 16 extra bytes
    AC_Error("cannot assign memory for compressed data buffer");
  code_buffer = new_buffer;                  // set buffer for compressed data
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void Arithmetic_Codec::start_encoder(void)
{
  /*if (mode != 0) AC_Error("cannot start encoder");
  if (buffer_size == 0) AC_Error("no code buffer set");*/

  mode   = 1;
  base   = 0;            // initialize encoder variables: interval and pointer
  length = AC__MaxLength;
  ac_pointer = code_buffer;                       // pointer to next data byte
}


void Arithmetic_Codec::start_decoder(void)
{
  //if (mode != 0) AC_Error("cannot start decoder");
  if (buffer_size == 0) AC_Error("no code buffer set");

                  // initialize decoder: interval, pointer, initial code value
  mode   = 2;
  length = AC__MaxLength;
  ac_pointer = code_buffer + 3;
  value = (unsigned(code_buffer[0]) << 24)|(unsigned(code_buffer[1]) << 16) |
          (unsigned(code_buffer[2]) <<  8)| unsigned(code_buffer[3]);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int Arithmetic_Codec::read_from_file(FILE * code_file)
{
  unsigned shift = 0, code_bytes = 0;
  int file_byte;
                      // read variable-length header with number of code bytes
  do {
    if ((file_byte = getc(code_file)) == EOF)
      AC_Error("cannot read code from file");
    code_bytes |= unsigned(file_byte & 0x7F) << shift;
    shift += 7;
  } while (file_byte & 0x80);
                                                       // read compressed data
  if (code_bytes > buffer_size) AC_Error("code buffer overflow");
  if (fread(code_buffer, 1, code_bytes, code_file) != code_bytes)
    AC_Error("cannot read code from file");

  start_decoder();                                       // initialize decoder
  return code_bytes;
}


unsigned Arithmetic_Codec::stop_encoder(void)
{
  if (mode != 1) AC_Error("invalid to stop encoder");
  mode = 0;

  unsigned init_base = base;            // done encoding: set final data bytes

  if (length > 2 * AC__MinLength) {
    base  += AC__MinLength;                                     // base offset
    length = AC__MinLength >> 1;             // set new length for 1 more byte
  }
  else {
    base  += AC__MinLength >> 1;                                // base offset
    length = AC__MinLength >> 9;            // set new length for 2 more bytes
  }

  if (init_base > base) propagate_carry();                 // overflow = carry

  renorm_enc_interval();                // renormalization = output last bytes

  unsigned code_bytes = unsigned(ac_pointer - code_buffer);
  if (code_bytes > buffer_size) AC_Error("code buffer overflow");

  return code_bytes;                                   // number of bytes used
}


unsigned Arithmetic_Codec::write_to_file(FILE * code_file)
{
  unsigned header_bytes = 0, code_bytes = stop_encoder(), nb = code_bytes;

                     // write variable-length header with number of code bytes
  do {
    int file_byte = int(nb & 0x7FU);
    if ((nb >>= 7) > 0) file_byte |= 0x80;
    if (putc(file_byte, code_file) == EOF)
      AC_Error("cannot write compressed data to file");
    header_bytes++;
  } while (nb);
                                                      // write compressed data
  if (fwrite(code_buffer, 1, code_bytes, code_file) != code_bytes)
    AC_Error("cannot write compressed data to file");

  return code_bytes + header_bytes;                              // bytes used
}


void Arithmetic_Codec::stop_decoder(void)
{
  if (mode != 2) AC_Error("invalid to stop decoder");
  mode = 0;
}


Static_Data_Model::Static_Data_Model(void)
{
  data_symbols = 0;
  distribution = 0;
}

Static_Data_Model::~Static_Data_Model(void)
{
  delete [] distribution;
}


void Static_Data_Model::set_distribution(unsigned number_of_symbols,
                                         const double probability[])
{
  if ((number_of_symbols < 2) || (number_of_symbols > (1 << 11)))
    AC_Error("invalid number of data symbols");

  if (data_symbols != number_of_symbols) {     // assign memory for data model
    data_symbols = number_of_symbols;
    last_symbol = data_symbols - 1;
	//memset(distribution, 0, 100);
	distribution = 0;
    delete [] distribution;
                                     // define size of table for fast decoding
    if (data_symbols > 16) {
      unsigned table_bits = 3;
      while (data_symbols > (1U << (table_bits + 2))) ++table_bits;
      table_size  = 1 << table_bits;
      table_shift = DM__LengthShift - table_bits;
      distribution = new unsigned[data_symbols+table_size+2];
      decoder_table = distribution + data_symbols;
    }
    else {                                  // small alphabet: no table needed
      decoder_table = 0;
      table_size = table_shift = 0;
      distribution = new unsigned[data_symbols];
    }
    if (distribution == 0) AC_Error("cannot assign model memory");
  }
                             // compute cumulative distribution, decoder table
  unsigned s = 0;
  double sum = 0.0, p = 1.0 / double(data_symbols);

  for (unsigned k = 0; k < data_symbols; k++) {
    if (probability) p = probability[k];
    if ((p < 0.0001) || (p > 0.9999)) AC_Error("invalid symbol probability");
    distribution[k] = unsigned(sum * (1 << DM__LengthShift));
    sum += p;
    if (table_size == 0) continue;
    unsigned w = distribution[k] >> table_shift;
    while (s < w) decoder_table[++s] = k - 1;
  }

  if (table_size != 0) {
    decoder_table[0] = 0;
    while (s <= table_size) decoder_table[++s] = data_symbols - 1;
  }

  if ((sum < 0.9999) || (sum > 1.0001)) AC_Error("invalid probabilities");
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
