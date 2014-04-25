/*
Title       : Experiment 4
Project     : Mtk 5th generation communication system
Author      : Ming Jie Yang
Date        : 
Description : 
              
              
              
              
              
*/
#include <itpp/itcomm.h>
#include <ctime>
#include <cstring>

using namespace itpp;

//These lines are needed for use of cout and endl
using std::cout;
using std::endl;
using std::string;

//Function defined here
double eval_avg_power(const cvec& symbol_vec);
bvec encode(Convolutional_Code& nsc, int constraint_length, const bvec& encoder_input, int blockSize, bool verbose);
bvec decode(Convolutional_Code& nsc, int constraint_length, const bvec& decoder_input, int blockSize, bool verbose);
void zero_pad_back(cvec& vec_in, int m_in);

int main( int argc, char* argv[])
{
  //Arg initial here
  int Number_of_bits = 2048;

  //Read arg Simple
  if( argc == 2) Number_of_bits = atoi( argv[1]);

  //Declarations of scalars and vectors:
  int i;
  double Ps, N0, dist_1, dist_2, h1, h2, Eb;
  double EbN0;
  double rcvd_power_1, rcvd_power_2;
  int nFFT, nCylicPrefix;

  vec alpha;                                            //vec is a vector containing double
  vec bit_error_rate_1, bit_error_rate_2;
  vec ber_theo_1, ber_theo_2;                           // Theoretical results for multiple access

  bvec transmitted_bits_1, transmitted_bits_2;          //bvec is a vector containing bits
  bvec received_bits_1, received_bits_2;

  cvec transmitted_symbols_1, transmitted_symbols_2, transmitted_symbols; //cvec is a vector containing double_complex
  cvec received_symbols_1, received_symbols_2, feedback_symbols_1, feedback_symbols_2;
  cvec ofdm_symbols_1, ofdm_symbols_2;

  // Declarations of classes:
  QPSK qpsk;                     //The QPSK modulator class
  AWGN_Channel awgn_channel;     //The AWGN channel class
  it_file ff;                    //For saving the results to file
  BERC berc;                     //Used to count the bit errors
  Real_Timer tt;                 //The timer used to measure the execution time
  OFDM ofdm;                     //OFDM modulator class

  MA_Filter<std::complex<double>, std::complex<double>, std::complex<double> > multipath_channel; //Simulate multi-path enviornment
  int ch_nb_taps = 2;
  cvec ch_imp_response( ch_nb_taps);
  cvec ini_state = to_cvec( ones(ch_nb_taps));
  ch_imp_response = to_cvec( randray( ch_nb_taps));
  ch_imp_response /= sqrt(sum_sqr( ch_imp_response));//normalized power profile
  cout << ch_imp_response << endl;
  multipath_channel.set_coeffs( ch_imp_response);
  multipath_channel.set_state(ini_state);//inital state is zero

  //Reset and start the timer:
  tt.tic();

  //Init:
  Ps = 5 * pow(10, -3);          //The transmitted energy per QPSK symbol is 1, 5e-3 W/MHz
  N0 = 4 * pow(10, -15);         //Thermal noise -144dBm/MHz
  dist_1 = 150;                  //Distance form transmitter to receiver 1 (meter)
  dist_2 = 220;                  //                   *                  2
  h1 = pow(dist_1, -4.5);        //Channel gain for receiver 1
  h2 = pow(dist_2, -4.5);        //
  cout << "h1 = " << h1 << " h2 = " << h2 << endl;
  alpha = linspace(0.5, 0.0, 21);//Simulate for different weight on power
  nFFT = 2048;                   //FFT size, default is 2048 LTE-a
  nCylicPrefix = 144;            //Length of Prefix, standard 144 (first prefix is different in real case)

  // test upconv
  cvec ocsi;
  cvec t_spanned;
  double unit_time = 3.25520833 * pow(10,-8); // 1/(2048*15k)
  const double pi = 3.1415926535897;
  const std::complex<double> unit_img (0.0,1.0);
  double fc = 2.6 * pow(10,9); // 2.6GHz
  t_spanned.set_size(2048, false);
  for (int i = 0; i < 2048; i++) {
    t_spanned[i] = exp(2*pi*fc*unit_time*i*unit_img);
  }

  // Declarations for equalizer & NSC coding
  ivec gen = "07 05";                                          //octal notation
  int constraint_length = 3;                                   //constraint length
  int blockSize = 13;//permutation length                //encoder output block size
  // other parameters
  int nb_bits_tail = blockSize / gen.length();                  //encoder input size + tail size
  int nb_bits = nb_bits_tail - (constraint_length - 1);        //encoder block size
  int nb_blocks;                                               //number of blocks
  // Convolutional code
  Convolutional_Code nsc;
  nsc.set_generator_polynomials(gen, constraint_length);

  //Allocate storage space for the result vector.
  //The "false" argument means "Do not copy the old content of the vector to the new storage area."
  bit_error_rate_1.set_size(alpha.length(), false);
  bit_error_rate_2.set_size(alpha.length(), false);

  //Storage for theoretical ber
  ber_theo_1.set_size(alpha.length(), false);
  ber_theo_2.set_size(alpha.length(), false);

  //Randomize the random number generators in it++:
  RNG_randomize();

  //Set OFDM parameters
  ofdm.set_parameters(nFFT, nCylicPrefix);

  //Iterate over all EbN0dB values:
  for (i = 0; i < alpha.length(); i++) {

    //Show how the simulation progresses:
    cout << "Now simulating alpha value = " << alpha(i); 
    cout << " # " << i + 1 << "/" << alpha.length() << endl;

    //Generate a vector of random bits to transmit:
    transmitted_bits_1 = randb(Number_of_bits);
    transmitted_bits_2 = randb(Number_of_bits);

    //convolutional code encode
    bvec out_binary_1 = encode(nsc, constraint_length, transmitted_bits_1, blockSize, false);
    bvec out_binary_2 = encode(nsc, constraint_length, transmitted_bits_2, blockSize, false);
    
    //nsc.encode_tail(transmitted_bits_1, nsc_coded_bits);//tail is added here to information bits to close the trellis

    //Modulate the bits to QPSK symbols:
    transmitted_symbols_1 = qpsk.modulate_bits(out_binary_1);
    transmitted_symbols_2 = qpsk.modulate_bits(out_binary_2);

    //Multiplex two signals
    transmitted_symbols = transmitted_symbols_1 * pow(alpha(i), 0.5) + transmitted_symbols_2 * pow(1 - alpha(i), 0.5);
    transmitted_symbols = transmitted_symbols * pow(Ps, 0.5);
    //eval_avg_power(transmitted_symbols);

    //Fading
    transmitted_symbols_1 = transmitted_symbols * pow(h1, 0.5);
    transmitted_symbols_2 = transmitted_symbols * pow(h2, 0.5);

    //OFDM modulate
    zero_pad_back(transmitted_symbols_1, 2048);
    zero_pad_back(transmitted_symbols_2, 2048);
    ofdm.modulate(transmitted_symbols_1, ofdm_symbols_1);
    ofdm.modulate(transmitted_symbols_2, ofdm_symbols_2);

    //Set the noise variance of the AWGN channel:
    awgn_channel.set_noise(N0);

    //Up-conversion

    //Run the transmited symbols through the channel using the () operator:
//    ofdm_symbols_1 = awgn_channel( multipath_channel( ofdm_symbols_1));
//    ofdm_symbols_2 = awgn_channel( multipath_channel( ofdm_symbols_2));
    ofdm_symbols_1 = awgn_channel( ofdm_symbols_1);
    ofdm_symbols_2 = awgn_channel( ofdm_symbols_2);

    // test zone!! KEEP OUT!!
    




    //alpha no greater than 0.5
    //OFDM demodulate
    ofdm.demodulate(ofdm_symbols_1, received_symbols_1);
    ofdm.demodulate(ofdm_symbols_2, received_symbols_2);

    //Demodulate the received QPSK symbols into received bits: Layer 1
    received_bits_2 = qpsk.demodulate_bits(received_symbols_2);
    bvec out_binary_recover_2 = decode(nsc, constraint_length, received_bits_2, blockSize, false);
    
    //Demodulate the received QPSK symbols into received bits: Layer 2
    received_bits_1 = qpsk.demodulate_bits(received_symbols_1);
    feedback_symbols_2 = pow(Ps * (1-alpha(i)) * h1, 0.5) * qpsk.modulate_bits(received_bits_1);
    received_bits_1 = qpsk.demodulate_bits(received_symbols_1 - feedback_symbols_2);
    bvec out_binary_recover_1 = decode(nsc, constraint_length, received_bits_1, blockSize, false);

    //Calculate the bit error rate:
    berc.clear();                               //Clear the bit error rate counter
    berc.count(transmitted_bits_2, out_binary_recover_2); //Count the bit errors
    bit_error_rate_2(i) = berc.get_errorrate();   //Save the estimated BER in the result vector

    berc.clear();
    berc.count(transmitted_bits_1, out_binary_recover_1);
    bit_error_rate_1(i) = berc.get_errorrate();
  }

  tt.toc();

  // Theoretical results for multiple access
  for (size_t i = 0; i < alpha.length(); ++i) {        //BER for theo 1
    EbN0 = (Ps * 0.5 * h1 * alpha(i)) / (N0 + Ps * 0.5 * h1 * (1 - alpha(i)));
    ber_theo_1(i) = 0.5*erfc(pow(EbN0, 0.5));
  }
  for (size_t i = 0; i < alpha.length(); ++i) {        //BER for theo 2
    EbN0 = (Ps * 0.5 * h2 * (1-alpha(i))) / N0;
    ber_theo_2(i) = 0.5*erfc(pow(EbN0, 0.5));
  }

  //Print the results:
  cout << endl;
  time_t rawtime;
  time (&rawtime);
  std::string nowTime( ctime( &rawtime));
  cout << nowTime << endl;
  cout << "alpha = " << alpha << " " << endl;
  cout << "BER 1 = " << bit_error_rate_1 << endl;
  cout << "BER 2 = " << bit_error_rate_2 << endl;
  cout << "Theoretical BER 1 = " << ber_theo_1 << endl;
  cout << "Theoretical BER 2 = " << ber_theo_2 << endl;
  cout << endl;

  //Save the results to file:
  const char xFilename[] = "result";
  char oFilename[100];
  int fileIndex = 0;
  strcpy( oFilename, xFilename);
  while(std::ifstream( oFilename)){
    ++fileIndex;
    sprintf( oFilename, "%s%d", xFilename, fileIndex);
  }
  cout << "Saving results to " << oFilename << endl;
  std::ofstream oo;
  oo.open(oFilename);
    oo << alpha << endl;
    oo << bit_error_rate_1 << endl;
    oo << bit_error_rate_2 << endl;
  oo.close();

  //Exit program:
  return 0;

}

//Evaluate average power of given vector of symbols #M00
double eval_avg_power(const cvec& symbol_vec)
{
  //Overflow check is not implemented in this function
  bool verbose = true;
  double average_power = 0.0;
  double acculmulate_power = 0.0;
  for (size_t i = 0; i < symbol_vec.length(); ++i) {
    acculmulate_power += pow(std::abs(symbol_vec(i)), 2.0);
  }
  average_power = acculmulate_power / symbol_vec.length();

  //Result [verbose]
  if (verbose) {
    cout << "[M00] " << "Average power = " << average_power << endl;
  }
  return average_power;
}

bvec encode(Convolutional_Code& nsc, int constraint_length, const bvec& encoder_input, int blockSize, bool verbose)
{
  if (verbose) {cout << "input : " << encoder_input << endl;}
  
  int codedLen = 2 * (blockSize + (constraint_length - 1));
  int nBlocks = encoder_input.length() / blockSize;
  ivec window(blockSize);
  for (int j = 0; j < blockSize; j++) {
    window[j] = j;
  }
  
  bvec nsc_coded_bits(codedLen);
  bvec tr_coded_bits;
  for (int j = 0; j < nBlocks; j++) {
    nsc.encode_tail(encoder_input(window), nsc_coded_bits);
    window = window + blockSize;
    tr_coded_bits = concat(tr_coded_bits, nsc_coded_bits);
  }
  
  // Deal with residual sources if remainder exsists
  if (nBlocks*blockSize != encoder_input.length()) {
    bvec residual_bits = encoder_input.get(nBlocks*blockSize, encoder_input.length()-1);
    nsc.encode_tail(residual_bits, nsc_coded_bits);
    tr_coded_bits = concat(tr_coded_bits, nsc_coded_bits);
  }
  
  if (verbose) {cout << "encoder output: " << tr_coded_bits << endl;}
  return tr_coded_bits;
}

bvec decode(Convolutional_Code& nsc, int constraint_length, const bvec& decoder_input, int blockSize, bool verbose)
{
  BPSK mod;
  vec decoder_input_mod = mod.modulate_bits(decoder_input);
  int codedLen = 2 * (blockSize + (constraint_length - 1));
  int nBlock_rcvd = decoder_input_mod.length() / codedLen;
  if (verbose) {cout << "rcvd number of blocks : " << nBlock_rcvd << endl;}
  
  vec codedBlock(codedLen);
  bvec bit_rcvd_tmp(blockSize);
  bvec bit_decoded;
  
  for (int j = 0; j < nBlock_rcvd; j++) {
    for (int k = 0; k < codedLen; k++) {
      codedBlock[k] = decoder_input_mod[k + j*codedLen];
    }
    //cout << codedBlock << endl;
    bit_rcvd_tmp = nsc.decode_tail(codedBlock);
    bit_decoded = concat(bit_decoded, bit_rcvd_tmp);
  }
  
  // Deal with residual sources if remainder exsists
  vec residual_bits = decoder_input_mod.get(nBlock_rcvd*codedLen, decoder_input_mod.length()-1);
  nsc.decode_tail( residual_bits, bit_rcvd_tmp);
  bit_decoded = concat(bit_decoded, bit_rcvd_tmp);
  
  if (verbose) {cout << "decode output : " << bit_decoded << endl;}
  return bit_decoded;
}

void zero_pad_back(cvec& vec_in, int m_in)
{
  int quotient = vec_in.length() / m_in;
  int n_remainder = (quotient + 1) * m_in - vec_in.length();
  vec_in = concat(vec_in, zeros_c(n_remainder));
}





