/*
Title       : Experiment 2
Project     : Mtk 5th generation communication system
Author      : Ming Jie Yang
Date        : 
Description : 
              
              
              
              
              
*/
#include <itpp/itcomm.h>
#include <ctime>

using namespace itpp;

//These lines are needed for use of cout and endl
using std::cout;
using std::endl;

//Function defined here
double eval_avg_power(const cvec& symbol_vec);

int main()
{
  //Declarations of scalars and vectors:
  int i, Number_of_bits;
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
  cvec received_symbols_1, received_symbols_2, feedback_symbols_1;
  cvec ofdm_symbols_1, ofdm_symbols_2;

  //Declarations of classes:
  QPSK qpsk;                     //The QPSK modulator class
  AWGN_Channel awgn_channel;     //The AWGN channel class
  it_file ff;                    //For saving the results to file
  BERC berc;                     //Used to count the bit errors
  Real_Timer tt;                 //The timer used to measure the execution time
  OFDM ofdm;                     //OFDM modulator class

  //Reset and start the timer:
  tt.tic();

  //Init:
  Ps = 5 * pow(10, -3);          //The transmitted energy per QPSK symbol is 1, 5e-3 W/MHz
  N0 = 4 * pow(10, -15);         //Thermal noise -144dBm/MHz
  dist_1 = 150;                   //Distance form transmitter to receiver 1 (meter)
  dist_2 = 220;                  //                   *                  2
  h1 = pow(dist_1, -4.5);        //Channel gain for receiver 1
  h2 = pow(dist_2, -4.5);        //
  cout << "h1 = " << h1 << " h2 = " << h2 << endl;
  alpha = linspace(1.0, 0.0, 21);//Simulate for different weight on power

//  Eb = Ec / 2.0;                 //The transmitted energy per bit is 0.5.
//  EbN0dB = linspace(0.0, 9.0, 10); //Simulate for 10 Eb/N0 values from 0 to 9 dB.
//  EbN0 = inv_dB(EbN0dB);         //Calculate Eb/N0 in a linear scale instead of dB.
//  N0 = Eb * pow(EbN0, -1.0);     //N0 is the variance of the (complex valued) noise.
  Number_of_bits = 204800;       //One hundred thousand bits is transmitted for each Eb/N0 value
//  alpha = 0.8;                   //Ratio of allocated power
  nFFT = 2048;
  nCylicPrefix = 144;

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

    //Modulate the bits to QPSK symbols:
    transmitted_symbols_1 = qpsk.modulate_bits(transmitted_bits_1);
    transmitted_symbols_2 = qpsk.modulate_bits(transmitted_bits_2);

    //Multiplex two signals
    transmitted_symbols = transmitted_symbols_1 * pow(alpha(i), 0.5) + transmitted_symbols_2 * pow(1 - alpha(i), 0.5);
    transmitted_symbols = transmitted_symbols * pow(Ps, 0.5);
    //eval_avg_power(transmitted_symbols);

    //Fading
    transmitted_symbols_1 = transmitted_symbols * pow(h1, 0.5);
    transmitted_symbols_2 = transmitted_symbols * pow(h2, 0.5);

    //OFDM modulate
    ofdm.modulate(transmitted_symbols_1, ofdm_symbols_1);
    ofdm.modulate(transmitted_symbols_2, ofdm_symbols_2);

    //Set the noise variance of the AWGN channel:
    awgn_channel.set_noise(N0);

    //Run the transmited symbols through the channel using the () operator:
    ofdm_symbols_1 = awgn_channel(ofdm_symbols_1);
    ofdm_symbols_2 = awgn_channel(ofdm_symbols_2);

    //OFDM demodulate
    ofdm.demodulate(ofdm_symbols_1, received_symbols_1);
    ofdm.demodulate(ofdm_symbols_2, received_symbols_2);

    //Demodulate the received QPSK symbols into received bits: Layer 1
    received_bits_1 = qpsk.demodulate_bits(received_symbols_1);

    //Demodulate the received QPSK symbols into received bits: Layer 2
    received_bits_2 = qpsk.demodulate_bits(received_symbols_2);
    feedback_symbols_1 = pow(Ps * alpha(i) * h2, 0.5) * qpsk.modulate_bits(received_bits_2);
    received_bits_2 = qpsk.demodulate_bits(received_symbols_2 - feedback_symbols_1);

    //Calculate the bit error rate:
    berc.clear();                               //Clear the bit error rate counter
    berc.count(transmitted_bits_1, received_bits_1); //Count the bit errors
    bit_error_rate_1(i) = berc.get_errorrate();   //Save the estimated BER in the result vector

    berc.clear();
    berc.count(transmitted_bits_2, received_bits_2);
    bit_error_rate_2(i) = berc.get_errorrate();
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
  cout << "alpha = " << alpha << " " << endl;
  cout << "BER 1 = " << bit_error_rate_1 << endl;
  cout << "BER 2 = " << bit_error_rate_2 << endl;
  cout << "Theoretical BER 1 = " << ber_theo_1 << endl;
  cout << "Theoretical BER 2 = " << ber_theo_2 << endl;
  cout << "Saving results to ./qpsk_result_file.it" << endl;
  cout << endl;

  //Save the results to file:
  time_t rawtime;
  time (&rawtime);
  std::string nowTime( ctime( &rawtime));
  cout << nowTime << endl;
  std::ofstream oo;
  oo.open(ctime( &rawtime));
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





