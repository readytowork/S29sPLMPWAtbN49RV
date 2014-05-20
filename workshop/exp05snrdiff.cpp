/*
Title       : Experiment 2
Project     : Mtk 5th generation communication system
Author      : Ming Jie Yang
Date        : 
Description : 
              
              
              
              
              
*/
#include <itpp/itcomm.h>

using namespace itpp;

//These lines are needed for use of cout and endl
using std::cout;
using std::endl;

//Function defined here
double eval_avg_power(const cvec& symbol_vec);
void FDE(cvec& signal_in, cvec& signal_out, vec& imp_response, int ifftSize, bool verbose);

int main()
{
  //Declarations of scalars and vectors:
  int i, Number_of_bits;
  double Ps, N0, dist_1, dist_2, h1, h2, Eb;
  double EbN0;
  double rcvd_power_1, rcvd_power_2;

  vec alpha;                                            //vec is a vector containing double
  vec bit_error_rate_1, bit_error_rate_2;
  vec ber_theo_1, ber_theo_2;                           // Theoretical results for multiple access

  bvec transmitted_bits_1, transmitted_bits_2;          //bvec is a vector containing bits
  bvec received_bits_1, received_bits_2;

  cvec transmitted_symbols_1, transmitted_symbols_2, transmitted_symbols; //cvec is a vector containing double_complex
  cvec received_symbols_1, received_symbols_2, feedback_symbols_1;

  //Declarations of classes:
  QPSK qpsk;                     //The QPSK modulator class
  AWGN_Channel awgn_channel;     //The AWGN channel class
  it_file ff;                    //For saving the results to file
  BERC berc;                     //Used to count the bit errors
  Real_Timer tt;                 //The timer used to measure the execution time

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
  Number_of_bits = 1000000;       //One hundred thousand bits is transmitted for each Eb/N0 value
//  alpha = 0.8;                   //Ratio of allocated power

  //Allocate storage space for the result vector.
  //The "false" argument means "Do not copy the old content of the vector to the new storage area."
  bit_error_rate_1.set_size(alpha.length(), false);
  bit_error_rate_2.set_size(alpha.length(), false);

  //Storage for theoretical ber
  ber_theo_1.set_size(alpha.length(), false);
  ber_theo_2.set_size(alpha.length(), false);

  //Randomize the random number generators in it++:
  RNG_randomize();

  bool verbose = false;

  //Iterate over all EbN0dB values:
  for (i = 0; i < alpha.length(); i++) {

    //Show how the simulation progresses:
    cout << "Now simulating alpha value = " << alpha(i); 
    cout << " # " << i + 1 << "/" << alpha.length() << endl;

    //Generate a vector of random bits to transmit:
    transmitted_bits_1 = randb(Number_of_bits);
    transmitted_bits_2 = randb(Number_of_bits);

    ivec getValIndex = "05 11 12 07 03";
    cout << getValIndex << endl;
    cout << transmitted_bits_1(getValIndex) << endl;
    for (int j = 0; j < 5; j++) {
      cout << transmitted_bits_1(getValIndex(j)) << " ";
    }
    cout << endl;
    ivec sorted = sort_index(getValIndex);
    cout << sorted << endl;
    cout << getValIndex(sorted) << endl;
    cout << getValIndex(sorted) + 1 << endl;

    // convolutional code
    if (verbose) {
    
    ivec gen = "01 07";
    int constraint_length = 3;
    Convolutional_Code nsc;
    nsc.set_generator_polynomials(gen, constraint_length);
    int blockSize = 3; // input size
    bvec encoder_input = randb(200);
    cout << "input : " << encoder_input << endl;
    int codedLen = 2 * (blockSize + (constraint_length - 1));
    int nBlocks = encoder_input.length() / blockSize;
    ivec window(blockSize);
    for (int j = 0; j < blockSize; j++) {
      window[j] = j;
    }
    bvec nsc_coded_bits(codedLen);
    //cout << codedLen << endl;
    bvec tr_coded_bits;
    //cout << window << endl;
    for (int j = 0; j < nBlocks; j++) {
      cout << encoder_input(window) << endl;
      nsc.encode_tail(encoder_input(window), nsc_coded_bits);
      window = window + blockSize;
      tr_coded_bits = concat(tr_coded_bits, nsc_coded_bits);
    }
    //cout << tr_coded_bits.length() << endl;
    //cout << (Number_of_bits/256)*(512+4) << endl;
    /*for (int j = 0; j < Number_of_bits; j++) {
      transmitted_bits_2[j] = tr_coded_bits[j];
    }*/
    cout << "encoder output: " << tr_coded_bits << endl;
    
    BPSK mod;
    vec rec_sig = mod.modulate_bits(tr_coded_bits);
    
    vec decoder_input = rec_sig;
    int nBlock_rcvd = decoder_input.length() / codedLen;
    //cout << "rcvd number of blocks : " << nBlock_rcvd << endl;
    vec codedBlock(codedLen);
    bvec bit_rcvd_tmp(blockSize);
    bvec bit_decoded;
    for (int j = 0; j < nBlock_rcvd; j++) {
      for (int k = 0; k < codedLen; k++) {
        codedBlock[k] = decoder_input[k + j*codedLen];
      }
      cout << codedBlock << endl;
      bit_rcvd_tmp = nsc.decode_tail(codedBlock);
      bit_decoded = concat(bit_decoded, bit_rcvd_tmp);
    }
    cout << "decode output : " << bit_decoded << endl;
    double nscErr = 0;
    for (int j = 0; j < bit_decoded.length(); j++) {
      if (bit_decoded[j] != encoder_input[j]) {
        nscErr++;
      }
    }
    cout << nscErr / bit_decoded.length() << endl;
    cout << encoder_input.length() << " ; " << bit_decoded.length() << endl;
    }
    
    // linear mmse test
    if (verbose) {
    
    cvec sig_test = "9.07032e-08-1.03701e-07i 3.13527e-07-7.27227e-08i -1.67824e-07-4.22301e-07i 1.3437e-07-2.14683e-08i -1.5889e-07+1.45796e-07i -2.84262e-07+4.55597e-08i 2.67125e-07+5.13786e-07i -2.85874e-07+1.4919e-07i";
    sig_test = concat(sig_test, zeros_c(4));
    vec ht = "1.0 1.0 0.0";
    int ch_tap = ht.length();
    MA_Filter<std::complex<double>, std::complex<double>, std::complex<double> > multipath_channel;
    multipath_channel.set_coeffs( to_cvec(ht));
    multipath_channel.set_state(to_cvec( zeros(ch_tap)));
    cvec sig_rcvd = multipath_channel(sig_test);
    cout << sig_rcvd << endl;
    // input ht vec
    double noisevar_ = 4 * pow(10, -15);
    int eq_tap = 7;
    cout << "ht : " << ht << endl;
    
    vec hAutoCorr(ch_tap * 2 -1);
    vec tmp(ch_tap);
    for (int j = 0; j < ch_tap; j++) {
      tmp = ht * ht[ch_tap - j - 1];
      for (int k = 0; k < ch_tap; k++) {
        hAutoCorr[j+k] += tmp[k];
      }
    }
    cout << "hAutoCorr : " << hAutoCorr << endl;
    
    mat hM(eq_tap, eq_tap);
    int idx;
    for (int j = 0; j < eq_tap; j++) {
      for (int k = 0; k < eq_tap; k++) {
        idx = ch_tap + k - j -1;
        if (idx >= ch_tap * 2 -1 || idx < 0) {
          hM(j, k) = 0;
        } 
        else {
          hM(j, k) = hAutoCorr(ch_tap + k - j -1);
        }
      }
    }
    cout << "hM : \n" << hM << endl;
    
    mat eye_7 = eye(eq_tap);
    hM += eye_7 * noisevar_;
    cout << "puls N0: \n" << hM << endl;
    
    vec c_mmse(eq_tap);
    vec d = concat(zeros(2), reverse(ht), zeros(2));
    #if defined(HAVE_LAPACK)
    cout << "YES\n";
    #endif
    mat inv_hM = inv(hM);
    cout << "puls N0: \n" << hM << endl;
    cout << "INV:\n" << inv_hM << endl;
    for (int j = 0; j < eq_tap; j++) {
      c_mmse[j] = dot(inv_hM.get_col(j),d);
    }
    cout << "c MMSE: " << c_mmse << endl;
    
    MA_Filter<std::complex<double>, std::complex<double>, std::complex<double> > equalizer;
    equalizer.set_coeffs( to_cvec(c_mmse));
    equalizer.set_state(to_cvec( zeros(eq_tap)));
    cvec sig_eqed = equalizer(sig_rcvd);
    cvec ans = sig_eqed.get(4, sig_eqed.length()-1);
    cout << "ANS: " << ans << endl;
    
    } // END
    
    //FDE
    if (verbose) {
    
    MA_Filter<std::complex<double>, std::complex<double>, std::complex<double> > testchannel;
    vec ht = "1 0.9";
    testchannel.set_coeffs( to_cvec(ht));
    testchannel.set_state(to_cvec( zeros(2)));
    cvec testCase = "1 1 1 1 1 0";
    cout << testchannel(testCase) << endl;
    
    OFDM ofdm;
    ofdm.set_parameters(64,16,1);
    cvec input_symbol = "0 0 0 0 0 0 1 1 1 -1 -1 -1 -1 1 -1 -1 -1 -1 -1 1 -1 -1 -1 -1 -1 -1 1 1 1 -1 -1 1 0 -1 -1 1 1 1 1 -1 1 -1 -1 1 1 1 1 1 -1 -1 -1 1 1 -1 -1 1 1 1 1 0 0 0 0 0";
    cvec symbol_t;
    cout << "input_symbol:\n" << input_symbol << endl;
    //input_symbol = concat(input_symbol.mid(input_symbol.length()/2, input_symbol.length()/2), input_symbol.mid(0, input_symbol.length()/2)); // fftshift
    ofdm.modulate(input_symbol, symbol_t);
    symbol_t = concat(symbol_t, zeros_c(1));
    cout << "symbol_t:\n" << symbol_t << endl;
    cout << symbol_t.length() << endl;
    cvec sig_rcvd_t = testchannel(symbol_t);
    //sig_rcvd_t = sig_rcvd_t.get(0, sig_rcvd_t.length()-2);  //MARk
    cout << "channel:\n" << sig_rcvd_t << endl << sig_rcvd_t.length() << endl;
    sig_rcvd_t = sig_rcvd_t.get(0, 79);
    cout << "channel:\n" << sig_rcvd_t << endl;
    
    cvec para = fft(to_cvec(concat(ht, zeros(62))));
    //para = concat(para.mid(32, 32), para.mid(0, 32));
    cout << "para\n" << para << endl;
    cvec symbol_t_de;
    ofdm.demodulate(sig_rcvd_t, symbol_t_de);
    //symbol_t_de = concat(symbol_t_de.mid(32, 32), symbol_t_de.mid(0, 32));
    cout << "demdul:\n" << symbol_t_de << endl;
    for (int j = 0; j < 64; j++) {
      symbol_t_de[j] /= para[j];
    }
    cout << "ans:" << symbol_t_de << endl;
    
    } // End FDE
    
    // ITU Pedestrian A channel
    if (verbose)
    {
    
    double Ts = 3.25521e-08;
    Channel_Specification channel_spec(ITU_Pedestrian_A);
    TDL_Channel my_channel(channel_spec, Ts);
    my_channel.set_norm_doppler(0.01);
    my_channel.set_doppler_spectrum(0, Jakes);
    cmat coeff;
    cout << my_channel.taps() << endl;
    my_channel.generate(1,coeff);
    cout << coeff << endl;
    
    } // end ITU channel model
    
    // FDE fit in
    if (verbose)
    {
    
    cout << "FDE fit in\n";
    cvec a = "1 2 3 4 5 6 7 8";
    a = randn_c(2048);
    cvec a_ofdm;
    OFDM ofdm;
    ofdm.set_parameters(2048, 144);
    cout << a << endl;
    ofdm.modulate(a, a_ofdm);
    cout << a_ofdm << endl;
    vec ht = "1.000000000000000 0.436232649307735 0.198694326007661";
    int ch_tap = ht.length();
    MA_Filter<std::complex<double>, std::complex<double>, std::complex<double> > multipath_channel;
    multipath_channel.set_coeffs( to_cvec(ht));
    multipath_channel.set_state(to_cvec( zeros(ch_tap)));
    cout << a_ofdm << endl;
    cvec a_ofdm_mult = multipath_channel(a_ofdm);
    cout << a_ofdm_mult << endl;
    cvec b;
    ofdm.demodulate(a_ofdm_mult, b);
    cout << b << endl;
    cvec signal_out;
    FDE(b, signal_out, ht, 2048, true);
    cout << eval_avg_power(a-signal_out)/eval_avg_power(a) << endl;
    } // END FDE fit in
    
    
	Punctured_Convolutional_Code code;
	ivec generator(2);
	generator(0)=05;
	generator(1)=02;
	//generator(2)=0171;
	code.set_generator_polynomials(generator, 2);
	bmat puncture_matrix = "1 0;1 1";
	code.set_puncture_matrix(puncture_matrix);
	//code.set_truncation_length(30);
    cout << code.get_rate() << "OKAY\n";
    bvec testb = "0 1 0 1 0 1";
    bvec testb_en = code.encode_tail(testb);
    cout << testb_en << endl;
    cout << code.decode_tail(to_vec(testb_en)) << endl;
    
    QAM qam16(16);
    bvec testbin = "0 0 0 1 0 1 1 0 0 0 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1";
    cvec testo = qam16.modulate_bits(testbin);
    bvec rcvdtest = qam16.demodulate_bits(testo);
    cout << rcvdtest << endl;
    return 0;

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

    //Set the noise variance of the AWGN channel:
    awgn_channel.set_noise(N0);

    //Run the transmited symbols through the channel using the () operator:
    received_symbols_1 = awgn_channel(transmitted_symbols_1);
    received_symbols_2 = awgn_channel(transmitted_symbols_2);

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
    EbN0 = (Ps * h1 * alpha(i)) / (N0 + Ps * h1 * (1 - alpha(i)));
    ber_theo_1(i) = 0.5*erfc(pow(EbN0, 0.5));
  }
  for (size_t i = 0; i < alpha.length(); ++i) {        //BER for theo 2
    EbN0 = (Ps * h2 * (1-alpha(i))) / N0;
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
  std::ofstream oo;
  oo.open("exp05result.txt");
  for (size_t i = 1; i < alpha.length(); i++) {
    oo << alpha << " ";
    oo << bit_error_rate_1 << " ";
    oo << bit_error_rate_2 << endl;
  }
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

void FDE(cvec& signal_in, cvec& signal_out, vec& imp_response, int ifftSize, bool verbose)
{
  cvec para = fft(to_cvec(concat(imp_response, zeros(ifftSize-imp_response.length()))));
  if (verbose) {cout << "para: \n" << para << endl;}
  
  int nSymbols = signal_in.length() / ifftSize;
  it_assert(nSymbols * ifftSize == signal_in.length(), "warning");
  
  signal_out = zeros_c(signal_in.length());
  int IdxSig = 0;
  for (int i = 0; i < nSymbols; i++) {
    for (int j = 0; j < ifftSize; j++) {
      signal_out[IdxSig] = signal_in[IdxSig] / para[j];
      IdxSig++;
    }
  }
}






