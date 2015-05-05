#include <itpp/itcomm.h>
#include <itpp/stat/misc_stat.h>
#include <math.h>

// OFDM constants
#define FFTSIZE 64
#define CP 7
#define NSYMBOLS 5

// Will run for near user, by default.
#ifndef NEAR
    #ifndef FAR
        #define NEAR
    #endif
#endif

// Modulation constants
#define MAX_QAM_BITS_PER_SYMBOL 4
#define nModulationSchemes 3

//Channel parameters
#define nTaps 5
#define conerence_nrof_symbol 40

//Deciding the number of bits to be generated.
//FFTSIZE * NSYMBOLS must be a factor of the no. of bits.
//For 16-Qam to work exhaust all these bits, 4 * FFTSIZE * NSYMBOLS must be a factor of the no. of bits.
//To exhaust all the bits generated after convolutional coding,
// 1) 2 * no of bits should be a multiple of 4 * FFTSIZE * NSYMBOLS and
// 2) 3/2 * no of bits should be a multiple of 4 * FFTSIZE * NSYMBOLS and
// 3) 4/3 * no of bits should be a multiple of 4 * FFTSIZE * NSYMBOLS and
// 4) 6/5 * no of bits should be a multiple of 4 * FFTSIZE * NSYMBOLS.

// Using specific values, 4 * 64 * 50 = 12800
// Choose the number to be (LCM of 1,2,3,5)*12800, i.e., 12800 * 30 = 384000

#define nrof_samples  (FFTSIZE * NSYMBOLS * MAX_QAM_BITS_PER_SYMBOL * 30)/conerence_nrof_symbol

// Punctured convolutional code constants
#define nConvoCodes 4

using namespace itpp;
using namespace std;

int main(void)
{
    if(NSYMBOLS>100)
    {
        cout << "Too many symbols" << endl;
        exit(1);
    }

    // Punctured convolutional code
    Punctured_Convolutional_Code conv_code[nConvoCodes];
    int constraint_length[nConvoCodes] = {7, 7, 7, 7};
    ivec generators[nConvoCodes] = {"0133 0171", "0133 0171", "0133 0171", "0133 0171"};
    bmat puncture_matrix[4] = {"1; 1", "1 1; 1 0", "1 1 0; 1 0 1","1 1 0 1 0; 1 0 1 0 1"};

    //OFDM
    OFDM ofdm(FFTSIZE, CP);
    cvec ofdm_modulated_symbols, channel_fft_estimate;

    //Modulation schemes
    BPSK_c bpsk;
    QAM qam_modulate;

    //Simulation parameters
    int MaxIterations;
    vec Ec, Eb;

    // Pilots
    cvec ones_encode;

    //Transmission-reception
    bvec uncoded_bits, coded_bits, decoded_bits, rx_bits;
    vec EbN0dB, EbN0, build, soft_demod_result, soft_demod_result1;
    double N0;
    cvec trans_symbols, rec_symbols, convolved_symbols, ofdm_demodulated, ofdm_demodulated_coherence;

    //Channels and Coherence time parameters
    AWGN_Channel channel;
    vec my_multipath_channel(nTaps);
    cvec my_multipath_channel_fft;
    double a_F_dB = 2;

    long int coherence_time = FFTSIZE * conerence_nrof_symbol;
    long int channel_number = 0;
    cvec channel_fft;
    cmat ch_coeffs;
    TDL_Channel my_channel;
    double Ts = 0.0000001;

    //Bit-error-rates and their holders
    BERC berc;
    vec ber;

    //iterators
    int p, i, j, modulations_count, nMod, nCode, t=0, k;
    const char *modulations[3];
    modulations[0] = "BPSK";
    modulations[1] = "QPSK";
    modulations[2] = "16QAM";

    //Output
    vec P_min;
    P_min.set_size(12, false);

    //Setting the simulation parameters
    EbN0dB = linspace(-6, 40, 40); // Set SNR steps to plot the BERvsSNR curve.
    EbN0 = inv_dB(EbN0dB);
    RNG_randomize();
    N0 = 1;

    // If it is near user,
    // generate and save a new channel from a Standard Normal distribution.
    #ifdef NEAR
        RNG_randomize();
	    my_channel.set_channel_profile(ITU_Vehicular_A, Ts);
	    my_channel.set_norm_doppler(0.01);
	    my_channel.generate(nrof_samples, ch_coeffs);
	    it_file channel_file("multipath_channel.it");
	    channel_file << Name("multipath_channel") << ch_coeffs;
	    channel_file.close();
	// Else if it is far user,
	// load the already generated channel
    #else
	    it_file channel_file("multipath_channel.it");
	    channel_file >> Name("multipath_channel") >> ch_coeffs;
	    channel_file.close();
    #endif

    Eb = N0 * EbN0;
    ber.set_size(EbN0dB.length(), false);

    channel.set_noise(N0 / 2.0); //Set the noise value of the AWGN channel.

    // Prepare the pilot.
    ones_encode.set_size(FFTSIZE);
    ones_encode = ones_c(FFTSIZE);

    my_multipath_channel_fft = fft_real(my_multipath_channel, FFTSIZE);

    // Simulation starts here
    for(nCode = 0; nCode < nConvoCodes; nCode++)
    {
        // Set up nCode'th Convolutional code.
        conv_code[nCode].set_generator_polynomials(generators[nCode], constraint_length[nCode]);
        conv_code[nCode].set_puncture_matrix(puncture_matrix[nCode]);

        // Normalize the energy
        Ec = Eb * conv_code[nCode].get_rate();
        modulations_count = 0;

        for(i = 0; i < 5; i += 2)
        {
            // Set up (clear) the BER vector to hold the fresh values.
            ber.clear();

            // If either QPSK or 16-QAM,
            if(i)
            {
                qam_modulate.set_M(pow2i(i));
                Ec = Ec * qam_modulate.get_k(); // Normalize energy
            }

            for(p = 0; p < EbN0dB.length(); p++) // Run for all the steps in SNR.
            {
                // Clear the counter.
                berc.clear();

                // Generate the data.
                uncoded_bits = randb(Nobits);         // The uncoded bits.
                coded_bits = conv_code[nCode].encode(uncoded_bits);

                // BPSK case
                if(!i)
                {
                    // Modulation
                    bpsk.modulate_bits(coded_bits, trans_symbols); //The BPSK modulator
                    channel_number = 0;
                    ofdm_modulated_symbols.clear();

                    convolved_symbols.clear();
                    convolved_symbols.set_size(0);
                    for(j = 0; j < trans_symbols.length(); j = j + coherence_time)
                    {
                        // Building the j'th slot of the coherence time.
                        ofdm_modulated_symbols = sqrt(Ec(p)) * ofdm.modulate(concat(ones_encode, trans_symbols(j, j + coherence_time - 1)));
                        // Passing it through the multipath channel and collecting the output.
                        convolved_symbols = concat(convolved_symbols, filter(ch_coeffs.get_row(channel_number++), 1, ofdm_modulated_symbols));
                    }
                    
                    #ifdef FAR
                        convolved_symbols /= sqrt(inv_dB(a_F_dB)); // Relative attenuation
                    #endif
                    rec_symbols = channel(convolved_symbols);

                    // Demodulation
                    ofdm_demodulated.set_size(rec_symbols.length());
                    ofdm_demodulated = ofdm.demodulate(rec_symbols) / sqrt(Ec(p));
                    soft_demod_result.clear();
                    soft_demod_result.set_size(0,false);
                    channel_number = 0;

                    for(k = 0; k < ofdm_demodulated.length(); k = k + coherence_time + FFTSIZE)
                    {
                        // Obtain the channel estimate for the k'th slot of the coherence time.
                        channel_fft_estimate = (ofdm_demodulated(k, k + FFTSIZE - 1)) / (1 + (N0 / (2 * Ec(p)))); 
                        ofdm_demodulated_coherence = ofdm_demodulated(k + FFTSIZE, k + coherence_time + FFTSIZE - 1);
 
                        // Generate soft demodulated bits for all the FFT blocks separately and concatenate them.
                        for(j = 0; j < coherence_time; j += FFTSIZE)
                        {
                            build = bpsk.demodulate_soft_bits(ofdm_demodulated_coherence(j, j + FFTSIZE - 1), channel_fft_estimate, N0 / (2.0 * Ec(p)));
                            soft_demod_result = concat(soft_demod_result, build);
                        }
                  }
                }
                
                // QPSK or 16-QAM
                else
                {
                    // Modulation                
                    qam_modulate.modulate_bits(coded_bits, trans_symbols);
                    channel_number = 0;

                    ofdm_modulated_symbols.clear();

                    convolved_symbols.clear();
                    convolved_symbols.set_size(0);
                    for(j = 0;j < trans_symbols.length(); j = j + coherence_time)
                    {
                        // Building the j'th slot of the coherence time.                    
                        ofdm_modulated_symbols = sqrt(Ec(p)) * ofdm.modulate(concat(ones_encode, trans_symbols(j, j + coherence_time - 1)));
                        // Passing it through the multipath channel and collecting the output.                       
                        convolved_symbols = concat(convolved_symbols, filter(ch_coeffs.get_row(channel_number++), 1, ofdm_modulated_symbols));
                    }

                    #ifdef FAR
                        convolved_symbols /= sqrt(inv_dB(a_F_dB)); // Relative attenuation
                    #endif
                    rec_symbols = channel(convolved_symbols);

                    // Demodulation
                    ofdm_demodulated.set_size(rec_symbols.length());
                    ofdm_demodulated = ofdm.demodulate(rec_symbols) / sqrt(Ec(p));
                    soft_demod_result.clear();
                    soft_demod_result.set_size(0, false);
                    channel_number = 0;

                    for(k = 0; k < ofdm_demodulated.length(); k = k + coherence_time + FFTSIZE)
                    {
                        // Obtain the channel estimate for the k'th slot of the coherence time.                    
                        channel_fft_estimate = (ofdm_demodulated(k, k + FFTSIZE - 1)) / (1 + (N0 / (2 * Ec(p))));
                        ofdm_demodulated_coherence = ofdm_demodulated(k + FFTSIZE, k + coherence_time + FFTSIZE - 1);
                        // Generate soft demodulated bits for all the FFT blocks separately and concatenate them.                        
                        for(j = 0; j < coherence_time; j += FFTSIZE)
                        {
                            build=qam_modulate.demodulate_soft_bits(ofdm_demodulated_coherence(j , j + FFTSIZE - 1), channel_fft_estimate, N0 / (2.0 * Ec(p)));
                            soft_demod_result=concat(soft_demod_result, build);
                        }
                    }
                }
                
                decoded_bits = conv_code[nCode].decode(soft_demod_result); // The Viterbi decoder  function.
                berc.count(uncoded_bits, decoded_bits); // Count the errors.
                ber(p) = berc.get_errorrate();

                if(ber(p) < 0.001) // Check feasibility.
                {
                    cout << "Clear point"<< endl;
                    P_min[t++] = Ec(p); // Store as P_min
                    cout << "P_min for rate " << conv_code[nCode].get_rate() << " and " <<  modulations[modulations_count++] << " modulation scheme" << " = "  << P_min[t - 1] << endl;
                    break;
                }
            }
        }
    }

    // Sort the P_min vector in ascending order.
    double temp;

	for(j = 1; (j <= P_min.length() - 1); j++)
	{
	    temp = P_min[j];
	    i = j - 1;
        while((i >= 0) && (P_min[i] > temp))
        {
		    P_min[i + 1] = P_min[i];
 		    i--;
        }
 		P_min[i + 1] = temp;
    }

    cout << "The 12 minimum powers to support the corresponding rates are " << endl << P_min << endl;
    cout << "The multipath channel used was " << endl << my_multipath_channel << endl;

    // Save the P_min vector appropriately.
    #ifdef FAR
        it_file Pmin_file("P_min_far.it");
        Pmin_file << Name("P_min_far") << P_min;
    #else
        it_file Pmin_file("P_min_near.it");
        Pmin_file << Name("P_min_near") << P_min;
    #endif
    Pmin_file.close();
    
    //Exit program:
    return 0;
}

