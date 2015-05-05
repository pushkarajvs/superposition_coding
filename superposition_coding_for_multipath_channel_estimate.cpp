#include <itpp/itcomm.h>
#include <itpp/stat/misc_stat.h>
#include <math.h>
#include<iostream>
using namespace itpp;
using namespace std;
#define nModulationSchemes 3
#define nConvoCodes 4
#define nUsers 2
#define nTaps 5
//OFDM constants
#define FFTSIZE 64
#define CP 7
#define NSYMBOLS 5
#define a_F_dB 2 
#define conerence_nrof_symbol 40
#define nrof_samples  (FFTSIZE * NSYMBOLS * MAX_QAM_BITS_PER_SYMBOL * 30)/conerence_nrof_symbol
//Modulation constants
#define MAX_QAM_BITS_PER_SYMBOL 4
#define Nobits ((FFTSIZE * NSYMBOLS * MAX_QAM_BITS_PER_SYMBOL * 30)-6)
//vec h_near="1.0 0.8 0.6 0.4";
//vec h_far="1.0 0.8 0.7 0.6 0.4";
double P_N=1;
#define nRATES 12


cvec Encode_users(int mod_number, int conv_number, double input_power, double alpha,  bvec block_input);
bvec Decode_users(int mod_number, int conv_number, double input_power, double alpha,  cvec block_input);
void calculate_BER(double alpha, int k, int i_F, double pmin_near, double pmin_far, cmat ch_coeffs, cmat ch_coeffs_fft, double *ber_near, double * ber_far);

int main(void)
{
    int i,k,l,found_pair=0,it;
    double power_for_far_user, ber_for_far_user, ber_for_near_user,alpha_dB,eta=1, a_F;
    double  ber_far, ber_near ;

    vec Pmin(nRATES),alpha(nRATES),beta(nRATES),i_F(nRATES),Pmin_far(nRATES);
    vec i_F_star="0 0 0 0 0 0 0 0 0 0 0 0";
    vec near_user_rates="0.50 0.67 0.75 0.83 1.00 1.33 1.50 1.67 2.00 2.67 3.00 3.33";
    mat rate_pairs(2,nRATES);
    cmat ch_coeffs_fft, ch_coeffs;

    //Taking minimum power vector of all modulation
    it_file pmins_file("P_min_near.it");
    pmins_file >> Name("P_min_near") >> Pmin;
    pmins_file.close();
  
    //Taking channel coefficient values from single user characterization
    it_file channel_file("multipath_channel.it");
    channel_file >> Name("multipath_channel") >> ch_coeffs;
    channel_file.close();

      //Calculating fraction of power allocated to near user
    for(i=0;i<nRATES;i++)
        beta[i]=Pmin[i]/Pmin[nRATES-1];

    it=0;
    
    for(k=0;k<nRATES;k++)
    {
        alpha[k]=beta[k]; //Near user fraction
        power_for_far_user=((1-beta[k])*Pmin[nRATES-1])/inv_dB(a_F_dB);

        //Determine i_F[k]
        for(l=0;l<nRATES;l++)
        {
            if(Pmin[l] > power_for_far_user)
            {
                if(l==0)
                    i_F[k]=-1;
                else
                    i_F[k]=l-1;
                break;
            }
        }

        while(1)
        {
            if(i_F[k]<0)
            {
                cout << "Cannot support far user at near user rate number " << k << endl;
                break;
            }
 
	    calculate_BER(alpha[k], k, i_F[k], Pmin[nRATES-1], Pmin[nRATES-1],  ch_coeffs, ch_coeffs_fft,  &ber_for_near_user, &ber_for_far_user);
     
           if((ber_for_near_user<0.001) && (ber_for_far_user<0.001))
            {
                i_F_star[k]=i_F[k];
                cout << "Found rate region point index " << "(" << k <<" "<< i_F[k]<< ")" << endl;
                rate_pairs(0,it)=near_user_rates[k];
                rate_pairs(1,it++)=near_user_rates(i_F_star[k]);
                break;
            }
            else if((ber_for_near_user<0.001) && (ber_for_far_user>0.001))
                i_F[k]--;
            else if((ber_for_near_user>0.001) && (ber_for_far_user<0.001))
            {
                alpha_dB=dB(alpha[k]);
                alpha_dB+=eta;
                alpha[k]=inv_dB(alpha_dB);
            }
            else
            {
                i_F[k]--;
                alpha_dB=dB(alpha[k]);
                alpha_dB+=eta;
                alpha[k]=inv_dB(alpha_dB);
            }
        }
    }

    cout << "The rate_pairs matrix is " << rate_pairs << endl;
    it_file rate_pairs_file("rate_pairs_due_to_spc.it");
    rate_pairs_file << Name("rate_pairs_due_to_spc") << rate_pairs;
    rate_pairs_file.close();
    
    return 0;

}

void calculate_BER(double alpha, int k, int i_F, double pmin_near, double pmin_far, cmat ch_coeffs, cmat ch_coeffs_fft,  double *ber_near,  double *ber_far)
{ 

    //calculating modulation number (0 --> BPSK, 1 --> QPSK, 2--> 16QAM)
    //calculating Convolutional code number (0 --> rate 1/2, 1 -->rate 2/3, 2 --> rate 3/4, 3 --> rate 5/6)
    int nMod_scheme_near_user = k/nConvoCodes;
    int nConvo_code_near_user=k%nConvoCodes;
    int nMod_scheme_far_user=i_F/nConvoCodes;
    int nConvo_code_far_user=i_F%nConvoCodes;
    
    
    long int  length_ofdm_near ,  length_ofdm_far;  
    int p, i,nMod, channel_number, j;
    long int coherence_time = FFTSIZE *conerence_nrof_symbol;

    AWGN_Channel channel_near, channel_far;

    //Bit-error-rates and their holders
    BERC berc;
    vec ber;
    
    // Transmission-Reception and pilot parameters 
    cvec zeros_encode, ones_encode;
    vec soft_demod_result;
    ber.set_size(1);
    bvec uncoded_bits_near_user, coded_bits_near_user, decoded_bits_near , decoded_bits_far_at_near;
    bvec uncoded_bits_far_user,coded_bits_far_user,decoded_bits_far;
    cvec ofdm_modulated_symbols_near, rec_symbols_far, rec_symbols_near, ofdm_modulated_symbols_far, Encode_far_again ;
    cvec convolved_symbols_near, convolved_symbols_far, Encode_far_again_filtered;
    
    //Setting the simulation parameters
    double N0 = 1.0 ;
    RNG_randomize();
    channel_near.set_noise(N0 / 2.0); //Set the noise value of the AWGN channel.
    channel_far.set_noise(N0); 
    ber.clear();
    //Clear the counter.
    berc.clear();

    //Generate the data.
    uncoded_bits_near_user =  randb(Nobits);         //The uncoded bits.
    uncoded_bits_far_user =  randb(Nobits);    
    
    ofdm_modulated_symbols_near =  Encode_users(nMod_scheme_near_user, nConvo_code_near_user, pmin_near, alpha, uncoded_bits_near_user);
    ofdm_modulated_symbols_far  =   Encode_users(nMod_scheme_far_user, nConvo_code_far_user, pmin_far, 1-alpha, uncoded_bits_far_user);
     
    convolved_symbols_near.set_size(0);
    convolved_symbols_far.set_size(0);
    
    //Conclutional operation with slow fading channel
    channel_number = 0;
    for(j=0;j<ofdm_modulated_symbols_near.length();j = j+coherence_time+FFTSIZE+CP*(conerence_nrof_symbol+1))
        convolved_symbols_near=concat(convolved_symbols_near, filter(ch_coeffs.get_row(channel_number++),1,ofdm_modulated_symbols_near(j,j+coherence_time+FFTSIZE+CP*(conerence_nrof_symbol+1) -1))) ;
       
        
    channel_number = 0 ;
    for(j=0;j<ofdm_modulated_symbols_far.length();j = j+coherence_time+FFTSIZE+CP*(conerence_nrof_symbol+1))
       convolved_symbols_far=concat(convolved_symbols_far, filter(ch_coeffs.get_row(channel_number++),1,ofdm_modulated_symbols_far(j,j+coherence_time+FFTSIZE+CP*(conerence_nrof_symbol+1) -1))) ;
    
          
    length_ofdm_near   =  ofdm_modulated_symbols_near.length();
    length_ofdm_far   =  ofdm_modulated_symbols_far.length();
   
    if(length_ofdm_near >  length_ofdm_far)
    {
        zeros_encode.set_size(length_ofdm_near -  length_ofdm_far);
        zeros_encode = zeros_c(length_ofdm_near -  length_ofdm_far);
        convolved_symbols_far = concat(convolved_symbols_far, zeros_encode);
    }
    else
    {
        zeros_encode.set_size(length_ofdm_far - length_ofdm_near);
        zeros_encode = zeros_c(length_ofdm_far - length_ofdm_near );
        convolved_symbols_near = concat(convolved_symbols_near, zeros_encode);
    }
     
   
    rec_symbols_near = channel_near(convolved_symbols_near + convolved_symbols_far);
    rec_symbols_far = channel_near(convolved_symbols_near + convolved_symbols_far)/sqrt(inv_dB(a_F_dB));
        
 
    if (length_ofdm_near > length_ofdm_far)
    {   
        rec_symbols_far =  rec_symbols_far.left(length_ofdm_far);
        decoded_bits_far =  Decode_users(nMod_scheme_far_user, nConvo_code_far_user, pmin_far, 1-alpha, rec_symbols_far);
        decoded_bits_far_at_near = Decode_users(nMod_scheme_far_user, nConvo_code_far_user, pmin_near, 1-alpha, rec_symbols_near.left(length_ofdm_far));
        Encode_far_again = Encode_users(nMod_scheme_far_user, nConvo_code_far_user, pmin_far, 1-alpha, decoded_bits_far_at_near) ;
        Encode_far_again_filtered.set_size(0);
        
        channel_number = 0 ;
        for(j=0;j<ofdm_modulated_symbols_far.length();j = j+coherence_time+FFTSIZE+CP*(conerence_nrof_symbol+1))
        Encode_far_again_filtered=concat(Encode_far_again_filtered, filter(ch_coeffs.get_row(channel_number++),1,Encode_far_again(j,j+coherence_time+FFTSIZE+CP*(conerence_nrof_symbol+1) -1))) ;
        
        rec_symbols_near  =    rec_symbols_near - concat( Encode_far_again_filtered, zeros_encode);
       
     }
     else
     {
         decoded_bits_far =  Decode_users(nMod_scheme_far_user, nConvo_code_far_user, pmin_far, 1-alpha, rec_symbols_far);
         decoded_bits_far_at_near = Decode_users(nMod_scheme_far_user, nConvo_code_far_user, pmin_far, 1-alpha, rec_symbols_near);
         Encode_far_again = Encode_users(nMod_scheme_far_user, nConvo_code_far_user, pmin_far, 1-alpha, decoded_bits_far_at_near) ;
         Encode_far_again_filtered.set_size(0);
         channel_number = 0 ;
         for(j=0;j<ofdm_modulated_symbols_far.length();j = j+coherence_time+FFTSIZE+CP*(conerence_nrof_symbol+1))
     	 Encode_far_again_filtered=concat(Encode_far_again_filtered, filter(ch_coeffs.get_row(channel_number++),1,Encode_far_again(j,j+coherence_time+FFTSIZE+CP*(conerence_nrof_symbol+1) -1))) ;
          
		rec_symbols_near    =    rec_symbols_near - Encode_far_again_filtered.left(length_ofdm_near);
		rec_symbols_near = rec_symbols_near.left(length_ofdm_near);
      }

       

    decoded_bits_near =  Decode_users(nMod_scheme_near_user, nConvo_code_near_user, pmin_near, alpha,  rec_symbols_near);     
    berc.count(uncoded_bits_near_user,decoded_bits_near);         //Count the errors.
    *(ber_near) =  berc.get_errorrate();
    berc.clear();                
    berc.count(uncoded_bits_far_user, decoded_bits_far);         //Count the errors.
    *(ber_far) =  berc.get_errorrate();
                    
                
                 
                
}
      
            



cvec Encode_users(int mod_number, int conv_number, double input_power, double alpha,  bvec block_input)
{
	    
    double N0 = 1.0;
    long int coherence_time = FFTSIZE *conerence_nrof_symbol;

    //punctured convolutional code parameters
    int constraint_length=7, j;
    Punctured_Convolutional_Code conv_code;
    ivec generators[nConvoCodes]={"0133 0171","0133 0171","0133 0171","0133 0171"};
    bmat puncture_matrix[nConvoCodes]={"1;1","1 1;1 0","1 1 0;1 0 1","1 1 0 1 0;1 0 1 0 1"};
    conv_code.set_generator_polynomials(generators[0], constraint_length);
    conv_code.set_puncture_matrix(puncture_matrix[conv_number]);

    //Transmission-Reception 
    cvec  trans_symbols, ofdm_modulated_symbols;
    cvec ones_encode;
    bvec coded_bits, decoded_bits,   rx_bits_near_user;
    vec Ec, Eb, soft_demod_result;

    //pilots formation
    ones_encode.set_size(FFTSIZE);
    ones_encode = ones_c(FFTSIZE);

    //Modulation schemes and OFDM 
    OFDM ofdm(FFTSIZE, CP);
    BPSK_c bpsk;
    QAM qam_modulate;

    coded_bits = conv_code.encode(block_input);
    if(!mod_number)
    {  
              bpsk.modulate_bits(coded_bits, trans_symbols); //The BPSK modulator
              ofdm_modulated_symbols.clear();
              ofdm_modulated_symbols.set_size(0, false);
              for(j=0;j<trans_symbols.length();j = j+coherence_time)
                ofdm_modulated_symbols =concat(ofdm_modulated_symbols, ofdm.modulate(concat(ones_encode, trans_symbols(j,j+coherence_time -1)))); //ofdm modulation with pilots

              return sqrt(input_power* alpha)*ofdm_modulated_symbols;
    }
    else
    {
            qam_modulate.set_M(pow2i(2*mod_number)); //QAM modulator 
            ofdm_modulated_symbols.set_size(0);
            qam_modulate.modulate_bits(coded_bits, trans_symbols);
            for(j=0;j<trans_symbols.length();j = j+coherence_time)
               ofdm_modulated_symbols =concat(ofdm_modulated_symbols, ofdm.modulate(concat(ones_encode, trans_symbols(j,j+coherence_time -1)))); 
            
           return sqrt(input_power* alpha)*ofdm_modulated_symbols;
                          
    }

}


bvec Decode_users(int mod_number, int conv_number, double input_power, double alpha,  cvec block_input)
{
    	
        
    //punctured convolutional code parameters
    int constraint_length=7, j, k=0;
    Punctured_Convolutional_Code conv_code; 
    ivec generators[nConvoCodes]={"0133 0171","0133 0171","0133 0171","0133 0171"};
    bmat puncture_matrix[nConvoCodes]={"1;1","1 1;1 0","1 1 0;1 0 1","1 1 0 1 0;1 0 1 0 1"};
    conv_code.set_generator_polynomials(generators[0], constraint_length);
    conv_code.set_puncture_matrix(puncture_matrix[conv_number]);

    //Noise variance
    double N0 = 1.0;

    // Transmission-Reception
    long int coherence_time = FFTSIZE *conerence_nrof_symbol;
    int channel_number = 0;
    OFDM ofdm(FFTSIZE, CP);
    vec Ec, Eb, soft_demod_result, build;
    bvec coded_bits, decoded_bits,   rx_bits_near_user;
    cvec  trans_symbols, ofdm_modulated_symbols, ofdm_demodulated, ofdm_demodulated_coherence, channel_fft_estimate;

    //Modulation schemes
    BPSK_c bpsk;
    QAM qam_modulate;

    soft_demod_result.clear();
    soft_demod_result.set_size(0,false);
    block_input /= sqrt(input_power * alpha) ;
    ofdm_demodulated  = ofdm.demodulate(block_input);
                        
    if(!mod_number) 
    {   
        for(k=0;k<ofdm_demodulated.length();k = k+coherence_time+FFTSIZE)
        {
            channel_fft_estimate = (ofdm_demodulated(k,k+FFTSIZE-1))/(1+ (N0/(2*input_power * alpha)));
            ofdm_demodulated_coherence = ofdm_demodulated(k+FFTSIZE,k+coherence_time+FFTSIZE-1);
            
            for(j=0;j<coherence_time;j+=FFTSIZE)
            {
                 build=bpsk.demodulate_soft_bits(ofdm_demodulated_coherence(j,j+FFTSIZE-1),channel_fft_estimate,N0/(2*input_power * alpha));              
                 soft_demod_result=concat(soft_demod_result,build);
            }
         } 
    }
    else               
    { 
        qam_modulate.set_M(pow2i(2*mod_number));         
        
        for(k=0;k<ofdm_demodulated.length();k = k+coherence_time+FFTSIZE)
        {
            channel_fft_estimate = (ofdm_demodulated(k,k+FFTSIZE-1))/(1+ (N0/(2*input_power * alpha)));
            ofdm_demodulated_coherence = ofdm_demodulated(k+FFTSIZE,k+coherence_time+FFTSIZE-1);
                 
            for(j=0;j<coherence_time;j+=FFTSIZE)
            {
                    build=qam_modulate.demodulate_soft_bits(ofdm_demodulated_coherence(j,j+FFTSIZE-1),channel_fft_estimate,N0/(2*input_power * alpha));
                    soft_demod_result=concat(soft_demod_result,build);
             }
        } 
    }
        
   decoded_bits = conv_code.decode(soft_demod_result);
   return decoded_bits;

}




