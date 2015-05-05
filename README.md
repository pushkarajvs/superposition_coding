#Superposition coding

- This is the implementation of superposition coding scheme to expand the rate-region in a broadcast wireless communication scenario.

- To use this implementation, 
    1) Install ITPP library from [here](http://itpp.sourceforge.net/4.3.1/installation.html) which is a standard C++ library available for communication engineers containing implementations of many standard algortihms.
    
    2) Install Octave by typing the following command. sudo apt-get install octave
    
    3) To build the project, clone it by typing git clone https://github.com/
    
    4) Go to directory superposition_coding and run build.sh by typing ./build.sh. It will generate minimum power vector for the near and the far users using single_user_12rates.cpp one by one and then locate the rate-pairs outside the region spanned by the conventional orthogonal schemes such as TDMA, CDMA etc. by running superposition_coding.cpp file. Finally, the octave script plot_rate_region.m is run to plot the results.
    
    5) The algorithm used to implement this system was referred from,
    
    S. Vanka et al., 2012, Superposition Coding Stategies: Design and Experimental Evaluation, *IEEE Transactions on Wireless Communications*, Vol. 11, No. 7, July 2012.
