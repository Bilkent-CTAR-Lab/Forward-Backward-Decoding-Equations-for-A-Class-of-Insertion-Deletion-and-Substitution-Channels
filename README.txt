This simulation package contains MATLAB codes that simulate watermark coding, standard marker coding, 
half-marker coding and concatenation of LDPC-marker coding, for transmission over channels with
insertion-deletion-substitution (IDS) noise.

Definition of main variables follows the notes referenced below, which are also included in this 
package (in PDF format):

[1] J.Haghighat, "Forward-Backward Decoding Equations for A Class of Insertion Deletion 
and Substitution Channels"

half-marker codes are proposed in the following work:

[2] J. Haghighat, and T. M. Duman, "Half-Marker Codes for Deletion Channels with Applications 
in DNA Storage", Submitted to IEEE Communications letters

The codes included in this package are:

-FB_sim_half_marker.m : 
This code gives the achievable rate, as defined in [2], for standard marker codes and half-marker codes

-FB_sim_watermark.m :
This code simulates the performance of watermark codes over the IDS channels

-LDPC_and_FB_iterative.m :
This code simulates a concatenated error correcting scheme where the outer code is an LDPC code
and the inner code is a standard marker code or a half-marker code.

In order to gain a better understanding of these simulations, it is recommended to study [1] which
is included in this package

For further questions, please contact Javad Haghighat at:
javad.haghighat@tedu.edu.tr
