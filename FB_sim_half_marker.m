%This code simulates encoding and decoding of standard marker
%and half-marker codes over insertion-deletion-and-substitution (IDS)
%channels.

%The main variables are named according to the following notes,
%which are accompanying these simulation codes:

%[1] J.Haghighat, "Forward-Backward Decoding Equations for A Class of Insertion Deletion 
% and Substitution Channels"

%Half-marker codes are proposed in the following work:

%[2] J. Haghighat, and T. M. Duman, "Half-Marker Codes for Deletion Channels 
% with Applications in DNA Storage", Submitted to IEEE Communications
% letters

%In order to support large block lengths, the decoder 
% is implemented in logarithmic domain, using the log-likelihood algebra
%For details on log-likelihood algebra, see relevant references including:

%[3] J. Hagenauer, E. Offer and L. Papke, "Iterative decoding of binary block and convolutional codes," 
% in IEEE Transactions on Information Theory, vol. 42, no. 2, pp. 429-445, March 1996.



clear;

%defining the function required to implement log-domain decoding [3]
delta_step = 0.01;
xx_vec = 0:delta_step:10;
log_map_vec = log(1+exp(-xx_vec));


%these vectors will save the LLR values corresponding to
%transmitted 0 and 1 bits, respectively
%this is required to estimate conditional entropy values
%and the mutual information between transmitted bits and their
%corresponding LLRs after decoding
llr0 = [];
llr1 = [];


% the following 4-ary mapping is used: [0,1,2,3] = [00,01,10,11]

%block length, according to the notation in [1]
T = 1000;

%insertion, deletion, and substitution probabilities
pi = 0.010;
pd = 0.010;
ps = 0.020;

%the maximum number of inserted symbols, see [1] for details
l_max = 3;

if(pi == 0)
    l_max = 0;
end


pil_vec = transpose(pi.^(0:l_max));

%joint insertion-deletion probability, as defined in [1]
mu = [pd*pil_vec, (1-pd)*pil_vec];
mu = mu/sum(sum(mu));

%We define two marker patterns, by introducing two vectors mp1 and mp2
%the vector of 4-ary symbols is equivalent to two binary vectors 
%mp1 (mp2) denotes marker bits of the first (second) binary vector
%if the j-th bit is a marker with a 0 or 1 value, mp1(j) denotes that value
%otherwise, mp1(j) is -1, meaning that the j-th bit is 
% an information bit which will be determined later

%when we employ standard marker symbols, both mp1(j) and mp2(j) are
%either marker bits or information bits

%for half-marker symbols, all marker bits are defined on mp1, see [2]

%the period by which markers are inserted, according to [2]
Np = 9;

mp1 = -1*ones(1,T);
mp2 = -1*ones(1,T);

mp1(1:Np:end) = 1;
mp1(2:Np:end) = 0;
%mp1(3:Np:end) = 1;
%mp1(4:Np:end) = 0;

%mp2(1:Np:end) = 1;
%mp2(2:Np:end) = 0;
%mp2(1+Np:2*Np:end) = 1;
% mp2(2+Np:2*Np:end) = 0;


%rho is the a priori information vector [2]
%f is the substitution channel model ([2], Equation (4))
%zeta is the joint probability function defined in ([2], Equation (7))
[rho f zeta] = rfz(mp1,mp2,ps,T);

%the total number of information bits in one block
data_length = sum(mp1 == -1) + sum(mp2 == -1);

%marker code rate
rM = (data_length/(2*T))

%The number of simulated blocks
%a larger number of simulated blocks gives a more accurate results
sent_blocks = 1000;

bit_errors = 0;

symbol_errors = 0;

for sb=1:sent_blocks,

    sb = sb

    %randomly generate data bits
    ub = round(rand(1,data_length));
    
    %apply marker coding
    x = marcode(ub,mp1,mp2);

    %pass through the IDS channel
    y = ids_channel(x,T,pi,pd,ps,l_max);

    %decode using the forward-backward (FB) decoder
    %the output, p_ub_1, is a vector where p_ub_1(j) 
    % denotes the a posteriori probability estimated by the FB decoder
    % given the j-th transmitted bit is 1
    p_ub_1 = FB_decode(y,T,mu,rho,f,zeta,mp1,mp2,log_map_vec,delta_step,l_max);

    %hard-decision on transmitted bits
    ub_dec = round(p_ub_1 >= 0.5);
    
    %encoding the hard decisions bits using the marker code
    %to find hard decisions on transmitted symbols
    x_dec = marcode(ub_dec,mp1,mp2);


    %updating the count of bit errors
    bit_errors = bit_errors + sum(ub_dec ~= ub);
    
    %updating the count of symbol errors
    symbol_errors = symbol_errors + sum(x_dec ~= x);
    
    %The LLR values found by FB decoding
    llr_bench = log10(p_ub_1 ./ (1-p_ub_1));

    
    %separatly record the LLR values corresponding 
    % to 0 and 1 transmitted bits
    [i0, j0] = find(ub == 0);
    [i1, j1] = find(ub == 1);

    llr0 = [llr0, llr_bench(j0)];
    llr1 = [llr1, llr_bench(j1)];

end

%The bit error and symbol error rate achieved by hard decision
BER = bit_errors / sent_blocks / data_length
SER = symbol_errors / sent_blocks / T

%Approximating H(LLR|u=0) using the histogram as explained in [2]
[x0, y0] = hist(llr0);
px0 = (x0+1e-10)/sum(x0);
h0 = sum(-px0 .* log2(px0))+log2(mean(diff(y0)));

%Approximating H(LLR|u=1) using the histogram as explained in [2]
[x1, y1] = hist(llr1);
px1 = (x1+1e-10)/sum(x1);
h1 = sum(-px1 .* log2(px1))+log2(mean(diff(y1)));

%Approximating H(LLR) using histogram
[x01, y01] = hist([llr0 llr1]);
px01 = (x01+1e-10)/sum(x01);
h01 = sum(-px01 .* log2(px01))+log2(mean(diff(y01)));

%Since H(LLR|u) = P(u=0)*H(LLR|u=0)+P(u=1)*H(LLR|u=1)
%the mutual information, I(LLR;u) is estimated as below
Mutual_inf = h01 - 0.5*(h0+h1)

%Achievable rate, denoted by Ra in [2]
achievable_rate = Mutual_inf * rM


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function calculates the a priori probabilities (rho)
%channel transition probabilities (f)
%and the joint probability distribution function, zeta 
%zeta is defined in ([1]-Equation (7))
%The variables are named according to the notation given in [1]
function [rho f zeta] = rfz(mp1,mp2,ps,T),

rho = [];

for k=1:T,

    if(mp1(k)==0 & mp2(k)==0),      rhok = [1 0 0 0];
    elseif(mp1(k)==0 & mp2(k)==1),  rhok = [0 1 0 0];
    elseif(mp1(k)==1 & mp2(k)==0),  rhok = [0 0 1 0];
    elseif(mp1(k)==1 & mp2(k)==1),  rhok = [0 0 0 1];
    elseif(mp1(k)==0 & mp2(k)==-1), rhok = [0.5 0.5 0 0];
    elseif(mp1(k)==1 & mp2(k)==-1), rhok = [0 0 0.5 0.5];
    elseif(mp1(k)==-1 & mp2(k)==0), rhok = [0.5 0 0.5 0];
    elseif(mp1(k)==-1 & mp2(k)==1), rhok = [0 0.5 0 0.5];
    else                            rhok = [0.25 0.25 0.25 0.25];
    end

    rho = [rho; rhok];

end

q = 4;

f = (ps/(q-1))*ones(q,q);

for i=1:q,
    f(i,i) = 1-ps;
end

qo = 1;

zeta = zeros(T,q);

for j=1:T,

    for i=0:q-1,

        for iprim = 0:q-1,

            zeta(j,i+qo) = zeta(j,i+qo) + rho(j,iprim+qo)*f(iprim+qo,i+qo);

        end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function encodes the information bits, ub, into
%a 4-ary marker codeword, according to the marker patterns
%mp1 and mp2
function x = marcode(ub,mp1,mp2)

[i1, j1] = find(mp1 == -1);

[i2, j2] = find(mp2 == -1);

xb1 = mp1;
xb2 = mp2;

xb1(j1) = ub(1:length(j1));

xb2(j2) = ub(length(j1)+1:length(j1)+length(j2));

x = xb1*2+xb2;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function simulates a realization of the IDS channel
%as formulated in [1]
%The notation is according to [1]
%y is the channel output
%x is the 4-ary marker codeword input to the channel
%T is the block length
%pi, pd, ps are insertion, deletion, and substitution probabilities
%l_max is the maximum number of inserted symbols, defined in [1]
function y = ids_channel(x,T,pi,pd,ps,l_max),

q = 4;

y = [];

for k=1:T,

    ins_num = 0;

    rpi = round(rand < pi);

    while (ins_num < l_max & rpi == 1),

        ins_sym = ceil(rand*q) - 1;

        y = [y, ins_sym];

        ins_num = ins_num + 1;

    end

    rpd = round(rand < pd);

    rps = round(rand < ps);

    if rpd == 1,

        y=y;

    elseif (rpd == 0 & rps==0),

        y=[y, x(k)];

    else
        sn = ceil(q*rand)-1;
        y = [y, rem(x(k)+sn,q)];
    end


end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function performs forward-backward decoding
%It follows the notation given in [1]
%the output, p_ub_1, is a vector denoting the probability that
%a transmitted bit is equal to 1, given the received vector, y
%T is the block length, mu is the joint insertion-deletion probability
%defined in [1]

function p_ub_1 = FB_decode(y,T,mu,rho,f,zeta,mp1,mp2,log_map_vec,delta_step,l_max),

%if an LLR value is larger than Mval, it is truncated to Mval
%This is to avoid overflows
Mval = 1e7;

R = length(y);

q = 4;

qo = 1; mo = 1; ko = 1; no = 2+l_max;


%generate alfa and beta matrices
%l_max = 0;

log_alf = -Mval*ones(T+1,R+2*l_max+3);
log_bet = -Mval*ones(T+1,R+2*l_max+3);


pd = mu(1,1);
mu00 = mu(1,1);

%initialize alfa

log_alf(0+ko,0+no) = 0;

if mu00 == 0,
    log_alf((1:T)+ko,0+no) = -Mval;
else
    log_alf((1:T)+ko,0+no) = log(mu00)*(1:T);
end

%initialize beta

log_bet(T+ko,R+no) = 0;

if mu00 > 0,
    log_bet((0:T-1)+ko,R+no) = 0;
else
    log_bet((0:T-1)+ko,R+no) = -Mval;
end


%recursive equations for alfa
for k = 1:T,
    for n = 1:R,
        for l=0:l_max,
            for b=0:1,

                coef1 = mu(l+mo,b+mo)*(q^(-l))*((zeta(k,y(n)+qo))^b);

                if (coef1 == 0),
                    calf = -Mval;
                else
                    calf = log(coef1) + log_alf(k-1+ko,n-l-b+no);
                end

                tca = abs(calf-log_alf(k+ko,n+no));

                tca_loc = min(length(log_map_vec),1+floor(tca/delta_step));

                log_alf(k+ko,n+no) = max(log_alf(k+ko,n+no),calf) + log_map_vec(tca_loc);


            end
        end
    end
end

%recursive equations for beta
for k = T-1:-1:0,
    for n = R-1:-1:0,
        for l=0:l_max,
            for b=0:1,

                zky = zeta(k+1,y(min(R,n+l+1))+qo);
                %%%%%%%%%%%%%%%%%%%
                % the term min(R,n+l+1) is for cases where the index
                %of y exceeds the boundary limit R
                %in such cases, the coefficient below is not important,
                %since bet(k+1+ko,n+l+b+no) will be zero

                coef1 = mu(l+mo,b+mo)*(q^(-l))*(zky^b);

                if(coef1 == 0),
                    cbet = -Mval;
                else
                    cbet = log(coef1)+log_bet(k+1+ko,n+l+b+no);
                end

                tcb = abs(cbet-log_bet(k+ko,n+no));

                tcb_loc = min(length(log_map_vec),1+floor(tcb/delta_step));

                log_bet(k+ko,n+no) = max(log_bet(k+ko,n+no),cbet) + log_map_vec(tcb_loc);



            end
        end
    end
end


%Run forward-backward decoding and find a posteriori probabilities


log_ap_prob = -Mval * ones(T,q);
joint_prob = zeros(T,q);


for k=1:T,
    for a=0:q-1,
        for l=0:l_max,
            for b=0:1,

                Upper_val = min(R,(l_max+1)*(k-1));

                for n=0:Upper_val,

                    fbb = (f(a+qo,y(min(R,n+l+1))+qo))^b;
                    %min(R,n+l+1) is for when the index exeeds the
                    %boundary limit, R

                    coef1 = mu(l+mo,b+mo)*(q^(-l))*fbb;

                    if(coef1 == 0),
                        log_add_term = -Mval;
                    else
                        log_add_term = log(coef1)+log_alf(k-1+ko,n+no)+log_bet(k+ko,n+l+b+no);
                    end



                    dpa = abs(log_add_term - log_ap_prob(k,a+qo));

                    dpa_loc = min(length(log_map_vec),1+floor(dpa/delta_step));

                    log_ap_prob(k,a+qo) = max(log_ap_prob(k,a+qo),log_add_term) + log_map_vec(dpa_loc);

                end

            end
        end



    end
end

mmap = max(max(log_ap_prob));

joint_prob = rho .* exp(log_ap_prob - mmap);



bit_prob1 = zeros(2,T);

bit_prob2 = zeros(2,T);

for k=1:T,

    bit_prob1(0+1,k) = joint_prob(k,0+qo)+joint_prob(k,1+qo);
    bit_prob1(1+1,k) = joint_prob(k,2+qo)+joint_prob(k,3+qo);

    bit_prob2(0+1,k) = joint_prob(k,0+qo)+joint_prob(k,2+qo);
    bit_prob2(1+1,k) = joint_prob(k,1+qo)+joint_prob(k,3+qo);

end

[i1, j1] = find(mp1 == -1);

[i2, j2] = find(mp2 == -1);

bit_llrs = [log(bit_prob1(1+1,j1)./bit_prob1(0+1,j1)), log(bit_prob2(1+1,j2)./bit_prob2(0+1,j2))];

p_ub_1 = exp(bit_llrs) ./ (1+exp(bit_llrs));


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%