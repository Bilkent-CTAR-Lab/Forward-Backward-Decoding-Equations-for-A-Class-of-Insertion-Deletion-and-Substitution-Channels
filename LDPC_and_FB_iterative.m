%This code simulates encoding and decoding of a concatenated
% coding scheme, where the outer code is a regular LDPC code
% and the inner code is either a standard marker or a half-marker code
% The channel is an insertion-deletion-and-substitution (IDS) channel.

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

% the following 4-ary mapping is used: [0,1,2,3] = [00,01,10,11]

%block length, according to the notation in [1]
T = 180;
%T = 1000;

%insertion, deletion, and substitution probabilities
pi = 0.000;
pd = 0.030;
ps = 0.020;

%the maximum number of inserted symbols, see [1] for details
l_max = 1;

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

mp1 = -1*ones(1,T);
mp2 = -1*ones(1,T);

%the period by which markers are inserted, according to [2]
Np = 6;

mp1(1:Np:end) = 1;
mp1(2:Np:end) = 0;
%mp1(3:Np:end) = 1;
%mp1(4:Np:end) = 0;

%mp2(1:Np:end) = 0;
%mp2(2:Np:end) = 0;
%mp2(1+Np:2*Np:end) = 1;
% mp2(2+Np:2*Np:end) = 0;

%rho is the a priori information vector [2]
%f is the substitution channel model ([2], Equation (4))
%zeta is the joint probability function defined in ([2], Equation (7))
%Note that since we allow exchange of information between the LDPC decoder
%and the forward-backward (FB) decoder of marker codes,
%the a priori information vector, rho, will be updated at each 
% information exchange iteration
%since zeta is defined based on rho, zeta will be updated as well
[rho f zeta] = rfz(mp1,mp2,ps,T);

rho_orig = rho;

zeta_orig = zeta;

%the total number of non-marker bits in a
%codeword of the marker code
%note that these bits correspond to a codeword of the
%outer code (the LDPC code)
data_length = sum(mp1 == -1) + sum(mp2 == -1);

%marker code rate
rM = (data_length/(2*T));

%variable and check node degrees of the LDPC code
dv = 3;
dc = 6;

%The number of information bits input to the LDPC encoder
message_length = 150;

%Constructing the parity-check matrix of the LDPC code
%and its corresponding generator matrix
[H_orig, G_sys] = build_incomplete_Gallager(dv,dc,data_length,message_length);

sent_blocks = 10000;

bit_errors = 0;

message_bit_errors = 0;

correct_frames = 0;

symbol_errors = 0;

for sb=1:sent_blocks,

    if(rem(sb,100)==0),
        sb_e_c=[sb, message_bit_errors, correct_frames]
    end

    %input('press a key');
    
    %randomly generate the information bits to be encoded by
    %the LDPC code
    u_message = round(rand(1,message_length));
    
    %generate an LDPC codeword
    ub = rem(u_message*G_sys, 2);

     %apply marker coding
    x = marcode(ub,mp1,mp2);

    %pass through the IDS channel
    y = ids_channel(x,T,pi,pd,ps,l_max);
    
    %This is the number of times for which the information is exchanged
    %between the forward-backward (FB) decoder and the LDPC decodrer
    %Example: max_iters = 1 implies that FB decoding is followed by LDPC decoding
    %and then making a decision
    %Example: max_iter = 2 means that after the first round of FB and LDPC decoding, 
    % LLR values obtained by the LDPC decoder are passed back to the FB 
    %decoder and then the LLRs updated by the FB decoder are passed to
    %the LDPC decoder and a decision is made after this LDPC decoding
    max_iters = 5;

    %During the information exchange process, the vectors
    %rho and zeta which are employed by the FB decoder, will be updated
    %at each exchange
    %we begin with their original values found above
    rho = rho_orig;
    zeta = zeta_orig;

    iter = 1;

    while (iter <= max_iters),

        %decode using the forward-backward (FB) decoder
        %the output, p_ub_1, is a vector where p_ub_1(j) 
        % denotes the a posteriori probability estimated by the FB decoder
        % given the j-th transmitted bit is 1
        p_ub_1 = FB_decode(y,T,mu,rho,f,zeta,mp1,mp2,log_map_vec,delta_step,l_max);
        
        %The LLR values found by FB decoding
        llr_bench = log((1-p_ub_1)./(p_ub_1));
        
        %Find the LLRs with absolute values greater than 100
        %and truncate them to avoid overflows 
        [pi0, pj0] = find(p_ub_1 == 0);
        [pi1, pj1] = find(p_ub_1 == 1);
        llr_bench(pj0) = 100;
        llr_bench(pj1) = -100;

        %Total number of decoding iterations
        %this is the product of the number of times
        %that information is exchanged between the FB and the LDPC
        %decoders, multipled by the number of LDPC decoding iterations
        %following each FB decoding iteration
        total_iters = 20;
        
        %This is the number of LDPC decoding iteration following
        %each FB decoding
        ldpc_iters = total_iters/max_iters;

        %Running LDPC decoding and finding the LLRs corresponding
        %to the codeword bits (ldpc_llrs), 
        % a tentative hard decision on message bits (u_message_dec)
        %and the syndrome corresponding to the hard-decision (synd)
        [ldpc_llrs,u_message_dec,synd] = ldpc_dec_soft(llr_bench, H_orig, ldpc_iters);
        
        
        p1_ldpc_vec = 1./(1+exp(ldpc_llrs));

        %dating rho and zeta vectors to be employed by the 
        %FB decoder at next decoding iteration
        [rho zeta] = update_rho_zeta(mp1,mp2,ps,T,p1_ldpc_vec);

        sum_synd = sum(synd);
        
        %if the syndrome of the hard-decision vector
        %is is an all-zero vector, we accept this codeword
        %as the final decoding result and do not iterate
        %the LLRs back to the FB decoder
        if(sum_synd == 0),
            iter = max_iters+1;

        else 
            iter = iter + 1;
        end

        
    end

    %Updating the count on bit errors and frame errors
    message_bit_errors = message_bit_errors + sum(u_message_dec ~= u_message);
    correct_frames = correct_frames + isequal(u_message_dec,u_message);


end



%We only print the Bit Error Rate
%Since the frame errors are also counted, one may also
%print the Frame Error Rate
message_BER = message_bit_errors / sent_blocks / message_length

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

bit_llrs1 = [log(bit_prob1(1+1,j1)./bit_prob1(0+1,j1))];

[pi0, pj0] = find(bit_prob1(1+1,j1) == 0);

bit_llrs1(pj0) = -100;

[pi0, pj0] = find(bit_prob1(0+1,j1) == 0);

bit_llrs1(pj0) = 100;

bit_llrs2 = [log(bit_prob2(1+1,j2)./bit_prob2(0+1,j2))];

[pi0, pj0] = find(bit_prob2(1+1,j2) == 0);

bit_llrs2(pj0) = -100;

[pi0, pj0] = find(bit_prob2(0+1,j2) == 0);

bit_llrs2(pj0) = 100;

bit_llrs = [bit_llrs1, bit_llrs2];

p_ub_1 = exp(bit_llrs) ./ (1+exp(bit_llrs));


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho, zeta] = update_rho_zeta(mp1,mp2,ps,T,p1_vec),

[i1, j1] = find(mp1 == -1);
[i2, j2] = find(mp2 == -1);

prob_mat = [mp1;mp2];

prob_mat(1,j1) = p1_vec(1:length(j1));
prob_mat(2,j2) = p1_vec(length(j1)+1:length(j1)+length(j2));

pb1 = prob_mat(1,:);
pb2 = prob_mat(2,:);

rho = [(1-pb1).*(1-pb2); (1-pb1).*pb2; pb1.*(1-pb2); pb1.*pb2];

rho =rho ./ (ones(4,1)*sum(rho));

rho = transpose(rho);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function performs LDPC decoding using the sum-product update rules
%final_LLRs are the LLRs calculated for each codeword bit
%u_dec is a tentative hard decision 
%synd is the syndrome corresponding to the tentative hard decision
%Li is the vector of a priori information received through the
%inner decoder (FB decoder)
%H_orig is the sparse parity-check matrix of the LDPC code
function [final_LLRs,u_dec,synd] = ldpc_dec_soft(Li, H_orig, Max_iter),

%If LLR values are larger than this threshold (in absolute value), 
% they are truncated to this absolute value to avoid overflows
Large_LLR = 1e2;

[m, n] = size(H_orig);

mv = H_orig .* (ones(m,1)*Li);

mvc = mv;



for t = 1:Max_iter,

    %Check node to variable node message passing

    mcv = zeros(size(H_orig));

    for cn = 1:m,

        [dum_vec, j_vec] = find(H_orig(cn,:));

        for l = 1:length(j_vec),

            loc = j_vec(l);

            prod1 = 1;

            for lt = 1:length(j_vec),

                if lt ~= l,

                    loct = j_vec(lt);

                    prod1 = prod1 * tanh(mvc(cn,loct)/2);

                 end

            end

            if abs(prod1) ~= 1,
                
                mcv(cn,loc) = log((1+prod1)/(1-prod1));
            else

                mcv(cn,loc) = prod1 * Large_LLR;
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Variable node to check node message passing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mvc = zeros(m,n);

    for vn = 1:n,

        sv = sum(mcv(:,vn));

        [i_vec, dum_vec] = find(H_orig(:,vn));

        for l=1:length(i_vec),

            loc = i_vec(l);

            mvc(loc,vn) = sv - mcv(loc,vn);

        end

    end

    mvc = mvc + mv;
    %%%%

    temp_LLRs = sum(mcv) + Li;

    x_hat_temp = round(temp_LLRs < 0);

    synd_temp = mod(H_orig*transpose(x_hat_temp),2);

    if(sum(synd_temp) == 0),
        %check = 1
        t = Max_iter + 1;
    end


end

final_LLRs = sum(mcv) + Li;

x_hat = round(final_LLRs < 0);

synd = mod(H_orig*transpose(x_hat),2);

u_dec = x_hat(m+1:n);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Constructing a sparse parity check matrix for a regular LDPC 
%code based on Gallager's method
%and finding its corresponding systematic generator matrix
%using Gaussian elimination
%dv, dc, n, ki are variable node degree, check node degree
%codeword length and message length
function [H_orig, G_sys] = build_incomplete_Gallager(dv,dc,n,ki),

l = n/dc;

inv_flags = [];

force_flip = 0;

A = zeros(l, dc*l);

for i = 1:l,
    A(i,dc*i-dc+1:dc*i) = ones(1,dc);
end
%%%%%
H_orig = A;

for j=2:dv,

    H_orig = [H_orig; A(:,randperm(dc*l))];

end

H_orig = H_orig(1:n-ki,:);

%%%%%
[m, n] = size(H_orig);

H_sys = H_orig;

inv_flag = 1;

rni = 0;

for i = 1:m,

    inv_flag_i = 1;

    if(H_sys(i,i) == 0),

        inv_flag_i = 0;

        for t = i+1:m,

            if(H_sys(t,i) == 1),

                inv_flag_i = 1;

                temp_vec = H_sys(i,:);
                H_sys(i,:) = H_sys(t,:);
                H_sys(t,:) = temp_vec;

                temp_vec = H_orig(i,:);
                H_orig(i,:) = H_orig(t,:);
                H_orig(t,:) = temp_vec;
            end

        end

        
    end

    if inv_flag_i == 0,

        force_flip = force_flip + 1;
        
        H_sys(i,i) = 1;

       
        H_orig(i,i) = mod(H_orig(i,i)+1,2);

        
        inv_flag_i = 1;

    end

    if inv_flag_i ~= 0,
    
        for j = 1:size(H_sys,1),

            if j ~= i,
    
                if H_sys(j,i) == 1
                
                    V = H_sys(j,:) + H_sys(i,:);
    
                    H_sys(j,:) = V - 2*floor(V/2);
    
                end
    
            end

        end
    end

    inv_flags = [inv_flags, inv_flag_i];

    rni = rni + inv_flag_i;

    ir = [i, rni];

    inv_flag = inv_flag * inv_flag_i;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
end


mf = size(H_sys,1);

G_sys = [transpose(H_sys(:,mf+1:n)), eye(n-mf)];

end

