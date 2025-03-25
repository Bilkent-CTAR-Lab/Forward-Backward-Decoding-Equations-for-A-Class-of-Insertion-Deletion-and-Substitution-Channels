%This code simulates encoding and decoding of watermark codes 
% over insertion-deletion-and-substitution (IDS) channels.

%The main variables are named according to the following notes,
%which are accompanying these simulation codes:
%[1] J.Haghighat, "Forward-Backward Decoding Equations for A Class of Insertion Deletion 
% and Substitution Channels"

%watermark codes are given in the following paper:
%[2] M. C. Davey and D. J. Mackay, "Reliable communication over channels with 
% insertions, deletions and substitutions,” IEEE Trans. Inf. Theory, 
% vol. 47, no. 2, pp. 687–698, Feb. 2001.

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

%insertion, deletion, and substitution probabilities
pi = 0.0040;
pd = 0.0040;
ps = 0.0080;

%the maximum number of inserted symbols, see [1] for details
l_max = 2;

if(pi == 0)
    l_max = 0;
end

pil_vec = transpose(pi.^(0:l_max));

%joint insertion-deletion probability, as defined in [1]
mu = [pd*pil_vec, (1-pd)*pil_vec];

mu = mu/sum(sum(mu));

%(n,k) watermark code parameters, as defined in [2] 
k_wat = 6;
n_wat = 8;
%Finding the sparsification lookup table, as proposed in [2]
%v_orig_wat lists all 2^k_wat binary vectors of length k_wat
%v_map_wat list the 2^k_wat binary vectors of length k_wat with smallest
%possible Hamming weights
[v_orig_wat, v_map_wat, p_wat] = watermark_mapping(k_wat,n_wat);

%Length of the original information sequence, in bits
data_length = 360;

%length of the block transmitted over the IDS channel, 
%according to the notation in [1]
T = (n_wat/k_wat)*(data_length/2)

%The watermark sequence is taken as a randomly generated sequence
water_seq_bin = round(rand(1,2*T));

%rho is the a priori information vector [2]
%f is the substitution channel model ([2], Equation (4))
%zeta is the joint probability function defined in ([2], Equation (7))
[rho f zeta] = rfz_wat(p_wat,water_seq_bin,ps,T);

sent_blocks = 100;

bit_errors = 0;

symbol_errors = 0;

for sb=1:sent_blocks,

    sb = sb

    %randomly generate data bits
    ub = round(rand(1,data_length));

    %apply sparsifying, as given in [2], and using the previously found
    %lookup table (i.e., the mapping from v_orig_wat to v_map_wat)
    x = wat_marcode(ub,v_orig_wat,v_map_wat);
    
    %dividing the 4-ary marker codeword into two
    %binary sequences
    xb1 = floor(x/2);
    xb2 = rem(x,2);
    
    %adding watermark sequence to the codeword
    %the addition is defined as bit-wise X-oR
    %e.g. 11+01 = 10
    xw1 = rem(xb1+water_seq_bin(1:2:end),2);
    xw2 = rem(xb2+water_seq_bin(2:2:end),2);
    
    %re-combining the watermarked bit streams as 4-ary
    x_wat = 2*xw1+xw2;
    
    %pass through the IDS channel
    y = ids_channel(x_wat,T,pi,pd,ps,l_max);

    %decode using the forward-backward (FB) decoder
    %the output, joint_prob, is a matrix where joint_prob(j) 
    % denotes the a posteriori probability estimated by the FB decoder
    % given the j-th transmitted 4-ary symbol
    joint_prob = FB_decode_wat(y,T,mu,rho,f,zeta,log_map_vec,delta_step,l_max);

    %Now we make a hard decision by picking the value (between 0 and 3)
    %for which the a posteriori probability is maximized
    [p_x_dec, loc_x_dec] = max(transpose(joint_prob));
    x_dec = loc_x_dec - 1;
    
    %updating the count of symbol errors
    %note that we are directly comparing the hard decision, x_dec
    %with the watermarked transmitted sequence, x_wat
    %alternatively, one may subtract the watermark sequence from
    %x_dex (using bit-wise X-oR), and then compare it with the
    %sequence x
    symbol_errors = symbol_errors + sum(x_dec ~= x_wat);

end

%printing the symbol error rate
SER = symbol_errors / sent_blocks / T


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function calculates the a priori probabilities (rho)
%channel transition probabilities (f)
%and the joint probability distribution function, zeta 
%zeta is defined in ([1]-Equation (7))
%The variables are named according to the notation given in [1]
function [rho f zeta] = rfz_wat(p_wat,water_seq_bin,ps,T),

rho = [];

for k=1:T,
    
    w = water_seq_bin(2*k-1:2*k);

    if(w(1)==0)
        rb1 = p_wat;
    else
        rb1 = 1-p_wat;
    end

    if(w(2)==0)
        rb2 = p_wat;
    else
        rb2 = 1-p_wat;
    end

    rhok = [(1-rb1)*(1-rb2), (1-rb1)*rb2, rb1*(1-rb2), rb1*rb2];

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
%This function implements the watermark coding as explained in [2]
%for this, the previously found sparsification lookup table (i.e.,
%the mapping from v_orig_wat to v_map_wat) is used
%the watermark sequence is added to x after running this function and
%in the main body
function x = wat_marcode(ub,v_orig_wat,v_map_wat)

v_map_qarry = 2*v_map_wat(:,1:2:end) + v_map_wat(:,2:2:end);

k = size(v_orig_wat,2);

l = ceil(length(ub)/k);

x = [];

for i=0:l-1,

    ub_seg = ub(i*k+1:(i+1)*k);

    for j=1:size(v_orig_wat,1),

        if(isequal(ub_seg,v_orig_wat(j,:))),

            x = [x,v_map_qarry(j,:)];

            break;
        end
    end

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function implements the watermark decoding as explained in [2]
%for this, the previously found sparsification lookup table (i.e.,
%the mapping from v_orig_wat to v_map_wat) is used
function x_dec = wat_mardec(p_ub_1,v_orig_wat,v_map_wat),

v_map_qarry = 2*v_map_wat(:,1:2:end) + v_map_wat(:,2:2:end);

n = size(v_map_wat,2);

l = ceil(length(p_ub_1)/n);

x_dec = [];

for i=0:l-1,

    p_seg = p_ub_1(i*n+1:(i+1)*n);

    [u_dec, dum, j_index] = watermark_demapping(v_orig_wat,v_map_wat,p_seg);

    x_dec = [x_dec, v_map_qarry(j_index,:)];

end

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
%the output, joint_prob, is a matrix where joint_prob(j) 
% denotes the a posteriori probability estimated by the FB decoder
% given the j-th transmitted 4-ary symbol
%T is the block length, mu is the joint insertion-deletion probability
%defined in [1]

function joint_prob = FB_decode_wat(y,T,mu,rho,f,zeta,log_map_vec,delta_step,l_max),

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



bit_llrs = [log(bit_prob1(1+1,:)./bit_prob1(0+1,:)), log(bit_prob2(1+1,:)./bit_prob2(0+1,:))];

p_ub_1 = exp(bit_llrs) ./ (1+exp(bit_llrs));


sum_jp = sum(transpose(joint_prob));

s_mat = transpose(sum_jp)*ones(1,4);

joint_prob = joint_prob ./ s_mat;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the mapping to be used in watermark coding
%the mapping is defined according to [2]
function [vo vw p] = watermark_mapping(k,n),

i_vec = 0:2^n-1;

vwn = [];

for j=1:n,

    vwn = [vwn; rem(i_vec,2)];

    i_vec = floor(i_vec/2);

end

i_vec = 0:2^k-1;

vo = [];

for j=1:k,

    vo = [vo; rem(i_vec,2)];

    i_vec = floor(i_vec/2);

end

[iw, jw] = sort(sum(vwn));

vw = transpose(vwn(:,jw(1:2^k)));

vo = transpose(vo);

p = sum(sum(vw))/prod(size(vw));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the demapping to be used in watermark decoding
%the mapping is defined according to [2]
function [y x j] = watermark_demapping(vo,vw,p1_vec),

l = size(vw,1);

p1_vec_ext = ones(l,1)*p1_vec;

h = vw.*p1_vec_ext + (1-vw).*(1-p1_vec_ext);

[i, j] = max(sum(transpose(h)));

x = vo(j,:);

y = vw(j,:);

end