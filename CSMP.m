% Writen by Shi-Wen Deng, 2015/10/4
% Update: 2015/12/28
% function: implement 
% 

function [Xqs, Eqs, Err, Alpha] = CSMP(x, Qmax, nMP)
% Input paras:
%   x: input signal 
%   Qmax: max period that for scaning
%   nMP: order of the match pursuit
% Output paras:
    %Alpha: complex coefficient vector of length M
    
    r = x;
    N = length(x);
    qs=1:Qmax;
    W = (N+qs)./(2*qs);  % weight defined in paper, from energy to periodicity meature
%     W=1./(qs.^2);
    
    Xqs = sparse(N, Qmax); % periodic components
    Eqs = 0;
    Err = 0;
    
    %% coefficients vector setting
    Mqs = 0; % vector, the numbers of complex exp basis in each Sq
    for q = 1:Qmax
        Mq = calMq(q);
        Mqs(q) = Mq;
    end
    M = sum(Mqs);
    Alpha = zeros(M,1);    %complex coefficient vector of lenth M
    
    
    
    %% construct Ramanujan subspace
%     tic
    %% perform MP alg
    EX = norm(x)^2;
    iQ = 1;
    V=[];
    P = zeros(N,N);
    for n = 1:nMP
        %% detect and extract the q-periodic subsequence
        EPE = ESPD_normal(r, Qmax);
        WEPE = W.*EPE/(N);
        [v, q] = max(WEPE);
      
        %% generate Sq and project x into it with OMP
        [Vq, c] = Generate_ComplexExp(q,N);
        
       %% Select dominant frequency bin
        vr = Vq'*r;
        alpha = 1./(1-abs(c).^2).*(vr-c.*conj(vr));
        [v, k] = max(abs(alpha));
        if q>2
            xk = 2*real(alpha(k)*Vq(:,k));
        else
            xk = real(alpha(k)*Vq(:,k));
        end
        Xqs(:, q) = Xqs(:, q) + xk;
        m = q_i2n(q, k, Mqs);
        Alpha(m) = Alpha(m) + alpha(k);
        Err(n) = norm(r)^2/EX;
        r = r - xk;
    end
    Err(n)=norm(r)^2/EX;
    %% calculate energy of subsequence
    U=Xqs.^2;
    E=sum(U,1);
    Eqs=E;
    
% 	toc
end


function ExactlyPeriodicEnergy = ESPD_normal(x, Pmax)
% written by Dengswen 2015/5/28  
% implementation of the alg in 2003, Orthogonal, Exactly Periodic Subspace Decomposition
    K0 = length(x);
%% two corr ---------
%1)
    Phi = xcorr(x);
    [v, pos] = max(Phi);
    Phi = Phi(pos:end);
%2)    has problem
%     Phi=autocorr(x);
%% ------------------    
    ExactlyPeriodicEnergy(1)=g_fun(1,K0,Phi);
    
    for P=2:Pmax
        facts=AllFactors(P);
        g_P = g_fun(P, K0,Phi);   % current energy in space P
        g = 0;
        for i=1:length(facts)
            fact=facts(i);
            g = g+ExactlyPeriodicEnergy(fact);
        end
        g_P = g_P - g;
%         g_P=g_P/P;    % this is a modification for the original alg, to elimate the effect of P
%         g_P=g_P/calEuler(P);  

        ExactlyPeriodicEnergy(P)=g_P;
    end
    
    
end

function g = g_fun(P, K0, Phi)
    M=floor(K0/P);
    r = 0;
    for i=1:(M-1)
        r=r+Phi(i*P+1);
    end
    %1)
    g = Phi(1)+2*r;  % Eq.(3)
    %2)
%     g = 2*r;    % Eq. (4)

    g= g*P/K0;
end

function fact=AllFactors(P)
    j=1;
    for i=1:P-1
        if mod(P,i)==0
            fact(1,j)=i;
            j=j+1;
        end
    end
end


% Author: Deng Shiwen
%2015/10/12
% function: generate the Ramanujan subspace from complex exp basis
function [Vq, c] = Generate_ComplexExp(q, N) 
% Parameters:
% Input:
%   q: period used to generate the Ramanujan subspace Sq
%   N, signal length
% Output:
%   Vq: half of the complex exp basis of Sq: 1, ..., floor(q)/2
%   c: correlation coefficients of <g, g*>

    n=1:N;
    if q==1 || q==2
            v=exp(-2*pi*j*n.*1./q)';
            v = v/norm(v);
            Vq=v;
            c = 0;
            return;
    end
    M = floor(q/2); 
    i=1; 
    for k=1:M
        if gcd(k,q)==1
            v=exp(-2*pi*j*n.*k./q)';
            v = v/norm(v);
            Vq(:, i)=v;
            i=i+1;
        end
    end
    C=Vq'*conj(Vq);
    c =diag(C);
end

%% calculate the number of complex exp basis, Mq, in each Ramanujan subspace Sq
function Mq = calMq(q)
    if q <=2
        Mq = 1;
    else
        Mq = round(calEuler(q)/2);
    end
end

%% transform (q,i) to the 1-dim form
function n = q_i2n(q,i, Mqs)
    Mc = 0;
    if q==1
        Mc = 0;
    else
        Mc = sum(Mqs(1:q-1));
    end
    n = Mc + i; 
end

