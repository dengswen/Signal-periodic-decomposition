% Generate several period signal for evaluating several algs.
% ref:
% [1] 2003, Orthogonal, Exactly Periodic Subspace Decomposition.pdf
% [2] 2012, A comparison of periodicity profile methods for sequence analysis.pdf
% [3] 2014, The farey-dictionary for sparse representation of periodic signals.pdf
% [4] 2010, Sparse signal decomposition for periodic signal mixtures.pdf
% [5] 2015, Intrinsic Integer-Periodic Functions for Discrete Periodicity Detection
function [x,X] = GenerateForpaper(nType)
    switch nType
      
        case 2  % used for my paper
%             N = 4320;
            N = 3060;
            x6=gen_cos_signal(17, N);
            x36=gen_cos_signal(36, N);
            x45=gen_cos_signal(45, N);
            x=x6+x36+x45;
            x=x-mean(x);            
            x=x(1:N);            
            X(1,:)=x6;
            X(2,:)=x36;
            X(3,:)=x45;
        
            
    end
    
    
end

function  x = gen_cos_signal(P, N)
    %generate period signal with period P
    n=0:1:N-1;
    x=cos(2*pi*n./P);  
end

function [x]= Generate_signal(FS,f0, N)

    %采样频率(HZ)
    sample_F=FS;                 
    fs = sample_F;
    %采样周期
    sample_T=1/sample_F;

    %采样点个数
%     N=600;

    %采样点标号
    n=0:(N-1);           

    %每个采样点的时间轴坐标
    t=n*sample_T;  


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %2. 构造信号
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %信号频率(HZ)
    w0=2*pi*f0;

    %构造含有两个频率成分的正弦信号：50HZ 100HZ

    %x=sin(2*pi*f0*t)+sin(2*pi*f1*t);
    x=sin(w0*t);
   
end

function s = Signal_amp_varying(N,lab)
	n=1:N;
	if lab==1
		s=exp(-(2*n/N).^2);
	else 
		s=exp(-(2*(n-N)/N).^2);
	end
end


%%  -----------------------------------
% Writen by Shi-Wen Deng, 2015/3/15
% function:
    % return matrix B whose first Phi(q) columns are the integeral number
    % basis of the Ramanujan space. and the Euler number Phi(q)  is equal
    % the rank of B
function [B, C] = Generate_Sq(q)
    % q: input number
    % C: the first Phi(q) columns of B which is an integer basis of the
    % Ramanujan subspace
    % B: matrix of RSq
    if q==1
        C=1;
        B=1;
        return;
    end
    
    Rs= Generate_RS(q)';
    B=zeros(q, q);
    B(:,1) = Rs;
    for i=2:q
        B(:, i) = circshift(Rs, [i-1,0]);
    end
    R = rank(B);
    C=B(:, 1:R);    
end

function Rs = Generate_RS(q)
    % Return RS of q
    %% 1) factor q base on the fundamental theorem of arithmetic
    f = factor(q);
    L = length(f);
    s = sparse(q,1,0);
    for i =1:length(f)
        p = f(i);
        s(p) = s(p)+1;
    end
    
    %% 2 generate RS for prime p
    Rs = ones(1, q);
    ps = find(s>0);
    for i=1:length(ps)
        p=ps(i);
        m=s(p);
        c = Generate_prime_power_RS(p, m);
        n = q/(p^m);
        ce = repmat(c, 1, n);
        Rs =Rs.*ce;
    end
    
end

function Rs = Generate_prime_power_RS(p, m)
    % generate Rs of p^m
    % p is a prime
    % m is q's power
    % L is totoal number periods of the return RS
    if ~isprime(p)
        Rs = -1;
        disp('Input p is not a prime!');
        return;
    end
    for i = 1:p^m
        n = i-1;
        if mod(n, p^(m-1))
            c = 0;
        else
            if mod(n, p^m)
                c = -p^(m-1);
            else
                c = p^(m-1)*(p-1);
            end
        end
        Rs(i) = c;
    end
end
function P = GenerateOrthRSP(q,N)
    Pq = Generate_Sq(q)/q;    
    M = N/q;
    P=repmat(Pq,M,M);
    P=P/M;
end
%% ---------------------------------------------


function [x, Eqs,  Qs]= GenerateExactlyMulPeriodicSig(Qmax, Qrange, K, N)
    
    Qs=[];
    x=0;
    Eqs=sparse(1,Qmax);
    for k=1:K
        randn('seed',k);
        q=Getq(Qs, Qrange);
        Qs=[Qs, q];
        
        xq = randn(1,q);
        Pq = Generate_Sq(q)/q;
        xq = xq*Pq;
        M=ceil(N/q);
        xq=repmat(xq,1,M);
        xq=xq(1:N);
        
        Eqs(1, q) = norm(xq)^2;
        x=x+xq;
    end
end

function q=Getq(Qs, Qrange)
    for i=1:max(Qrange)
        q=randi(Qrange);
        p=find(Qs==q);
        if isempty(p)
            break;
        end
        
    end
    
    
end

function x = gen_square_sig(N)
%     P = 1/2048;
    P = 1/4096;
    t = 0:P:P*N;
    t = 2*pi*t;
    x=square(2*pi*30*t,50);
    x = x(1:N);
end
function x = gen_chirp(N)
    t = 0:0.01:N*0.01;            % 2 secs @ 1kHz sample rate
    x = chirp(t,0,1,150);
    x=x(1:N);
end