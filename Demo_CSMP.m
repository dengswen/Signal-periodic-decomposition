%auther:	dengswen (dengswen@gmail.com)
%date:       2018/3/13
%function:  Demo for paper in [1]

% [1] Shiwen Deng, Jiqing Han, ''Ramanujan subspace pursuit for signal periodic decomposition'', Mechanical Systems and Signal Processing, vol. 90, June 2017, pp. 79-96.

%### License
%-------
% The analysis work performed with the program(s) must be non-proprietary work. 
% Licensee and its contract users must be or be affiliated with an academic facility. 
% Licensee may additionally permit individuals who are students at such academic facility 
% to access and use the program(s). Such students will be considered contract users of licensee. 
% The program(s) may not be used for commercial competitive analysis 
% (such as benchmarking) or for any commercial activity, including consulting.

function Demo_CSMP
    clc;clear;
   
    Color1=[0.929411764705882 0.694117647058824 0.125490196078431];
    Color2=[0 0.447058826684952 0.74117648601532];
    Color3=[0.39215686917305 0.474509805440903 0.635294139385223];
    Color_Selected=Color3;    
    
    addpath('..\');
 

    [x] = GenerateForpaper(2)';

    Qmax = 50; 
    nMP = 20;
   

     [Xqs, Eqs, Err, Alpha]  = CSMP(x,Qmax, nMP);

    xr=sum(Xqs,2);
    err = norm(x-xr)^2/norm(x)^2
 
    figure;
    subplot(211)
    plot(x,'linewidth',1);
    hold on;
    plot(xr,'r');axis tight
    subplot(212)
    stem(Eqs,'Marker','none','LineWidth',1,...
     'Color', Color3);  
    title('Periodic spectral');
  

    rmpath('..\');


end

function x = GeneratePeriodSignal(M,N)
    n=0:1:N-1;
    x=cos(2*pi*n./M);    
end

function x=GenerateRampSignal
    K0=180*24;
    s6=[1:6];M6=K0/length(s6);
    s36=[1:36];M36=K0/length(s36);
    s45=[1:45];M45=K0/length(s45);
    
    x6=repmat(s6,1,M6);
    x36=repmat(s36,1,M36);
    x45=repmat(s45,1,M45);
    x=x6+x36+x45;
end

function x=GenerateBinSignal
    s=zeros(1,21);
    s(1,21)=1;
    M = 50;
    x = repmat(s,1,M*length(s));
    x = x';
end





