function [BM_out] =carfac(sig,fs,nsec, duration, num_rep)
%CARFAC  bioplausible cochear response to the incoming audio
%   sig: audio after the required processings
%   fs:  audio sampling frequency
%   nsec: number of sections in the cochlea between apex and base;
%   determines the number of cascaded filters

%BM_hpf_out, BM_out,,BM_li_out IHC_li_out, IHC_out, 

%duration =1; 
%num_rep = 5;                % number of repetition 
dur = duration * num_rep;       % Total simulation time

sig=lowpass(sig,1000,fs);

stimulus=repmat(sig,num_rep,1);
npoints = dur*fs;           % Stimulus points' number

%% CAR

%%%%%%%%%%%% BM parameters %%%%%%%%%%%%%%%%
xlow = 0.1;                % lowest frequency position along the cochlea 
xhigh = 0.9;              % highest frequency position along the cochlea
 
x = linspace(xhigh,xlow,nsec)';% position along the cochlea 1 = base, 0 = apex
f =165.4*(10.^(2.1*x)-1); %165.4*(10.^(2.1*x)-1);      % Greenwood for humans (Charactristic frequencies)
%f= 90*(10.^(0.85*x)+9); % Greenwood for songbirds

a0 = cos(2*pi*f/fs);            % a0 and c0 control the poles and zeros
c0 = sin(2*pi*f/fs);


damping = 0.16;               % damping factor
r  = 1 - damping*2*pi*f/fs;  % pole & zero actual radius
%r0=r;
r1 = 1 - damping*2*pi*f/fs;  % pole & zero minimum radius (r = r1 + (delta_r>0))
h = c0;                      % p302 h=c0 puts the zeros 1/2 octave above poles
g = (1-2*a0.*r+r.*r)./(1-(2*a0-h.*c0).*r+r.*r);   % p303 this gives 0dB DC gain for BM
%g0=g;

%%%%%%%%%%%%% Digital IHC parameters %%%%%%%%%%%%
f_hpf = 200;                     % p328 20Hz corner for the BM HPF
q = 1/(1+(2*pi*f_hpf/fs));      % corresponding IIR coefficient 

tau_in = 10e-3;                 % p329 transmitter creation time constant
c_in = 1/(fs*tau_in);           % p329 corresponding IIR coefficient    
tau_out = 0.5e-3;               % p329 transmitter depletion time constant
c_out = 1/(fs*tau_out);         % p329 corresponding IIR coefficient 
tau_IHC = 80e-6;                % p329 ~8kHz LPF for IHC output
c_IHC = 1/(fs*tau_IHC);         % corresponding IIR coefficient 

%%%%%%%%%%% Digital OHC parameters %%%%%%%%%%%%%%%%
scale = 0.1;                    % p313 NLF parameter
offset = 0.04;                  % p313 NLF parameter
b = 1.0;                        % automatic gain loop feedback (1=no undamping).
d_rz = 0.7*(1 - r1);            % p310 relative undamping

%%%%%%%%%%%%%% AGC parameters %%%%%%%%%%%%%
tau_AGC = .002 * 4.^(0:3); % p336

% The AGC filters are decimated, i.e., running at a lower sample rate so fs
% must be reduced
Nsf = (8*2.^(0:3));             % decimation rate: 8,16,32,64 %Q: what happens if these parameters are changed?
c_AGC = Nsf./(fs*tau_AGC);

% spatial filtering
shift_AGC = 0.65 * c_AGC .* (sqrt(2).^(0:3));
spread_sq_AGC = (1.65^2 + 1) * c_AGC .* (2.^(0:3));
cr = (spread_sq_AGC + shift_AGC.^2 - shift_AGC)/2;  % sa in schaik
cl = (spread_sq_AGC + shift_AGC.^2 + shift_AGC)/2;  % sb in schaik
cc = 1 - cl - cr;   % sc in schaik

% initialise internal states
W0 = zeros(nsec,1);             % BM filter internal state
W1 = zeros(nsec,1);             % BM filter internal state
W1old = zeros(nsec,1);          % BM filter internal state at t-1 for velocity calculation
BM = zeros(nsec,npoints);       % BM displacement
BM_hpf = zeros(nsec,npoints);   % BM displacement high-pass filtered at 20Hz
trans = ones(nsec,1);           % transmitter available
IHC = zeros(nsec,npoints);      % IHC output
IHCa = zeros(nsec,npoints);     % IHC filter internal state
In8 = zeros(nsec,1);            % Accumulator for ACG4
In16 = zeros(nsec,1);           % Accumulator for AGC3
In32 = zeros(nsec,1);           % Accumulator for AGC2
In64 = zeros(nsec,1);           % Accumulator for AGC1
AGC = zeros(nsec,npoints);      % AGC filter internal state
AGC0 = zeros(nsec,1);           % AGC filter internal state
AGC1 = zeros(nsec,1);           % AGC filter internal state
AGC2 = zeros(nsec,1);           % AGC filter internal state
AGC3 = zeros(nsec,1);           % AGC filter internal state

flag=false;
%%

for t = 1:npoints
    % BASILAR MEMBRANRE
    % for the first section
    W0new = stimulus(t) + r(1)*(a0(1)*W0(1) - c0(1)*W1(1));
    if (isnan(W0new))&&~flag
        disp(t);
        flag=true;
    end
    W1old(1) = W1(1);   %=====================
    W1(1) = r(1)*(a0(1)*W1(1) + c0(1)*W0(1));
    W0(1) = W0new;
    BM(1,t) = g(1)*(stimulus(t) + h(1)*W1(1));    
    % all other sections
    for s = 2:nsec
        W0new = BM(s-1,t) + r(s)*(a0(s)*W0(s) - c0(s)*W1(s));
        W1old(s) = W1(s);    %=====================
        W1(s) = r(s)*(a0(s)*W1(s) + c0(s)*W0(s));
        W0(s) = W0new;
        BM(s,t) = g(s)*(BM(s-1,t) + h(s)*W1(s));
    end
    
%     dlmwrite ( 'Carlog/W0.txt' , W0,'delimiter','\t','-append');
%     dlmwrite ( 'Carlog/W1.txt' , W1,'delimiter','\t','-append');
%     dlmwrite ( 'Carlog/g.txt' , g,'delimiter','\t','-append');

    % DIHC
    if t==1
        BM_hpf(:,t) = q*(zeros(length(BM_hpf(:,t)),1) + BM(:,t) - zeros(length(BM(:,t)),1));
    else
        BM_hpf(:,t) = q*(BM_hpf(:,t-1) + BM(:,t) - BM(:,t-1));
    end
    z = (BM_hpf(:,t)+0.175);    % p328
    z(z<0)=0;                   % clipping in MATLAB
    v_mem = (z.^3)./(z.^3+z.^2+0.1);    % p328
%     dlmwrite ( 'Carlog/v_mem.txt' , v_mem,'delimiter','\t','-append');
    IHC_new = v_mem.*trans;             % p329
    trans = trans + c_in*(1-trans) - c_out*IHC_new;
%     dlmwrite ( 'Carlog/trans.txt' , trans,'delimiter','\t','-append');
    if t==1
        IHCa(:,t) = (1-c_IHC)*zeros(length(IHCa(:,t)),1) + c_IHC*IHC_new;
        IHC(:,t) = (1-c_IHC)*zeros(length(IHCa(:,t)),1) + c_IHC*IHCa(:,t);
    else
        IHCa(:,t) = (1-c_IHC)*IHCa(:,t-1) + c_IHC*IHC_new;
        IHC(:,t) = (1-c_IHC)*IHC(:,t-1) + c_IHC*IHCa(:,t);
    end

    % DOHC
    v_OHC = W1-W1old;
%     sqr=(vOHC*scale+offset).^2; % NOT USED??
    NLF = 1./(1+(scale*v_OHC + offset).^2);
    % the parameter b will be updated in AGC-SF1
%     dlmwrite ( 'Carlog/NLF.txt' , NLF,'delimiter','\t','-append');
    % AGC
    In8 = In8 + IHC(:,t)/8;

    if ~rem(t-1,64)
%         disp(["here1",t]);
        AGC3 = (1-c_AGC(4))*AGC3 + c_AGC(4)*In64;
        % it seems that the cr and the cl should be replaced!!!
        AGC3 = cr(4)*circshift(AGC3,1) + cc(4)*AGC3 + cl(4)*circshift(AGC3,-1);  % LPF in spatial domain
        In64 = 0*In64;
        
    end
    if ~rem(t-1,32)
%         disp(["here2",t]);
        AGC2 = (1-c_AGC(3))*AGC2 + c_AGC(3)*(In32 + 2*AGC3);
        AGC2 = cr(3)*circshift(AGC2,1) + cc(3)*AGC2 + cl(3)*circshift(AGC2,-1);
        In64 = In64 + In32;
        In32 = 0*In32;
    end
    if ~rem(t-1,16)
%         disp(["here3",t]);
        AGC1 = (1-c_AGC(2))*AGC1 + c_AGC(2)*(In16 + 2*AGC2);
        AGC1 = cr(2)*circshift(AGC1,1) + cc(2)*AGC1 + cl(2)*circshift(AGC1,-1);
        In32 = In32 + In16;
        In16 = 0*In16;
    end
    if ~rem(t-1,8)
%         disp(["here4",t]);
        AGC0 = (1-c_AGC(1))*AGC0 + c_AGC(1)*(In8 + 2*AGC1);
        AGC0 = cr(1)*circshift(AGC0,1) + cc(1)*AGC0 + cl(1)*circshift(AGC0,-1);
        AGC(:,t) = AGC0;                                % store AGC output for plotting
        b = AGC0;
        r = r1 + d_rz.*(1-b).*NLF;                        % feedback to BM    
        g = (1-2*a0.*r+r.*r)./(1-(2*a0-h.*c0).*r+r.*r); % gain for BM
        In16 = In16 + In8;
        In8 = 0*In8;
    else
        if t~=1
            AGC(:,t) = AGC(:,t-1);
        end
    end

    
end

%%

% start_p = duration*fs*(num_rep-1);
%IHC_li = zeros(nsec, npoints);

% for n=1:nsec - 1
%         IHC_li(n, :) = IHC(n, :) - IHC(n + 1, :);
% end

end_p = duration*fs*num_rep;
start_p = duration*fs*(num_rep-1);

BM_out = (BM(1:nsec, start_p+1:end_p));
%BM_hpf_out=(BM_hpf(1:nsec, start_p+1:end_p));
% IHC_li_out=(IHC_li(1:nsec, start_p+1:end_p));
% IHC_out=(IHC(1:nsec, start_p+1:end_p));
%BM_li_out=(BM_li(1:nsec, start_p+1:end_p));
% BM_show = my_BM ./ max(max(my_BM));

% figure
% plot(BM_show');
% 
% figure
% plot(IHC(1:6, start_p+1:200:end_p));
% 
% figure
% plot(AGC(1:6, start_p+1:200:end_p));
% 
% fig=figure;
% fig.Units               = 'centimeters';
% fig.Position(3)         = 16;
% fig.Position(4)         = 8;
% imshow(abs(BM_show(:,1:20:end)))
% colormap('hot')

end