%% Reservoir of Spiking neurons with FORCE learning and SSTE
%	
%
%
%% Provide Input file

t_sig_file='sample_sequence.mp3';

%%%*** before calling the InputRead! The sylBoundaries function
%		provides approximate whereabouts of the syllable onset/offsets
%		for words longer than single syllable. Provide a manually set sylBounds
%		vector in InputRead for better performance.

%%%*** for sizeHDTS, make sure to choose a divisor of audio length.
%      chech with these codes first: 
% [temp_audio,temp_fs]=audioread(t_sig_file);
% itrLen=(ceil(length(temp_audio)/temp_fs)*temp_fs);
% divisors(itrLen) %% choose from this list
% clear temp_audio temp_fs itrLen

%% Read supervisor audio file and construct its corresponding encoders
wordSyls=1; %number of syllables in each word in the sequence
wordCnt=6;	% number of words in the provided audio sequence file
batchSize=10;	% Number of TE inputs assigned to each distinct syllable
sizeHDTS=500;	% Number of HDTS inputs

tic
[audio, fs, sste, hdts]=InputRead(t_sig_file, 'speech', wordSyls, wordCnt,sizeHDTS, batchSize, 0);
timer_inputRead=toc;

audioLen=length(audio);     % number of samples in the audio signal

dt=1000/fs;                 %time step in ms , will be used for plots, etc.
ttime=audioLen*dt;          %signal duration in ms

%% Initializing SIM params
T= ttime*3;         % Total simulation time = x-times the signal duration *** change x manually to reconfigure
nt = round(T/dt);   % number of time steps 

rng('default');     % seed varies based on system time.
%rng(seed);         %replace seed with a constant value to generate the same random vectors at all executions 


%% Initializing Izhikevich Parameters and Reservoir Neurons
N =  800;                      %number of neurons *** change manually to reconfigure
C = 250;
vr = -60; 
b = 0; 
ff = 2.5; 
vpeak = 30; 
vreset = -65;
vt = -40;
Er = 0;
u = zeros(N,1); 
a = 0.002;
d = 100; 
tr = 0;                         % rise time constant
td = 20;                        % decay time constant
p=0.1;
v = vr+(vpeak-vr)*rand(N,1);    % initial neuron voltages
v_ = v;                         % storing previous time step for Euler integration

%% Initializing  synaptic parameter
IPSC = zeros(N,1);              %post synaptic current 
h = zeros(N,1); 
r = zeros(N,1);                 %reservoir output current, input to the readout layer
hr = zeros(N,1);
JD = zeros(N,1);                %current from neurons that spike in this iteration

%% Calling CARFAC and producing cochleagram
nseCar=60;                      % number of CARFAC filters modeling the basilar membrain *** advised number between 60 to 100 (Lyon, 2019)
duration = ttime/1000;          % signal duration in seconds
num_rep=2;                      % number of CAR_FAC repetitions *** helps BM response converge to a more plausible result, always set to 2 and above.

tic
[sig] =carfac(audio,fs,nseCar, duration, num_rep);
timer_carfac=toc;

zx = abs(sig)/max(max(abs(sig)));   % normalize zx and make it positive
zx(isnan(zx)==1)=0;                 % check for nan values

%%  uncomment to display Teacher cochleagram
drawnow
figure('Name','Teacher Cochleagram','NumberTitle','off') 
colormap('hot')
imshow(zx(1:nseCar,1:200:end))
title("Teacher Cochleagram")
colormap('hot')

figure('Name','SSTE','NumberTitle','off') 
plot (0.001*(1:1:audioLen)*dt,sste(:,1:batchSize), 'r.'), hold on
plot (0.001*(1:1:audioLen)*dt,sste(:,batchSize+1:2*batchSize), 'b.'), hold on
plot (0.001*(1:1:audioLen)*dt,sste(:,2*batchSize+1:3*batchSize), 'g.'), hold on
plot (0.001*(1:1:audioLen)*dt,sste(:,3*batchSize+1:4*batchSize), 'c.'), hold on
plot (0.001*(1:1:audioLen)*dt,sste(:,4*batchSize+1:5*batchSize), 'y.'), hold on
plot (0.001*(1:1:audioLen)*dt,sste(:,5*batchSize+1:6*batchSize), 'm.'), hold on
plot (0.001*(1:1:audioLen)*dt,sste(:,6*batchSize+1:7*batchSize), '.' ,'MarkerFaceColor','#FF7800'), hold on
plot(0.001*(1:1:audioLen)*dt,zx(:,:)), hold off
title('Normalized');
title("SSTE")


%% Continute Initializing reservoir parameters
m1= nseCar;
m2= ((wordSyls*wordCnt)+1)*batchSize; %   for SSTE
%m2 = sizeHDTS;                %   for HDTS

G = 1.4*10^3;   % reservoir weights constant coefficient
Q = 1*10^2;     % feedback weights constant coefficient
WE2 = 8*10^2;   % input weights constant coefficient
OMEGA =  G*(randn(N,N)).*(rand(N,N)<p)/(p*sqrt(N)); %Random reservoir weight matrix, will remain unchanged throught training 

z2=sste';
z1 = zeros(m1,1);               %   readout layer output at each iteration
BPhi1 = zeros(N,m1);            %   initialize readout weights
E1 = (2*rand(N,m1)-1)*Q;        %   feedback weight encoders 
E2 = (2*rand(N,m2)-1)*WE2;      %   TE input encoders
BIAS = 1000; %Bias, at the rheobase current.  

Pinv1 = eye(N)*2;               %   The initial correlation weight matirx for RLS method 
step_choices= divisors(fs);                     % step must be a divisor of fs, so ther won't be signal mismatch between the same steps in different audio repetitions
step = max(step_choices(divisors(fs)<=100));    % Training step, *** suitable choice will enable faster convergence without the loss of accuracy
imin = round(10/dt);                            %First step to execute FORCE method, chaotic behavior observed before this 
icrit = round(ttime*2/dt);%round(0.75*T/dt);                       %Last step to execute FORCE method, test phase begins after this 

%% initialize QRD-RLS parameters
Nrls = N;
lambda = 0.1;

for i=1:m1
internalParameters(i).n = Nrls;
internalParameters(i).beta = 1;         % forgetting factor
internalParameters(i).alpha = lambda;   % regularization parametere
internalParameters(i).sigma = ones(Nrls+2,1);
internalParameters(i).r = zeros(Nrls+1,2*Nrls+2);
internalParameters(i).a = zeros(Nrls+1,1);
internalParameters(i).b = zeros(Nrls+1,1);
internalParameters(i).xprime = zeros(Nrls+1,Nrls+1);
internalParameters(i).w = zeros(Nrls, 1);

for j=1:Nrls+1
    internalParameters(i).r(j,j) = lambda;
    internalParameters(i).r(j,j+Nrls+1) = 1/lambda;
end
end

%% initialize storage variables
current =(zeros(nt/step,m1));   % Store the network output at every "step" interval
correl =(zeros(nt/step,1));       % Store the correlation of the previous 200 network outputs with their corresponding zx
RECB = (zeros(nt/step,10));     % Store the first 10 readout connection weights
%  RECV = (zeros(nt/step,10));         % Store some voltage traces 
RECmse=(zeros(nt/step,1));
 
%% initialize counter variables used in the test loop
ilast=1;    %   simulation start point, 
ss = 0;     %   Bphi logging step counter
qq = 1;     %   sample counter, counts over audio signal, resets each time the entire audio length has been processed
k2 = 0;     %   result logging step counter.
ns=0;       %   keeps total number of spikes

read_start=false;
read_end=false;
%% learning and recall loop

for i = ilast:1:nt 
    if i==1
        tic
    end
    if i==2
        timer_step=toc;
    end
%% EULER INTEGRATE
I = IPSC + E1*z1 + E2*z2(:,qq) +  BIAS; 
v = v + dt*(( ff.*(v-vr).*(v-vt) - u + I))/C ; % v(t) = v(t-1)+dt*v'(t-1)
u = u + dt*(a*(b*(v_-vr)-u)); %same with u, the v_ term makes it so that the integration of u uses v(t-1), instead of the updated v(t)
%% 
index = find(v>=vpeak);         %   find neurons elligible for firing
if ~isempty(index)              %   if there are any neurons firing at this iteration
JD = sum(OMEGA(:,index),2);     %   compute the increase in current due to spiking  
%tspike(ns+1:ns+length(index),:) = [index,0*index+dt*i];    %store spike times (used for raster plot) but takes longer to simulate.  
%ns = ns + length(index);       %   total number of spikes 

end

% implementing the synapse, either single or double exponential
if tr == 0 
    IPSC = IPSC*exp(-dt/td)+   JD*(~isempty(index))/(td);
    r = r *exp(-dt/td) + (v>=vpeak)/td;
else
    IPSC = IPSC*exp(-dt/tr) + h*dt;
    h = h*exp(-dt/td) + JD*(~isempty(index))/(tr*td);  %Integrate the current

    r = r*exp(-dt/tr) + hr*dt; 
    hr = hr*exp(-dt/td) + (v>=vpeak)/(tr*td);
end

%Compute the approximant
 z1 = BPhi1'*r;

%Compute the error (for RLS/ FORCE only)
 err = z1 - zx(:,qq); 

 %% Implement FORCE method  
 if mod(i,step)==1      % update weights once every [step] iterations
   
     if i > imin        % post chaotic activity
      if i < icrit      % pre test phase
      %% RLS  (uncomment if this is prefered. runs faster, not hardware friendly)
%        cd1 = Pinv1*r;
%        BPhi1 = BPhi1 - (cd1*err');        % update readout weights
%        Pinv1 = Pinv1 -((cd1)*(cd1'))/( 1 + (r')*(cd1));
    %% QRD_RLS (uncomment if this is your choice, slower but hardware friendly )
    if (~read_start)
        tic
        read_start=true;
    end
       for j=1:m1
        [internalParameters(j), w_out] = general_qr_opt_func(internalParameters(j), r', zx(j,qq));
        BPhi1(:,j)=w_out;                   % update readout weights
       end
     if (~read_end)
        timer_QRDRLS=toc;
        read_end=true;
    end
      end 
     end 
 end
 
%% reset sample counter 
 if (qq>=audioLen)      
     qq = 0;
 end 
 qq = qq + 1; 



%% apply voltage resets + store the current voltages for next iteration
u = u + d*(v>=vpeak);         %sets u to u+d if v>vpeak, applies to all components. 
v = v+(vreset-v).*(v>=vpeak); %sets v = c if v>=vpeak, add 0 if false, add c-v if true, v+c-v = c
v_ = v;                       % sets v(t-1) = v for the next itteration of loop


%% prepare to compute correlation
zx2= zx(:,1:step:audioLen);

 if mod(i,step) ==1  %don't store every time step, saves time. 
     k2 = k2 + 1;
     current(k2,:) = z1'; 
     RECmse(k2)= mean(err.^2);
 
    %  if (k2<=nt/100-1)
     k3=mod(k2, size(zx2,2));
     if (k3>499)&& (k2>499)
         zx3=zx2(:,k3-499:k3);

     elseif (k3<=499) && (k2>499)
         zx3=[zx2(:,end-500+k3+1:end),zx2(:,1:k3)];

     else
         zx3=zx2(:,(k2-500)*(k2>499)+1:k2);
     end
     cortemp=corrcoef(current((k2-500)*(k2>499)+1:k2,:)', zx3);
     correl(k2)=cortemp(1,2);
    %  end
    %% Record some of the readout weights
    RECB(k2,1:20)=BPhi1(1:20);
    %% And voltages (maybe)
    % RECV(k2,:) = [v(1:5)',u(1:5)']; 
 end
 
 
 %Plot progress
if mod(i,round(step/dt))==1 
    
drawnow

%% plot the spike rasterplot
figure(2) 
%plot(fig_handle2,tspike(1:1:ns,2),tspike(1:1:ns,1),'k.')
%ylim([0,100])

%% output as plotted cochleagram (wave shapes)
figure(3)
plot(0.001*(1:1:k2)*dt*i/k2,current(1:1:k2,1:m1)) 
xlabel('Time')
ylabel('Network Output')
xlim([dt*i/1000-5,dt*i/1000])
ylim([-1,1])

%% Output as cochleagram
figure(4)
ShowThis=current(((k2-1999)*(k2>2000))+1:2:k2,1:m1)';
ShowThis=ShowThis./max(max(abs(ShowThis)));
imshow(ShowThis)
colormap('hot')
title('Network Replay')

%% plot the correlation diagram
figure(5)
plot(0.001*dt*i*(1:1:(k2))/(k2),correl(1:1:(k2)))
title ('Accuracy with '+string(m2) + ' TE Inputs')

figure(8)
plot(0.001*dt*i*(1:1:(k2))/(k2),RECmse(1:1:(k2)))
title ('MSE with '+string(m2) + ' TE Inputs')

%% plot the evolution of readout weights
figure(6)
plot(0.001*dt*i*(1:1:k2)/k2,RECB(1:1:k2,1:10),'r.')
xlabel('Time (s)')
ylabel('Readout Weights')

%% plot voltage trace 
% figure(7)
% plot(0.001*dt*i*(1:1:k2)/k2,RECV(1:1:k2,1:10),'r.')
% xlabel('Time (s)')
% ylabel('Voltage Trace')
end   

end

%% save what you need here
