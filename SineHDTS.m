function hdts=SineHDTS(audioLen, enCount, fs)
%SINEHDTS Generate High dimensional temporal signals as temporal encoder
%inputs. This was proposed in (Nicola & Clopath 2017).
%
%   audioLen: the length of the target audio signal.
%
%   enCount: the number of TE inputs designated to the reservoir.
%
%   fs: the sampling frequency rate in the audio
%
%   hdts: a enCount* audioLen matrix

dt=1000/fs;
ttime= audioLen*dt; % audio length in miliseconds
%***generate a vector the same size as audio, mark each index with the
%equivalent time stamp (in ms). since we are calculating the sinusoid, the
%vector values must be normalized. Hence division by ttime (dt*audioLen)
%***each encoder input will correspond to one upstate, we need the sinusoid
%to complete enCound/2 full periods. Hence the multiplication by enCount
upstate = abs(sin(enCount*pi*((1:1:audioLen)*dt)/ttime)); 
for i = 1:1:enCount
hdts(i,:) = (upstate.*((1:1:audioLen)*dt<i*round(ttime/enCount)).*((1:1:audioLen)*dt>(i-1)*round(ttime/enCount)));
end

%% Uncomment to plot hdts for the first few seconds and the first 10 encoders
% y=100;
% time_lim= y*ttime/100; %display only the first y percent of the signal of the signal
% x=1:time_lim/dt;
% for i=1:time_lim/dt
%      hdts_show(:,i)=hdts(:,i)+cumsum(ones(enCount,1));
% end
% 
% figures
% plot(x.*dt/1000,hdts_show)
end