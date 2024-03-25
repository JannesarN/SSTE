function [syl_bounds] = SylBoundaries(audio, WordCt, WordSyls)
%SYL_BOUNDS determines the boundaries of syllables
%   *** This gives an approximate result, best to refine the output
%   *** manually.
%
%   audio: the intput audio, vector.
%
%   WordCt: the number of words in the audio, needed to determine the 
%   number of boundaries
%
%   WordSyls: The number of syllables in the words, must be equal for all
%   words in the input sequence.
%

if nargin < 3
    disp("WARNING! SylBoundaries takes 3 arguments.")
end

%% Determining word boundaries
thresh=mean([mean(abs(audio))/4,median(abs(audio))]); % a threshold to signify foreground sound
borders=find(abs(audio)>2*thresh);  % all the samples with a value above this threshold

bounds=zeros(2*WordCt*WordSyls,1); %each syllable has a lower and a higher bound
d=2;    %represents the fact that there are two borders to a syllable

for i=2:size(borders,1)
    if borders(i)-borders(i-1)>1000 % minimum distance allowed between significant samples to be considered from distinct syllables. 
        if (d+1<WordCt*2)
            bounds(d)=borders(i-1);
            bounds(d+1)=borders(i);
            d=d+2;
        else
            bounds(1)=borders(1);
            bounds(end)=borders(end);
            
        end
    end  
end
 bounds(1)=borders(1);
 bounds(end)=borders(end);
 
 %% Uncomment to plot word boundaries
%  figure
%  plot(audio)
%  xline(reshape(bounds,1,[]))
%  
 %% Determining syllable boundaries (approximate results.)
 % increase the number of conditional cases for more syllables,
 % or write a more generic code supporting all number of syllables in
 % words.
 if WordSyls==2
 syl_bounds=zeros(length(bounds)*2,1);
 
 for i=1:2:length(bounds)-1
    syl_bounds(2*(i-1)+1)=bounds(i);
    syl_bounds(2*(i-1)+2)=floor((bounds(i)+bounds(i+1))/2)-101;% 1
    syl_bounds(2*(i-1)+3)=ceil((bounds(i)+bounds(i+1))/2)-100;% 2
    % you may manually adjust 1 and 2 two to fit the boundaries better
    syl_bounds(2*(i-1)+4)=bounds(i+1);
 end
 else %for single syllable stuff
     syl_bounds=bounds;
 end
 %% Uncomment to view syllable boundaries
 figure
 plot(audio)
 xline(reshape(syl_bounds,1,[]),'-r')

end