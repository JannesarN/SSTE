function [audio]=evenize(audio, sylBounds, syls)
%%EVENIZE	normalize each syllable to a value below one
% helpful for audio signals including weak utterances
%	
%	audio: input audio
%
%	sylBounds: syllable boundaries
%
%	syls: number of syllables excluding the blank periods


 for i=1:syls
    audio(sylBounds(2*i-1):sylBounds(2*i)) = audio(sylBounds(2*i-1):sylBounds(2*i))./max(abs(audio(sylBounds(2*i-1):sylBounds(2*i))));
 end

end