function [audio, fs, sste, hdts]=InputRead(filename, label, wordSyls, wordCt,sizeHDTS, batchSize, numShared)
% INPUTREAD  Read the input file.
%
%   filename: file location. complete address or location in root
%   directory.
%
%   label: Type of audio sequence. e.g. words, alphabet, birdsong.
%   wordSyls: The number of syllables in each word. (should be the same for
%   every word in the sequence. default is 2.
%
%   wordCt: The number of words in the provided audio.
%
%   numShared: The number of shared syllables among the test words. e.g.
%   'EX' in 'EXPORT' and 'EXPERT'. default is 0 activating the the SSB
%   approach: the existence of any share syllables among words is ignored.
%   
%   Note1: Yes the number of syllables in the test audio must be known
%   prior to using. Current version does not support intelligent syllable
%   identification.
%
%   Note2: Use clean audio for accurate results. Since the syllable
%   detection is not intelligent any sound with hight intensity envelope
%   will be detected as a distinct syllable. 
%
%   Note3: If NumSyls does not match the number of syllabic audio events
%   detected in this function, learning will be incorrect.
%
%   Example: InputRead("./files/ABCD.wav", "alphabet", 4)

%% dafault values assignment
if nargin <7
    numShared=0;
end
if nargin <6
    batchSize=4;
end
if nargin <5
    %default is the largest possible HDTS encoder count not exceeding 550
    sizeHDTS=max(divisors(it_len).*(divisors(it_len)<550)); 
end


%% reading audio
[orig_audio,fs]=audioread(filename);

%% making audio length divisible to its sampling rate, for better learning.
itrLen=(ceil(length(orig_audio)/fs)*fs); %the ideal signal length
audio=[orig_audio;zeros(itrLen-length(orig_audio),1)];

%% filtering audio to make it smoother
%audio=lowpass(audio,log(orig_fs)*2, orig_fs);

%% uncomment to plot the original vs. the adjusted audio signal
% figure
% subplot(2,1,1)
% plot(orig_audio)
% subplot(2,1,2)
% plot(audio)

%% retrieving syllable boundaries
%%uncomment to run the approximate syllable boundary finder function
sylBounds= SylBoundaries(audio, wordCt, wordSyls);			% runs well with single syllable utterance sequences only

%%uncomment to fetch your manually constructed syllable boundaries for the provided input 
%syl/bounds=load('[file_path.mat]');


[audio]=evenize(audio, sylBounds, wordCt*wordSyls);

%% Uncomment to plot syllable boundaries
%  figure
%  plot(audio)
%  xline(reshape(sylBounds,1,[]),'-r')

%% Generating Temporal Encoders

hdts=SineHDTS(itrLen, sizeHDTS, fs);

sste=SSTEncode(sylBounds, itrLen, batchSize, numShared); 

%%% fucn3 for binary coding

end