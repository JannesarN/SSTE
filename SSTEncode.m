function [sstep,sstes]=SSTEncode(sylBounds, audioLen, batchSize, numShared, shareMtx)
%SINEHDTS Generate Syllable Specific Temporal Encoders
%
%   sylBounds: a vector containing syllable boundary indexes.
%
%   audioLen: the length of the target audio signal.
%
%   batchSize: the number of TE inputs corresponding to the same syllable.
%
%   numShared: the number of shared syllables (for SSS encoding)
%
%   shareMtx: if any shared syllables use the mapping in share matrix to
%   assign encoders.
%   EXAMPLE: for a sequence containing only the words export, expert,and
%   report share matrix is as follows:
%   shareMtx= {'1','2';
%              '1','0';
%              '0','2'}
%   *** each number correspond to a distinct shared syllable, and '0'
%   indicates a single instance of a syllable throughout the entire
%   sequence. The maximum number in shareMtx must equal "numShared", and
%   must occur more than once in the shareMtx, otherwise there must have
%   been a mistake.

if nargin<4
    numShared=0;
end

if numShared
    %% Shared Syllable Support
    % coming soon
else
    %% Shared Syllable Blind
    sylCount=length(sylBounds)/2;
    sstep= zeros(audioLen,batchSize*(sylCount+1)); %initializing sst encoders
    syls=batchSize*(sylCount+1);    
    
    for i=1:sylCount
        sstep(sylBounds(2*i-1):sylBounds(2*i), 1+(i-1)*batchSize:batchSize*(i))=1;       
    end
    
    for i=1:sylCount-1
        sstep(sylBounds(2*i)+1:sylBounds(2*i+1)-1,syls-batchSize+1:syls )=1;
    end
        sstep(1:sylBounds(1)-1,syls-batchSize+1:syls )=1;
        sstep(sylBounds(end)+1:audioLen,syls-batchSize+1:syls )=1;    
    
end
%% Uncomment to plot SSTE over audio
% will write when I feel like it
end