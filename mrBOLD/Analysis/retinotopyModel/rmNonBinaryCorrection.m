function stim = rmNonBinaryCorrection(stim, params)
%   A scan had 96 frames + 5 pre-scan (removed) frames.
%   The subject fixated 3 deg to the right and above fixation for the
%   entire scan.
%
%   nFrames = 101;
%   deg = 3;
%   x = ones(1, nFrames) * deg;
%   y = ones(1, nFrames) * deg;
%   save Stimuli/jitter3deg101frames x y;
%
% Example to simulate random eye movements with a standard deviation of 1
% deg:
%
%


%% If no eye movement data, return without doing anything
if ~checkfields(stim, 'jitterFile'), return; end
[p, n, e] = fileparts(stim.jitterFile);
if strcmpi(n,'none'), return; end  % A filename of None returns

%% If no nonbinary Eyetracker correction should be performed - do nothing
%if ~strcmp(stim.imFilter,'energy'); return; end
%% Build jittered images

load(stim.jitterFile)

% Parse eye eyePosition data
%[x, y] = parseEyePositionData(stim);

% Get stimulus mesh grid in visual angle
[m n step] = prfSamplingGrid(params);
nrows = size(m,1); ncols = size(m, 2);

%Normalize Data to grid!



%Normalize to [0 1]
rangex=1280-1;
rangey=1024-1;

x = (x-1)/rangex;
y = (y-1)/rangey;

%Normalize to [-50 50]

newrangex=nrows-1;
newrangey=ncols-1;

newminx=(nrows-1)/2;
newminy=(ncols-1)/2;

x = round((x*newrangex)-newminx);
y = round((y*newrangey)-newminy);

warning('[%s]: Attention: y will be inverted',mfilename)

y=-y;

% Initialize an image of the correct 2D dimenstions
im = zeros(size(m));
inStimWindow = params.stim(1).instimwindow;

% TR
TR=stim.framePeriod;

EyetrackerSampleEnd=0;

%% Jitter stimulus frame-by-frame to compensate for eye position
% for f = 1:stim.nFrames
%
%     % Which tracker samples correspond to nFrame f
%
%     EyetrackerSampleStart=EyetrackerSampleEnd+1;
%     % -1 to correct for one based
%     EyetrackerSampleEnd=(EyetrackerSampleStart+TR*1000)-1;
%
%     %check if we have to few samples
%     if EyetrackerSampleEnd>length(x)
%
%         %warning('[%s]: Attention: y will be inverted',mfilename)
%
%         EyetrackerSampleEnd=length(x);
%
%     end
%
%     %Map die Eye Position Distribution
%     EyePosDist{f} = zeros(size(m));
%
%     for i=EyetrackerSampleStart:EyetrackerSampleEnd
%
%         EyePosDist{f}(x(i)+51,y(i)+51)=EyePosDist{f}(x(i)+51,y(i)+51)+1;
%
%     end
%
%     EyePosDistVector{f}=EyePosDist{f}(:);
%
%     %pause(0.3)
%     %figure
%     %imagesc(EyePosDist)
%
%
%
%     %imWeighted=im;
%
%
%
%     % Reshape from 2D to 1D
%     %stim.images(:, f) = imWeighted(inStimWindow);
%
%     %clear EyePosDist
% end



for FrameNum=1:stim.nFrames
    
    %DETERMINE EYE POSITION
    
    % Which tracker samples correspond to nFrame f
    
    EyetrackerSampleStart=EyetrackerSampleEnd+1;
    % -1 to correct for one based
    EyetrackerSampleEnd=(EyetrackerSampleStart+TR*1000)-1;
    
    %check if we have to few samples
    if EyetrackerSampleEnd>length(x)
        
        %warning('[%s]: Attention: y will be inverted',mfilename)
        
        EyetrackerSampleEnd=length(x);
        
    end
    
    %Map die Eye Position Distribution
    EyePosDist{FrameNum} = zeros(size(m));
    
    for i=EyetrackerSampleStart:EyetrackerSampleEnd
        
        EyePosDist{FrameNum}(x(i)+51,y(i)+51)=EyePosDist{FrameNum}(x(i)+51,y(i)+51)+1;
        
    end
    
    
    
    % RESHAPE THE IMAGE
    
    % Reshape image from 1D to 2D
    
    im(inStimWindow) = stim.images(:, FrameNum);
   
    
    % Create all possible jitters for Frame (which are actually present, i.e. nonzero in EyePosDist)
    
    [xEyePos,yEyePos]=find(EyePosDist{FrameNum});
    
    yEyePosNeg=yEyePos-51;
    xEyePosNeg=xEyePos-51;
    
    for f = 1:length(xEyePosNeg)
        for g = 1:length(yEyePosNeg)
            
            
            % Jitter in opp direction to eye movement
            
            imShifted = zeros(size(im));
            
            
            % We can shift in 4 possible ways:
            if xEyePosNeg(f) >= 0
                xShifted = 1:ncols-xEyePosNeg(f); xUnshifted =  xEyePosNeg(f)+1:ncols;
            else
                xShifted = -xEyePosNeg(f)+1:ncols; xUnshifted =  1:ncols+xEyePosNeg(f);
            end
            
            if yEyePosNeg(g) >= 0
                yShifted = 1:nrows-yEyePosNeg(g); yUnshifted =  yEyePosNeg(g)+1:nrows;
            else
                yShifted = -yEyePosNeg(g)+1:nrows; yUnshifted =  1:nrows+yEyePosNeg(g);
            end
            
            imShifted(yShifted, xShifted) = im(yUnshifted, xUnshifted);
            
            OriginalShifted{xEyePos(g)}{yEyePos(f)} = imShifted;
            
            clear imShifted xShifted yShifted xUnshifted yUnshifted
            
        end
        
    end
    
    
    %Mean the images
    
    realstimimage=zeros(size(m));
    
    for p = 1:length(xEyePos)
        for q = 1:length(yEyePos)
            
            realstimimage=realstimimage+EyePosDist{FrameNum}(xEyePos(p),yEyePos(q))*OriginalShifted{xEyePos(p)}{yEyePos(q)};
            
        end
    end
    
    
    realstimimage=realstimimage./sum(sum(EyePosDist{FrameNum}));
    
    stim.images(:,FrameNum)=realstimimage(inStimWindow);
    
    clear realstimimage p q f g clear OriginalShifted
    
end

end

