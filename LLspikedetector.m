function [ets,ech]=LLspikedetector(d,sfx,llw,prc,badch)
%jon.kleen@ucsf.edu 2016-2021
% Transforms data into linelength then detects events (spikes) surpassing
% the designated percentile threshold. Note that this function assumes any 
% detections in any channel occurring simultaneously are involved in the 
% same spike event. 
% Based on Estellar et al 2001, DOI 10.1109/IEMBS.2001.1020545
%INPUTS
  % d: vector or matrix of ICEEG data and 
  % sfx: sampling frequency
  % llw: linelength window (in seconds) over which to calculate transform
  % prc: percentile to use as a threshold for detections
  % badch: logical index of bad channels (1=bad, 0=ok)
% OUTPUTS
  % ets: matrix of events (rows) and their on/off times (2 columns) in samples
  % ech: logical index of which channels are involved in each detection, 
%        thus having the same number of rows (spikes) as ets
  
%Example: [ets,ech]=LLspikedetector(d,512,.04,99.99)

if ~exist('llw','var'); llw=.04; end %default linelength window for transform is 40ms
if ~exist('prc','var'); prc=99.5; end %default percentile is 99.9%
if length(size(d))>2; error('Accepts only vector or 2-D matrix for data'); end
if size(d,1)>size(d,2); d=d'; end %flip if needed for loop (assumes longer dimension is time)
if ~exist('badch','var'); badch=false(1,size(d,1)); end %default: all channels ok


%%  1. LINE-LENGTH TRANSFORM
% Will be same size as d with tail end (<llw) padded with NaNs.
numsamples=round(llw*sfx); % number of samples in the transform window
if any(size(d)==1)   %if d is a vector
  L=nan(1,length(d)); % will fill this with transformed data in loop below
  for i=1:length(d)-numsamples
    L(i)=sum(abs(diff(d(i:i+numsamples-1))));
  end 
else                 %if d is a matrix
  L=nan(size(d));
  for i=1:size(d,2)-numsamples
  L(:,i)=sum(abs(diff(d(:,i:i+numsamples-1),1,2)),2);
  end
end


%%  2. DETECT EVENTS
%vectorized version of L for threshold computation below
Lvec=reshape(L,1,numel(L)); Lvec(isnan(Lvec))=[]; 
%get raw LL threshold
thresh=prctile(Lvec,prc); 
%copy of L as a thresholding index (logical index)
Li=L>thresh; 
% Consolidate (sum of index 1's) across channels into a vector,
% then index for any times when at least one channel is above threshold.
a=nansum(Li,1)>0; 
a=diff(a);
eON=find(a==1)+1; % add 1 because diff function shifts back by 1
eOFF=find(a==-1);
  if eOFF(1)<eON(1); eON=[1 eON]; end %correct for any that seem to start prior to beginning of the data
  if length(eOFF)<length(eON); eOFF=[eOFF length(a)]; end %correct for any that seem to end after the end of the data
  if length(eOFF)~=length(eON); error('start and end of events is not matching up, check your code'); end
ets(:,[1 2])=[eON(:) eOFF(:)];
% Now we have rows= events, and 1st column is onset and 2nd column is offset

%Channel indexing for each event
ech=false(size(ets,1),size(Li,1));
for i=1:size(ets,1)
  ech(i,:)=logical(nansum(Li(:,ets(i,1):ets(i,2)),2)); 
end


%%  3. Additional checks/corrections
%Add half of LL window (sec) to event times, centersing LL transform window
ets=round(ets+(sfx*llw)/2); 

% Bad channels shouldn't contribute to detections, mark as zero's, 
ech(:,badch)=0; 
% then remove events that consisted only of bad channels
idx=sum(ech,2)<1;    ech(idx,:)=[];    ets(idx,:)=[]; clear idx

% Merge any events close in time (offset to onset <300ms) as same event
s=size(ets,1); indx=false(s,1); 
for i=1:s-1
  if (ets(i+1,1)-ets(i,2))<sfx*.3
    ets(i+1,1)=ets(i,1); 
    ech(i+1,:)=logical(sum(ech(i:i+1,:),1)); 
    indx(i)=true;
  end
end
ets(indx,:)=[];  ech(indx,:)=[]; % remove merged instances

% Lastly, impose minimum total spike event detection duration
minL=.025; %default 250ms (adjust per preference)
tooshort=diff(ets,1,2)<(sfx*minL);
ets(tooshort,:)=[];  ech(tooshort,:)=[]; clear tooshort

% Tada
