function [gci,goi] = dypsa(s,fs)

%DYPSA   Derive glottal closure instances and openings from speech
%   [gci,goi] = dypsa(s,fs) returns vectors gci and goi indicating samples
%   when glottal closure and opening instances occur in the speech s
%
%   Inputs:
%   s     is the speech signal
%   fs    is the sampling frequncy
%
%   Outputs:
%   gci   is a vector of glottal closure instances
%   gco   is a vector of glottal opening instances derived from the gcis.  The parameter dy_cpfrac
%         determines the ratio of the closed phase and the gco is the corresponding time instance

% Algorithm Parameters
%       The following parameters are defined in voicebox()
%
%   dy_cpfrac=0.3;           % presumed closed phase fraction of larynx cycle
%   dy_cproj=0.2;            % cost of projected candidate
%   dy_cspurt=-0.45;         % cost of a talkspurt
%   dy_dopsp=1;              % Use phase slope projection (1) or not (0)?
%   dy_ewdly=0.0008;         % window delay for energy cost function term [~ energy peak delay from closure] (sec)
%   dy_ewlen=0.003;          % window length for energy cost function term (sec)
%   dy_ewtaper=0.001;        % taper length for energy cost function window (sec)
%   dy_fwlen=0.00045;        % window length used to smooth group delay (sec)
%   dy_fxmax=500;            % max larynx frequency (Hz) 
%   dy_fxmin=50;             % min larynx frequency (Hz) 
%   dy_fxminf=60;            % min larynx frequency (Hz) [used for Frobenius norm only]
%   dy_gwlen=0.0030;         % group delay evaluation window length (sec)
%   dy_lpcdur=0.020;         % lpc analysis frame length (sec)
%   dy_lpcn=2;               % lpc additional poles
%   dy_lpcnf=0.001;          % lpc poles per Hz (1/Hz)
%   dy_lpcstep=0.010;        % lpc analysis step (sec)
%   dy_nbest=5;              % Number of NBest paths to keep
%   dy_preemph=50;           % pre-emphasis filter frequency (Hz) (to avoid preemphasis, make this very large)
%   dy_spitch=0.2;           % scale factor for pitch deviation cost
%   dy_wener=0.3;            % DP energy weighting
%   dy_wpitch=0.5;           % DP pitch weighting
%   dy_wslope=0.1;           % DP group delay slope weighting
%   dy_wxcorr=0.8;           % DP cross correlation weighting
%   dy_xwlen=0.01;           % cross-correlation length for waveform similarity (sec)

%   Revision History: 
%
%   2.6 - 29 Jun 2006  - Tidied up algorithm parameters
%   2.4 - 10 Jun 2006  - Made into a single file aand put into VOICEBOX
%   2.3 - 18 Mar 2005  - Removed 4kHz filtering of phase-slope function 
%   2.2 - 05 Oct 2004  -  dpgci uses the slopes returned from xewgrdel
%                      -  gdwav from speech with fs<9000 is not filtered
%                      -  Various outputs and inputs of functions have been
%                         removed since now there is no plotting

%   Bugs:  Allow the projections only to extend to the end of the larynx cycle

%      Copyright (C) Tasos Kounoudes, Jon Gudnason, Patrick Naylor and Mike Brookes 2006
%      Version: $Id: dypsa.m,v 2.6 2006/06/25 20:54:14 dmb Exp $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   ftp://prep.ai.mit.edu/pub/gnu/COPYING-2.0 or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract algorithm constants from VOICEBOX

dy_preemph=voicebox('dy_preemph');
dy_lpcstep=voicebox('dy_lpcstep');
dy_lpcdur=voicebox('dy_lpcdur');
dy_dopsp=voicebox('dy_dopsp');              % Use phase slope projection (1) or not (0)?
dy_ewtaper=voicebox('dy_ewtaper');        % Prediction order of FrobNorm method  in seconds
dy_ewlen=voicebox('dy_ewlen');        % windowlength of FrobNorm method  in seconds
dy_ewdly=voicebox('dy_ewdly');        % shift for assymetric speech shape at start of voiced cycle
dy_cpfrac=voicebox('dy_cpfrac');        % presumed ratio of larynx cycle that is closed

dy_wpitch=voicebox('dy_wpitch');           % DP pitch weighting
dy_wener=voicebox('dy_wener');           % DP energy weighting
dy_wslope=voicebox('dy_wslope');           % DP group delay slope weighting
dy_wxcorr=voicebox('dy_wxcorr');           % DP cross correlation weighting
dy_DPWeights=[dy_wpitch dy_wener dy_wslope dy_wxcorr];
dy_lpcnf=voicebox('dy_lpcnf');          % lpc poles per Hz (1/Hz)
dy_lpcn=voicebox('dy_lpcn');            % lpc additional poles

lpcord=ceil(fs*dy_lpcnf+dy_lpcn);       % lpc poles

%PreEmphasis
s_used=filter([1 -exp(-2*pi*dy_preemph/fs)],1,s);

%temp
selected_psps = [];

% perform LPC analysis, AC method with Hamming windowing
[ar, e, k] = lpcauto(s_used,lpcord,floor([dy_lpcstep dy_lpcdur]*fs));

if any(any(isinf(ar)))    % if the data is bad and gives infinite prediction coefficients we return with a warning
    warning('No GCIs returned');
    gci=[];
    return;
end;

% compute the prediction residual
r = lpcifilt(s_used,ar,k); 

% compute the group delay function:  EW method described in "Group Delay Paper" : Brookes, Naylor and Gudnason
[zcr_cand,sew,gdwav,toff]=xewgrdel(r,fs); 
gdwav=-[zeros(toff,1); gdwav(1:end-toff)];
zcr_cand=[round(zcr_cand), ones(size(zcr_cand))];   %flag zero crossing candidates with ones

sew=0.5+sew';  %the phase slope cost of each candidate

% projected candidates derived from the smoothed gdwav function
% if fs>=9000
%     [b_lp,a_lp]=butter(12, 2*4000/fs);
%     gdwav = filtfilt(b_lp, a_lp, gdwav);
% end;

pro_cand=[];
if dy_dopsp ~= 0
    pro_cand = psp(gdwav,fs);
    pro_cand = [pro_cand, zeros(length(pro_cand),1)]; %flag projected candidates with zeros
    sew =      [sew zeros(1,size(pro_cand,1))];      %the phase slope cost of a projected candidate is zero
end;

%Sort the zero crossing and projected candidates together and remove any candidates that
%are within 200 samples from the speech boundary
[gcic,sin] = sortrows([zcr_cand; pro_cand],1);  
sew=sew(sin);
sin=find(and(200<gcic,gcic<length(gdwav)-200));
gcic=gcic(sin,:);
sew=sew(sin);

% compute the frobenious norm function used for a cost in the DP
fnwav=frobfun(s_used,dy_ewtaper*fs,dy_ewlen*fs,dy_ewdly*fs);

%Dynamic programming, picks the most likely candidates based on the pitch consistency, energy etc.
[gci] = dpgci(gcic, s_used, sew, fnwav, fs, dy_DPWeights);

%Evaluate goi ... determined to be dy_cpfrac percentage of the larynx cylce away from the last gci
goi=simplegci2goi(gci,dy_cpfrac);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = psp(g,fs)
%PSP  Calculates the phase slope projections of the group delay function
%   Z = PSP(G) computes the 

%   Author(s): P. A. Naylor
%   Copyright 2002 Imperial College of Science Technology and Medicine, London
%   Revision History: 
%       V1.0 July 12th 2002:
%            Nov. 6th 2002 : if statement added to remove "last" midpoint 

g = g(:);

gdot = [diff(g);0];
gdotdot = [diff(gdot);0];

% find the turning points  as follows: [tp_number, index_of_tp, min(1) or max(-1), g(index_of_tp)]
turningPoints = zcr(gdot);
turningPoints = [[1:length(turningPoints)]', turningPoints, sign(gdotdot(turningPoints)), g(turningPoints)];

% useful for debug/plotting
%tplot = zeros(length(g),1);
%tplot(turningPoints(:,1)) = turningPoints(:,2);

% find any maxima which are < 0 
%negmaxima = turningPoints(find(turningPoints(:,3) == -1 & turningPoints(:,4) < 0),:);
negmaxima = turningPoints(find(turningPoints(:,3) == -1 & turningPoints(:,4) < 0 & turningPoints(:,1)~=1),:);  %Change 01.05.2003 JG: The first row can't be included

% find the midpoint between the preceding min and the negative max
nmi = negmaxima(:,1);
midPointIndex = turningPoints(nmi-1,2) + round(0.5*(turningPoints(nmi,2) - turningPoints(nmi-1,2)));
midPointValue = g(midPointIndex);

% project a zero crossing with unit slope
nz = midPointIndex - round(midPointValue);

% find any minima which are > 0 
posminima = turningPoints(find(turningPoints(:,3) == 1 & turningPoints(:,4) > 0),:);

% find the midpoint between the positive min and the following max
pmi = posminima(:,1); 

%Remove last midpoint if it is the last sample
if ~isempty(pmi), if pmi(end)==size(turningPoints,1), pmi=pmi(1:end-1); end; end;

midPointIndex = turningPoints(pmi,2) + round(0.5*(turningPoints(pmi+1,2) - turningPoints(pmi,2)));
midPointValue = g(midPointIndex);

% project a zero crossing with unit slope
pz = midPointIndex - round(midPointValue);

z = sort([nz;pz]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function i = zcr(x, p)
%ZCR  Finds the indeces in a vector to  zero crossings
%   I = ZCR(X) finds the indeces of vector X which are closest to zero-crossings.
%   I = ZCR(X, P) finds indeces for positive-going zeros-crossings for P=1 and
%   negative-going zero-crossings for P=0.

x = x(:);

if (nargin==2)
    if (p==0) 
        z1 = zcrp(x);   % find positive going zero-crossings
    elseif (p==1) 
        z1 = zcrp(-x);  % find negative going zero-crossings
    else
        error('ZCR: invalid input parameter 2: must be 0 or 1');
    end
else
    z1 = [zcrp(x); zcrp(-x)];
end

% find crossings when x==0 exactly
z0 = find( (x(1:length(x))==0) & ([x(2:length(x));0] ~= 0));

% concatenate and sort the two types of zero-crossings
    i = sort([z0; z1]);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function zz = zcrp(xx)  %only used in zcr
    % find positive-going zero-crossing
    z1 = find(diff(sign(xx)) == -2);
    % find which out of current sample or next sample is closer to zero
    [m, z2] = min([abs(xx(z1)), abs(xx(z1+1))], [], 2);
    zz =  z1 -1 + z2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [frob]=frobfun(sp,p,m,offset)

% [frob]=frobfun(sp,p,m)
% 
% sp is the speech signal assumed to be preemphasised
% p  is the prediction order  : recomended to be 1 ms in above paper
% m  is the window length     : recomended to be 1 ms in above paper
% offset is shift for assymetric speech shape at start of voiced cycle -
% default 1.5ms.
%
% this function implements the frobenius norm based measure C defined in:
% C. Ma, Y. Kamp, and L. F. Willems. "A frobenius norm approach to glottal
% closure detection from the speech signal. IEEE Trans. Speech Audio Processing, 2:258 - 265, Apr 1994.
% This equals the squre of the Frobenius norm of the m by p+1 data matrix divided by p+1

%   Author(s): J. Gudnason, P
%   Copyright 2002 Imperial College of Science Technology and Medicine, London
%   Revision History: 
%       V1.0 July 12th 2002:
%            Nov. 6th 2002 : if statement added to remove "last" midpoint 

%force p m and offset to be integers
p=round(p);
m=round(m);
offset=round(offset);

w=(p+1)*ones(1,m+p);
w(1:p)=1:p;
w(m+1:p+m)=p:-1:1;

w=w./(p+1); 
frob=filter(w,1,sp.^2);
frob(1:(round((p+m-1)/2) + offset))=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function goi=simplegci2goi(gci,pr)

gci=round(gci);
maxpitch=max(medfilt1(diff(gci),7));

% calculate opening instants
for kg=1:length(gci)-1
    goi(kg)=gci(kg)+min(pr*(gci(kg+1)-gci(kg)),pr*maxpitch);
end;
kg=kg+1;
goi(kg)=round(gci(kg)+pr*(gci(kg)-gci(kg-1)));  %use the previous pitch period instead
goi=round(goi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tew,sew,y,toff]=xewgrdel(u,fs)

% implement EW group delay epoch extraction

dy_gwlen=voicebox('dy_gwlen');          % group delay evaluation window length
dy_fwlen=voicebox('dy_fwlen');          % window length used to smooth group delay

% perform group delay calculation

gw=2*floor(dy_gwlen*fs/2)+1;            % force window length to be odd
ghw=window('hamming',gw,'s');
ghw = ghw(:);                           % force to be a column (dmb thinks window gives a row - and he should know as he wrote it!)
ghwn=ghw'.*(gw-1:-2:1-gw)/2;            % weighted window: zero in middle

u2=u.^2;
yn=filter(ghwn,1,u2);
yd=filter(ghw,1,u2);
yd(abs(yd)<eps)=10*eps;                 % prevent infinities
y=yn(gw:end)./yd(gw:end);               % delete filter startup transient
toff=(gw-1)/2;
fw=2*floor(dy_fwlen*fs/2)+1;            % force window length to be odd
if fw>1
    daw=window('hamming',fw,'s');
    y=filter(daw,1,y)/sum(daw);         % low pass filter 
    toff=toff-(fw-1)/2;
end
[tew,sew]=zerocros(y,'n');              % find zero crossings

tew=tew+toff;                           % compensate for filter delay and frame advance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gci=dpgci(gcic, s, Ch, fnwav, fs, DPWeights)

%DPGCI   Choose the best Glottal Closure Instances with Dynamic Programming
%   gci=dpgci(gcic, s, Ch, fnwav, fs, DPWeights) returns vectors of sample indeces corresponding
%   to the instants of glottal closure in the speech signal s at sampling frequency fs Hz.
%
%   Inputs:
%   gcic    is a matrix whos first column are the glottal closure instance candidates and
%           the second column is 1 if the corresponding gci is derived from a zero crossing 
%           but zero if the gci is from a a projected zero crossing
%   s       is the speech signal
%   Ch      the phase slope cost of every candidate
%   fnwav   is the frobenious norm function of s
%   fs      is the sampling frequncy
% DPWeights weights of each term in the cost function
%
%   Outputs:
%   gci     is a vector of glottal closure instances chosen by the DP

%   Author(s): J. Gudnason, P.A. Naylor and D.M. Brookes
%   Copyright 2003 Imperial College London
%   Revision History: 
%   Bugs:  Constants are hardwired but defined in a structure like pv (defined in grpdelpv)
%         


%Cost function weights
w = DPWeights;

%Constants
Ncand=length(gcic);
Nbest=voicebox('dy_nbest');        % Number of NBest paths to keep
wproj=voicebox('dy_cproj');           % cost of projected candidate
dy_cspurt=voicebox('dy_cspurt');           % cost of a talkspurt
sv=voicebox('dy_spitch');        % scale factor for pitch deviation cost


%Limit the search:
fxmin=voicebox('dy_fxmin');        % min larynx frequency (Hz)
fxmax=voicebox('dy_fxmax');        % min larynx frequency (Hz)
qrmin=ceil(fs/fxmax);
qrmax=floor(fs/fxmin);

%Cost and tracking r = current, q = previous, p = preprevious
cost=zeros(Ncand, Nbest); cost(:,:)=inf;    %Cost matrix, one row for each candidate
maxcost=zeros(Ncand,1); maxcost(:,:)=inf;   %Maximum cost in each row
imaxcost=ones(Ncand,1);                     %Index of maximum cost

prev = ones(Ncand, Nbest);                  %index of previous, q candidates
ind = ones(Ncand, Nbest);                   %index of p in row q (from prev)
qbest = [zeros(Ncand,1), ones(Ncand,2)]; % the minimum cost in any previous q [cost,q,i]

Cfn=fnrg(gcic(:,1),fnwav,fs);  %Frob.Energy Cost

%Add start and end state
gcic=[[gcic(1,1)-qrmax-2 0];gcic;[gcic(end,1)+qrmax+2 0]];
Cfn=[0 Cfn 0];
Ch = [0 Ch 0];

for r=2:Ncand 
    
    %Penalty for being projected zero crossing
    ZCB=wproj*abs(gcic(r,2)-1);

    for q=r-1:-1:2
                
        if (gcic(r,1)-gcic(q,1)<qrmin) 
            continue;
        end
        if (gcic(r,1)-gcic(q,1)>qrmax)
            cost0=qbest(q,1)+ZCB;
            if cost0 < maxcost(r)   % If the new path is better change elements in cost, prev and ind
                j=imaxcost(r);
                cost(r,j)=cost0;
                prev(r,j)=qbest(q,2);
                ind(r,j)=qbest(q,3);
                
                %Replace maxcost(r) and imaxcost(r)
                [cc, j]=max(cost(r,:));
                maxcost(r)=cc;
                imaxcost(r)=j;
            end;
            break;
            
        else
                        
            %Speech Waveform Similarity Cost           
            Ca=swsc(gcic(q,1),gcic(r,1),s,fs);
            
            for i=1:Nbest

                p=prev(q,i);

                if (gcic(q,1)-gcic(p,1)>qrmax)
                    Cpitch=dy_cspurt;
                else
                    %Pitch Deviation Cost
                    Cpitch=pdc(gcic(p,1),gcic(q,1),gcic(r,1),sv);
                end
                %Cost in state r
                Cpqr=w*[Cpitch; Cfn(r); Ch(r); Ca;]+ZCB;
                
                %Cumulative cost 
                cost0=cost(q,i)+Cpqr;
                
%                 fprintf('%d %d %d %f %f\n',p,q,r,cost(q,i),Cpqr); pause;

                if cost0 < maxcost(r)   % If the new path is better change elements in cost, prev and ind
                    j=imaxcost(r);
                    cost(r,j)=cost0;
                    prev(r,j)=q;
                    ind(r,j)=i;
                    
                    %Replace maxcost(r) and imaxcost(r) in cost
                    [cc, j]=max(cost(r,:));
                    maxcost(r)=cc;
                    imaxcost(r)=j;
                end;
            end;
        end;
    end;
    [cc,j]=min(cost(r,:));
    if cc<qbest(r-1,1)
        qbest(r,:)=[cc,r,j];
    else
        qbest(r,:)=qbest(r-1,:);
    end
end;

%Trace back.
[ctmp,cl]=min(cost(Ncand,:)); %cl gives the total min cost at the last candidate.
row=Ncand;
cnt=1;

Cl = [];
rowpath = [Ncand, cl];
while row > 1
    gci(cnt)=gcic(row,1); 
    nrow=prev(row,cl);
    ncl=ind(row,cl);
    rowpath = [rowpath; [nrow,ncl]]; %save which rows are chosen
    row=nrow;
    cl=ncl;
    cnt=cnt+1;
end;
gci=fliplr(gci);



%%%%%%%%%%%  Cost Functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Cp=pdc(p,q,r,sv)

%PDC Computes the Pitch Deviation Cost
%
%  Cp = pdc(p,q,r)
%  
%  r is the current gci candidate
%  q and p are the previous two propsed gci candidates



mm=1;

y=min(r-q,q-p)/max(r-q,q-p);
Cp=0.5-exp(-(y-mm)^2/(2*sv^2));

%%%%%%%%%%%%%%%%%
function Cfn=fnrg(gcic,frob,fs)
%Frobenious Energy Cost

dy_fxminf=voicebox('dy_fxminf');
frob=frob(:)';
mm=round(fs/dy_fxminf);

mfrob=maxfilt(frob,1,mm);
mfrob=[mfrob(floor(mm/2)+1:end) max(frob(end-ceil(mm/2):end))*ones(1,floor(mm/2))];

rfr=frob./mfrob;
Cfn=0.5-rfr(round(gcic));

%%%%%%%%%%%%%%%%%

function Ca=swsc(q,r,s,fs)

%SWSC Computes the Speech Waveform Similarity Cost
%
%  Ca = swsc(q,r,s)
%  
%  r is the current gci candidate
%  q is the previous propsed gci candidate
%  s is the speech waveform
%  fs is the sample frequency

% fprintf('%d\t%d\n', r, length(s));

%window set to 10 ms

dy_xwlen=voicebox('dy_xwlen');        % window set to 10 ms by default
m=ceil(fs*dy_xwlen);

mmin=floor(m/2);
mmax=floor(m/2)+1;

if r+mmax > length(s)   
    sr=s(r-mmin:end);
    sq=s(q-mmin:q-mmin+length(sr)-1);
    if (isempty(sr))
        display('break');
    end
    if (isempty(sq))
        display('break');
    end
    Ca=-0.5*((sq-mean(sq))'*(sr-mean(sr)))/(std(sq)*std(sr)*m);
else
    sr=s(r-mmin:r+mmax);
    sq=s(q-mmin:q+mmax);
    
    Ca=-0.5*((sq-mean(sq))'*(sr-mean(sr)))/(std(sq)*std(sr)*m);
        
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

