function [tt,f,b]=spgrambw(data,fs,bw)
%SPGRAMBW Draw grey-scale spectrogram [T,F,B]=(DATA,FS,BW)


%      Copyright (C) Mike Brookes 1997
%      Version: $Id: spgrambw.m,v 1.3 2005/02/21 15:22:14 dmb Exp $
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

if nargin<1
	error('Usage: SPGRAMBW(data,fs,bw)');
end
if nargin<2 fs=11025; end
if nargin<3 bw=200; end
fftmin = 256;
winlen = fix(2*fs/bw);
fftlen = max([winlen fftmin]);
win = [hamming(winlen) ; zeros(fftlen-winlen,1)]; 
win = win/sum(win);
windel =  (0:(length(win)-1)) * win;
ntime = 200;
overlap = fix(max(fftlen/2, (ntime*fftlen-length(data))/(ntime-1)));
ntime=fix((length(data)-overlap)/(fftlen-overlap));
c1=(1:fftlen)';
r1=(0:ntime-1)*(fftlen-overlap);
b=rfft(data(c1(:,ones(1,ntime))+r1(ones(fftlen,1),:)).*win(:,ones(1,ntime)));
f=(0:fftlen/2)*fs/fftlen;
b = b.*conj(b);
t = (r1+windel)/fs;
lim = max(b(:))*0.0001;
b=10*log10(max(b,lim));
imh = imagesc(t,f/1000,b);
axis('xy');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
map = (0:63)'/63;
colormap([map map map]);
colorbar;
if(nargout>0) tt=t; end
