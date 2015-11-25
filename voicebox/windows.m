function w = windows(wtype,n,mode,p)
%WINDOWS Generate a standard windowing function (TYPE,N,MODE,P)
%
%	TYPE	is one of:
%
%			'blackman'
%			'kaiser'	   with parameter P (often called beta) (default P=8)
%			'gaussian'	truncated at +-P std deviations (default P=3)
%			'hamming'
%			'hanning'
%			'harris3'	3-term balckman-harris with 67dB sidelobes
%			'harris4'	4-term balckman-harris with 92dB sidelobes
%			'rectangle'
%			'triangle'
%
%	   N is either the number of points to generate. If N>1 is
%		a non-integer, then FLOOR(N) points will be generated.
%
%	   MODE determines the scaling and sampling of the window function and
%     is a text string with up to 3 characters whose meanings are given below. The
%     default is 'ubw' for window functions whose end points are non-zero and 'unw'
%     for window functions whose end points are zero (e.g. hanning window)
%
%         scaling:
%                  u = unscaled  with the peak of the underlying continuous
%                      window equalling unity. [default]
%                  p = scaled to make the actual peak unity
%                  d = scaled to make unity DC gain (summed sample values).
%                  e = scaled to make unity energy (summed squared sample values).
%
%         first and last samples (see note on periodicity below):
%                  b [both]    = The first and last samples are at the extreme ends of
%                                the window [default for most windows].
%                  n [neither] = The first and last samples are one sample away from the ends
%                                of the window [default for windows having zero end points].
%                  s [shifted] = The first and last samples are half a sample away from the
%                                ends of the window .
%                  l [left]    = The first sample is at the end of the window while the last
%                                is one sample away from the end .
%                  r [right]   = The first sample is one sample away from the end while the
%                                last is at the end of the window .
%       
%         whole/half window (see note on periodicity below):
%                  w = The whole window is included [default]
%                  c = The first sample starts in the centre of the window
%                  h = The first sample starts half a sample beyond the centre
%
% Periodicity:
%     The underlying period of the window function depends on the chosen mode combinations and
%     is given in the table below. For overlapping windows with perfect reconstruction choose
%     N to be an integer and modes 'ws', 'wl' or 'wr'.
%
%        Whole/half window -->     w         c         h
%
%        End points:       b      N-1      2N-1      2N-2
%                          n      N+1      2N+1       2N
%                          s       N        2N       2N-1
%                          l       N       2N+1       2N
%                          r       N       2N-1      2N-2
%                  

%      Copyright (C) Mike Brookes 2002-2005
%      Version: $Id: windows.m,v 1.4 2005/10/14 07:29:19 dmb Exp $
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

kk=[-1 1 1 -1; 0 0 2 -2; 0 1 2 -1;   % mode  w,  h,  c  [normal windows]
   -1 0 1 0; 0 0 2 0; 0 1 2 1;       % modes lw, lh, lc
   -1 2 1 0; 0 0 2 -2; 0 1 2 -1;     % modes rw, rh, rc
   -1 1 1 -1; 0 0 2 -2; 0 1 2 -1;    % modes bw, bh, bc
   -1 1 1 1; 0 0 2 0; 0 1 2 1;       % modes nw, nh, nc
   -1 1 1 0; 0 0 2 -1; 0 1 2 0;];    % modes sw, sh, sc

if nargin<3 | length(mode)==0 | ~ischar(mode)
   mode='uw';
end;
mm=zeros(1,length(mode)+1);
ll='hc lrbns';
for i=1:8
   mm(mode==ll(i))=i-3;
end
wtype=lower(wtype);
k=1+3*max(mm)-min(mm);
if k<4
   switch wtype
   case {'hanning','triangle','blackman'}
      k=k+12;
   end
end

% determine the sample points
fn=floor(n);
v=((0:2:2*fn-2)+(kk(k,1)*fn+kk(k,2)))/(kk(k,3)*n+kk(k,4));

% now make the window
np=0;
switch wtype
case 'hanning'
  w = 0.5+0.5*cos(pi*v);

case 'rectangle'
  w = ones(size(v));

case 'triangle'
  w = 1-abs(v);

case 'gaussian'
  if nargin<4, p=3; end;
  w=exp(-0.5*p(1)^2*(v.*v));
  np=1;

case 'kaiser'
  if nargin<4, p=8; end;
  w=besseli(0,p*sqrt(1-v.^2))/besseli(0,p(1));
  np=1;

case 'hamming'
  w = 0.54+0.46*cos(pi*v);

case 'blackman'
  w = 0.42+0.5*cos(pi*v) + 0.08*cos(2*pi*v);

case 'harris3'
  w = 0.42323 + 0.49755*cos(pi*v) + 0.07922*cos(2*pi*v);

case 'harris4'
  w = 0.35875 + 0.48829*cos(pi*v) + 0.14128*cos(2*pi*v) + 0.01168*cos(3*pi*v);
otherwise
   error(sprintf('Unknown window type: %s', wtype));
end;

% scale if required
if any(mode=='d')
   w=w/sum(w);
elseif any(mode=='e')
   w=w/sqrt(sum(w.^2));
elseif any(mode=='p')
   w=w/max(w);
end

if ~nargout
   plot(w);
   if np>0
         title(sprintf('%s window (%s ) - mode=''%s''',wtype,sprintf(' %g',p(1:np)),mode));
      else
         title(sprintf('%s window - mode=''%s''',wtype,mode));
         end
end;





