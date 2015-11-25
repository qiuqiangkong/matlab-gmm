function x=irfft(y,n)
%IRFFT    Inverse fft of a conjugate symmetric spectrum X=(Y,N)
% Y contains FIX(1+N/2) complex samples from the spectrum: if argument N
% is specified then Y will be truncated or padded accordingly
% IMPORTANT: If N is odd, it MUST be specified explicitly.
%
% See also RFFT



%      Copyright (C) Mike Brookes 1998
%      Version: $Id: irfft.m,v 1.3 2005/02/21 15:22:12 dmb Exp $
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

fl=size(y,1)==1;
if fl y=y(:); end
[m,k]=size(y);
if nargin<2 n=2*m-2;
else
  mm=1+fix(n/2);
  if mm>m y=[y; zeros(mm-m,k)];
  elseif mm<m y(mm+1:m,:)=[];
  end
  m=mm;
end
if rem(n,2)		% odd case
  x=real(ifft([y;conj(y(m:-1:2,:))]));
else			% even case
  y(m,:)=real(y(m,:));	% force nyquist element real
  w=ones(1,k);
%  t=[cumprod([-0.5i; exp(2i*pi/n)*ones(m-2,1)]); 0.5i];
  t=-0.5i* exp((2i*pi/n)*(0:m-1)).';
  z=(t(:,w)+0.5).*(conj(flipud(y))-y)+y;
  z(m,:)=[];
  zz=ifft(z);
  x=zeros(n,k);
  x(1:2:n,:)=real(zz);
  x(2:2:n,:)=imag(zz);
end

if fl x=x.'; end
