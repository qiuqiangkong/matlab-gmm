function x=usasi(n,fs)
%USASI generates N samples of USASI noise at sample frequency FS X=(N,FS)

% This is based on hearsay: I haven't read the actual standard: EIA-549-1988
% It may also be defined in ANSI standard S1.4-1961: General Purpose Sound Meter
% USASI = USA Standards Institute (now ANSI)
% This generates gaussian noise and filters is with zeros at 0Hz & fs/2 and poles at 100 and 320 Hz
% I need to get the scaling correct



%      Copyright (C) Mike Brookes 1997
%      Version: $Id: usasi.m,v 1.3 2005/02/21 15:22:14 dmb Exp $
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

if nargin<2 fs=8000; end
b=[1 0 -1];
a=poly(exp(-[100 320]*2*pi/fs));

x=randfilt(b,a,n);