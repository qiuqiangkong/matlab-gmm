function b = frq2bark(f)
%FRQ2BARK  Convert Hertz to BARK frequency scale BARK=(FRQ)
%       bark = frq2bark(frq) converts a vector of frequencies (in Hz)
%       to the corresponding values on the BARK scale.

%   There are many published formulae approximating the Bark scale.
%   We use the one from Traunmuller.

%   [1] H. Traunmuller, Analytical Expressions for the
%       Tonotopic Sensory Scale”, J. Acoust. Soc. Am. 88,
%       1990, pp. 97-100.

%      Copyright (C) Mike Brookes 2006
%      Version: $Id: frq2bark.m,v 1.1 2006/06/09 14:37:27 dmb Exp $
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

b=26.81*f./(1960+f)-0.53;
%b=26.28-52547.6*(1960+f).^(-1);
b=b+0.15*(2-b).*(b<2)+0.22*(b-20.1).*(b>20.1);
