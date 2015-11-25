function f = bark2frq(b)
%BARK2FRQ  Convert the BARK frequency scale to Hertz FRQ=(BARK)
%       frq = bark2frq(bark) converts a vector of frequencies (in BARK)
%       to the corresponding values in Hertz.

%   There are many published formulae approximating the Bark scale.
%   We use the one from Traunmuller.

%   [1] H. Traunmuller, Analytical Expressions for the
%       Tonotopic Sensory Scale”, J. Acoust. Soc. Am. 88,
%       1990, pp. 97-100.

%      Copyright (C) Mike Brookes 2006
%      Version: $Id: bark2frq.m,v 1.1 2006/06/09 14:33:39 dmb Exp $
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

f=52547.6*(26.28-b-(22.11-1.1*b).*(b>20.1)/6.1-(3*b-6).*(b<2)/17).^(-1)-1960;