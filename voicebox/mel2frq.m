function frq = mel2frq(mel)
%MEL2FRQ  Convert Mel frequency scale to Hertz FRQ=(MEL)
%	frq = mel2frq(mel) converts a vector of Mel frequencies
%	to the corresponding real frequencies.
%	The Mel scale corresponds to the perceived pitch of a tone

%	The relationship between mel and frq is given by:
%
%	m = ln(1 + f/700) * 1000 / ln(1+1000/700)
%
%  	This means that m(1000) = 1000
%
%	References:
%
%	  [1] S. S. Stevens & J. Volkman "The relation of pitch to
%		frequency", American J of Psychology, V 53, p329 1940
%	  [2] C. G. M. Fant, "Acoustic description & classification
%		of phonetic units", Ericsson Tchnics, No 1 1959
%		(reprinted in "Speech Sounds & Features", MIT Press 1973)
%	  [3] S. B. Davis & P. Mermelstein, "Comparison of parametric
%		representations for monosyllabic word recognition in
%		continuously spoken sentences", IEEE ASSP, V 28,
%		pp 357-366 Aug 1980
%	  [4] J. R. Deller Jr, J. G. Proakis, J. H. L. Hansen,
%		"Discrete-Time Processing of Speech Signals", p380,
%		Macmillan 1993
%	  [5] HTK Reference Manual p73
%	



%      Copyright (C) Mike Brookes 1998
%      Version: $Id: mel2frq.m,v 1.3 2005/02/21 15:22:13 dmb Exp $
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

frq=700*(exp(mel/1127.01048)-1);
