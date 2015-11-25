function figbolden(pos,pv)
% embolden the current figure
% Inputs: pos = [xmin ymin width height] gives the position of the lower left corner
%                                        and the window size (default = leave where it is already)
%         pv is a cell name containing attribute-value pairs. 
%            default = {'FontName' 'Ariel'; 'FontSize' 16; 'LineWidth' 2; 'MarkerSize' 8}

%      Copyright (C) Mike Brookes 2003
%      Version: $Id: figbolden.m,v 1.4 2005/04/21 07:12:43 dmb Exp $
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

ps={'Title' 'XLabel' 'YLabel' 'Children'};
if nargin<2
pv={'FontName' 'Ariel'; 'FontSize' 16; 'LineWidth' 2; 'MarkerSize' 8};
end
if nargin<1
    pos=[];
end
scsz=get(0,'screensize');
if length(pos)
    po=get(gcf,'position');
    if length(pos)>2
        po(1:3)=pos(1:3);
        if length(pos)>3
            po(4)=pos(4);
        else
            po(4)=0.75*po(3);
        end     
    else
        po(3)=pos(1);      
        if length(pos)>1
            po(4)=pos(2);
        else
            po(4)=0.75*po(3);
        end
    end
    set(gcf,'position',po);
end
hlist=get(gcf,'children');
while length(hlist)
    pl=get(hlist(1));
    %fprintf('list length = %d, handle = %f\n',length(hlist),hlist(1));
    for i=1:size(pv,1)
        if isfield(pl,pv{i,1})
            set(hlist(1),pv{i,1},pv{i,2})
            %fprintf('set %f %s\n',hlist(1),pv{i,1});
        end
    end
    for i=1:length(ps)
        if isfield(pl,ps{i})
            hlist=[hlist; get(hlist(1),ps{i})];
            %fprintf('add %f:%s\n',hlist(1),ps{i});
        end
    end
    hlist(1)=[];
end

