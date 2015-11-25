function [m,v,w,g,f,pp,gg]=gaussmix(x,c,l,m0,v0,w0)
%GAUSSMIX fits a gaussian mixture pdf to a set of data observations [m,v,w,g,f]=(x,xv,l,m0,v0,w0)
%
% Inputs: n data values, k mixtures, p parameters, l loops
%
%     X(n,p)   Input data vectors, one per row.
%     c(1,p)   Minimum variance (can be a scalar if all components are identical or 0 if no minimum).
%              Use [] to take default value of var(x)/n^2
%     L        The integer portion of l gives a maximum loop count. The fractional portion gives
%              an optional stopping threshold. Iteration will cease if the increase in
%              log likelihood density per data point is less than this value. Thus l=10.001 will
%              stop after 10 iterations or when the increase in log likelihood falls below
%              0.001.
%              As a special case, if L=0, then the first three outputs are omitted.
%              Use [] to take default value of 100.0001
%     M0(k,p)  Initial mixture means, one row per mixture.
%     V0(k,p)  Initial mixture variances, one row per mixture.
%     W0(k,1)  Initial mixture weights, one per mixture. The weights should sum to unity.
%
%     Alternatively, if initial values for M0, V0 and W0 are not given explicitly:
%
%     M0       Number of mixtures required
%     V0       Initialization mode:
%                'f'    Initialize with K randomly selected data points [default]
%                'p'    Initialize with centroids and variances of random partitions
%                'k'    k-means algorithm ('kf' and 'kp' determine initialization of kmeans)
%                'h'    k-harmonic means algorithm ('hf' and 'hp' determine initialization of kmeans)
%              Mode 'hf' generally gives the best results but 'f' [the default] is faster
%
% Outputs: (Note that M, V and W are omitted if L==0)
%
%     M(k,p)   Mixture means, one row per mixture. (omitted if L==0)
%     V(k,p)   Mixture variances, one row per mixture. (omitted if L==0)
%     W(k,1)   Mixture weights, one per mixture. The weights will sum to unity. (omitted if L==0)
%     G       Average log probability of the input data points.
%     F        Fisher's Discriminant measures how well the data divides into classes.
%              It is the ratio of the between-mixture variance to the average mixture variance: a
%              high value means the classes (mixtures) are well separated.
%     PP(n,1)  Log probability of each data point
%     GG(l+1,1) Average log probabilities at the beginning of each iteration and at the end

%  Bugs/Suggestions
%     (2) Allow processing in chunks by outputting/reinputting an array of sufficient statistics
%     (3) Implement full covariance matrices
%     (5) Should scale before finding initial centres
%     (6) Other initialization options:
%              's'    scale dimensions to equal variance when initializing
%              'l'    LBG algorithm
%              'm'    Move-means (dog-rabbit) algorithm

%      Copyright (C) Mike Brookes 2000-2006
%      Version: $Id: gaussmix.m,v 1.7 2006/09/04 20:36:14 dmb Exp $
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

[n,p]=size(x);
x2=x.^2;            % need x^2 for variance calculation
if ~length(c)
    c=var(x,1)/n^2;
end
if ~length(l)
    l=100+1e-4;
end
if nargin<6             % no initial values specified for m0, v0, w0
    k=m0;
    if n<=k             % each data point can have its own mixture
        m=x(mod((1:k)-1,n)+1,:);    % just include all points several times
        v=zeros(k,p);               % will be set to floor later
        w=zeros(k,1);
        w(1:n)=1/n;
        if l>0
            l=0.1;            % no point in iterating
        end
    else
        if nargin<5
            v0='f';         % default initialization mode
        end
        m=zeros(k,p);
        v=ones(k,p);
        w=repmat(1/k,k,1);
        if any(v0=='k')                     % k-means initialization
            if any(v0=='p')
                [m,e,j]=kmeans(x,k,'p');
            else
                [m,e,j]=kmeans(x,k,'f');
            end
            for i=1:k
                v(i,:)=var(x(j==i,:),1);
            end    
        elseif any(v0=='h')                     % k-harmonic means initialization
            if any(v0=='p')
                [m,e,j]=kmeanhar(x,k,'p');
            else
                [m,e,j]=kmeanhar(x,k,'f');
            end
            for i=1:k
                v(i,:)=var(x(j==i,:),1);
            end    
            
        elseif any(v0=='p')                  % Initialize using a random partition
            ix=ceil(rand(1,n)*k);       % allocate to random clusters
            ix(rnsubset(k,n))=1:k;      % but force at least one point per cluster
            for i=1:k
                m(i,:)=mean(x(ix==i,:),1);
            end
        else                                % Forgy initialization: choose k random points [default] 
            m=x(rnsubset(k,n),:);         % sample k centres without replacement
        end
    end
else
    k=size(m0,1);
    m=m0;
    v=v0;
    w=w0;
end
if length(c)>1          % if c is a row vector, turn it into a full matrix so it works with max()
    c=c(ones(k,1),:);
end
v=max(v,c);         % apply the lower bound



% If data size is large then do calculations in chunks

memsize=voicebox('memsize'); 
nb=min(n,max(1,floor(memsize/(8*p*k))));    % chunk size for testing data points
nl=ceil(n/nb);                  % number of chunks
jx0=n-(nl-1)*nb;                % size of first chunk

im=repmat(1:k,1,nb); im=im(:);
th=(l-floor(l))*n;
sd=(nargout > 3*(l~=0)); % = 1 if we are outputting log likelihood values
l=floor(l)+sd;   % extra loop needed to calculate final G value

lpx=zeros(1,n);             % log probability of each data point
wk=ones(k,1);
wp=ones(1,p);
wnb=ones(1,nb);
wnj=ones(1,jx0);

% EM loop

g=0;                           % dummy initial value for comparison
gg=zeros(l+1,1);
ss=sd;                       % initialize stopping count (0 or 1)
for j=1:l
    g1=g;                    % save previous log likelihood (2*pi factor omitted)
    m1=m;                       % save previous means, variances and weights
    v1=v;
    w1=w;
    vi=v.^(-1);                 % calculate quantities that depend on the variances
    vm=sqrt(prod(vi,2)).*w;
    vi=-0.5*vi;
    
    % first do partial chunk
    
    jx=jx0;
    ii=1:jx;
    kk=repmat(ii,k,1);
    km=repmat(1:k,1,jx);
    py=reshape(sum((x(kk(:),:)-m(km(:),:)).^2.*vi(km(:),:),2),k,jx);
    mx=max(py,[],1);                % find normalizing factor for each data point to prevent underflow when using exp()
    px=exp(py-mx(wk,:)).*vm(:,wnj);  % find normalized probability of each mixture for each datapoint
    ps=sum(px,1);                   % total normalized likelihood of each data point
    px=px./ps(wk,:);                % relative mixture probabilities for each data point (columns sum to 1)
    lpx(ii)=log(ps)+mx;
    pk=sum(px,2);                   % effective number of data points for each mixture (could be zero due to underflow)
    sx=px*x(ii,:);
    sx2=px*x2(ii,:);
    ix=jx+1;
    
    for il=2:nl
        jx=jx+nb;        % increment upper limit
        ii=ix:jx;
        kk=repmat(ii,k,1);
        py=reshape(sum((x(kk(:),:)-m(im,:)).^2.*vi(im,:),2),k,nb);
        mx=max(py,[],1);                % find normalizing factor for each data point to prevent underflow when using exp()
        px=exp(py-mx(wk,:)).*vm(:,wnb);  % find normalized probability of each mixture for each datapoint
        ps=sum(px,1);                   % total normalized likelihood of each data point
        px=px./ps(wk,:);                % relative mixture probabilities for each data point (columns sum to 1)
        lpx(ii)=log(ps)+mx;
        pk=pk+sum(px,2);                   % effective number of data points for each mixture (could be zero due to underflow)
        sx=sx+px*x(ii,:);
        sx2=sx2+px*x2(ii,:);
        ix=jx+1;
    end
    g=sum(lpx);                    % total log probability summed over all data points
    gg(j)=g;
    w=pk/n;                         % normalize to get the weights
    if pk                       % if all elements of pk are non-zero
        m=sx./pk(:,wp);
        v=sx2./pk(:,wp);
    else
        wm=pk==0;                       % mask indicating mixtures with zero weights
        [vv,mk]=sort(lpx);             % find the lowest probability data points
        m=zeros(k,p);                   % initialize means and variances to zero (variances are floored later)
        v=m;
        m(wm,:)=x(mk(1:sum(wm),:));                % set zero-weight mixture means to worst-fitted data points
        wm=~wm;                         % mask for non-zero weights
        m(wm,:)=sx(wm,:)./pk(wm,wp);  % recalculate means and variances for mixtures with a non-zero weight
        v(wm,:)=sx2(wm,:)./pk(wm,wp);
    end
    v=max(v-m.^2,c);                % apply floor to variances
    
    if g-g1<=th && j>1
        if ~ss break; end  %  stop 
        ss=ss-1;       % stop next time
    end
    
end
if sd  % we need to calculate the final probabilities
    pp=lpx'-0.5*p*log(2*pi);   % log of total probability of each data point
    gg=gg(1:j)/n-0.5*p*log(2*pi);    % average log prob at each iteration
    g=gg(end);
    %     gg' % *** DEBUG ***
    m=m1;       % back up to previous iteration
    v=v1;
    w=w1;
    mm=sum(m,1)/k;
    f=prod(sum(m.^2,1)/k-mm.^2)/prod(sum(v,1)/k);
end
if l==0         % suppress the first three output arguments if l==0
    m=g;
    v=f;
    w=pp;
end    
