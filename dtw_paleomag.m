function [rid,rii,xc,t_out,Dist,w_scale,overlap]=dtw_v4(x,r,tx,tr,g,edge,tr_fix,tx_fix,sd,flag_norm,flag_rmrpt,qplot)

%Dynamic Time Warping Algorithm
%
%rid:    the warped version of r (dec)
%rii:    the warped version of r (inc)
%xc1:    correlation between x and r after warping (dec)
%xc2:    correlation between x and r after warping (inc)
%Dist:  unnormalized distance between t and r
%w_scale: amount the t is warped
%overlapd:  the number of data points that overlap (dec)
%overlapd:  the number of data points that overlap (inc)
%x:     target vector
%r:     vector to be warped  --- needs to be the same length as x, can't have NaN's
%qplot: display output
%tx:    time vector for x
%tr:    time vector for r
%g:     weighting term for going off-diagonal
%edge:  weighting term for the edges
%tr_fix: points in the warping sequence to fix (global constraints)
%xr_fix: point in the target to fix it to (global constraints)
%flag_norm: 1 if you want to normalize the input sections
%flag_rmrpt: 1 if you want to remove repeating points within the warped section
%
% version as of July 31 2018


set(0,'DefaultAxesFontName','Sans Serif')
set(0,'DefaultAxesFontSize',14)

% if nargin<7,
%     g = 1.05;
%     edge = 1;
%     qplot=0;  
%     tr_fix = [];
%     tx_fix = [];
%     flag_norm = 0;
% end;

% save original versions of the time series, xo and ro

xo=x;
ro=r;


% demeans and normalize the input time series to unit std.
% if flag_norm
%     x=normPH(x(:)');
%     r=normPH(r(:)');
% end
x = x';
r = r';
% sort the time series to ensure input is monotonic
[tx,ind] = sort(tx);
x1 = x(1,ind);
x2 = x(2,ind);

[tr,ind] = sort(tr);
r1 = r(1,ind);
r2 = r(2,ind);


[~,N]=size(x1);
[~,M]=size(r1);

% create matrix of the squared difference b/w zero-mean x and r, r - columns, x - rows 
%d=(repmat(x(:),1,M)-repmat(r(:)',N,1)).^2; 
for n=1:N
   for m=1:M
       d(n,m)=1-cosd((angdiff(x1(n), x2(n), r1(m), r2(m))));
   end
end


d(1,:)=d(1,:)*edge;
d(:,end)=d(:,end)*edge;
d(end,:)=d(end,:)*edge;
d(:,1)=d(:,1)*edge;
%d(1,2:end)=d(1,2:end)*10;
%d(2:end,end)=d(2:end,end)*10;

% impose the time constraints 
if ~isempty(tr_fix)
    
    for ii=1:length(tr_fix)
        ind = abs(tx-tx_fix(ii))==min(abs(tx-tx_fix(ii)));
        %sd=100;% standard devation, you can vary this parameter based on confidence in time constraint
        norm = -normpdf(tx,tx(ind),sd(ii));
        
        % find scaling factor, want minimum to be equal to 95% critical value of d
        norm =  -(prctile(d(:),95)/min(norm))*norm;
        
        ind = find(abs(tr-tr_fix(ii))==min(abs(tr-tr_fix(ii))));
        d(:,ind) = d(:,ind) + norm';
    end
end


D=zeros(size(d));
D(1,1)=d(1,1);

% fill the the edges of the cummulative difference matrix, D
if 1,
  for n=2:N
    D(n,1)=d(n,1)+D(n-1,1);
  end
  for m=2:M
    D(1,m)=d(1,m)+D(1,m-1);
  end
end;

%%  adjust this so that you look at the  5 next squares %%

% % run b/w 2:N & 2:M if looking preceeding 3 squares
% for n=2
%   for m=2
%     % check the three preceeding squares to determine square with minimum
%     value 
%     D(n,m)=d(n,m)+min([g*D(n-1,m),D(n-1,m-1),g*D(n,m-1)]);  
%     %D(n,m)=d(n,m)+min([1.01*g*D(n,m-2) g*D(n-1,m-2),D(n-2,m-2),g*D(n-2,m-1) 1.01*g*D(n-2,m)]);  
%   end
% end

% look two steps ahead
n=2;
for m=2:M
    % check the three preceeding squares to determine square with minimum value
    D(n,m)=d(n,m)+min([g*D(n-1,m),D(n-1,m-1),g*D(n,m-1)]);
end

m=2;
for n=2:N
    % check the three preceeding squares to determine square with minimum value
    D(n,m)=d(n,m)+min([g*D(n-1,m),D(n-1,m-1),g*D(n,m-1)]);
end

for n=3:N
  for m=3:M
    % check the 8 preceeding squares to determine square with minimum value 
    D(n,m)=d(n,m)+min([g*D(n-1,m),D(n-1,m-1),g*D(n,m-1),1.1*g*D(n,m-2),1.05*g*D(n-1,m-2),D(n-2,m-2),1.05*g*D(n-2,m-1),1.1*g*D(n-2,m)]);  
  end
end

% start at the end point (N,M) and work your way back along the minimum path
Dist=D(N,M);
n=N;
m=M;
k=1;
w=[];
w(1,:)=[N,M];
while ((n+m)~=2)
    if (n-1)==0  % if you get to the start of the column, move over one row
        m=m-1;
    elseif (m-1)==0 % if you get to the start of the row, move over one column
        n=n-1;
    elseif (n-2)==0 && (m-2)==0
        n=n-1;
        m=m-1;
    elseif (n-2)==0
        [~,number]=min([D(n-1,m),D(n,m-1),D(n-1,m-1),D(n,m-2) D(n-1,m-2)]);
        switch number
            case 1
                n=n-1;  % follow the path up a row
            case 2
                m=m-1; % follow the path to the left
            case 3
                n=n-1; % follow the diagonal path
                m=m-1;
            case 4
                m=m-2;  % follow the path to the left 2
            case 5
                n=n-1; % follow the path to the left 2 and down 1
                m=m-2;
        end
    elseif (m-2)==0
        [~,number]=min([D(n-1,m),D(n,m-1),D(n-1,m-1),D(n-2,m) D(n-2,m-1)]);
        switch number
            case 1
                n=n-1;  % follow the path up a row
            case 2
                m=m-1; % follow the path to the left
            case 3
                n=n-1; % follow the diagonal path 
                m=m-1;
            case 4
                n=n-2;  % follow the path to the down 2
            case 5
                n=n-2; % follow the path to the left 1 and down 2
                m=m-1;
        end        
    else
        %[values,number]=min([D(n-1,m),D(n,m-1),D(n-1,m-1)]);
        [~,number]=min([D(n-1,m),D(n,m-1),D(n-1,m-1),D(n,m-2),D(n-1,m-2),D(n-2,m-2),D(n-2,m-1),D(n-2,m)]);
        switch number
            case 1
                n=n-1;  % follow the path up a row
            case 2
                m=m-1; % follow the path to the left
            case 3
                n=n-1; % follow the diagonal path
                m=m-1;
            case 4
               % m=m-2; % follow the path to the left 2
               n=[n;n]; 
               m=[m-1;m-2];           
            case 5
                %n=n-1; % follow the path to the left 2 and down 1
                %m=m-2;
                n=[NaN;n-1];
                m=[m-1;m-2];
            case 6
                %n=n-2; % follow the diagonal path 2
                %m=m-2;
                n=[n-1;n-2];
                m=[m-1;m-2];
            case 7
                n=n-2; % follow the path to the left 2 and down 2
                m=m-1;
            case 8
                %n=n-2; % follow the path down 2
                n=[n-1;n-2];
                m=[m;m];
        end
    end
    k=k+1;
    w=cat(1,w,[n,m]);
    
    % extract the final value to step forward in time
    n = n(end);
    m = m(end);
end

%% Extract the time and data along the warping path 


% Ensures warped time is monotonic.
% Option to have warped time starting points averaged or not

%Average values of y at repeated values of x, and ensure x is monotonic
if 1 % if you don't want to average over the first time points
    % if you want to do this, you also need to compute a warped time
    % note: currenty not interpolating to get warped time
    %pl=find( (w(:,1)>1) & (w(:,2)>1));  % criteria for not averaging first time point (along the edge)
    pl=find( (w(:,1)>1 & w(:,2)>1) | isnan(w(:,1)));  % also keep the NaN's -- needed to interp depth
    pl=[pl; pl(end)+1];
    
    % find the end indices that are repeating
    ind1 = find(w(:,2)==w(1,2));
    pl = pl(ind1(end):end);
    
    % get the extra points that pile up at the end
    ind2 = find(w(:,1)==w(1,1));
    pl = pl(ind2(end):end);
    
    % remove repeating warped points in the new time series
    if flag_rmrpt
        [~,ind3] = unique(w(:,2),'rows','stable'); % finds the unique rows in w(:,2)
        [~,~,ind] = intersect(ind3,pl);% finds the unique rows that are still included in pl
        pl = pl(ind);
    end
    
    % find the NaN's -- this indices need to have interpolated depths
    ireal = find(~isnan(w(:,1)));
    inan = find(isnan(w(:,1)));
    
    % extract the non-NaN depths
    [~,~,ind] = intersect(ireal,pl);
    
    tmp = NaN(length(pl),2);
    tmp(ind,1) = tx(w(pl(ind),1));
    ro1=ro(:,1);
    ro2=ro(:,2);
    tmp(ind,2) = ro1(w(pl(ind),2));
    tmp(ind,3) = ro2(w(pl(ind),2));

    
    % find the NaN points, interpolate between the adjacent depths
    [~,~,ind] = intersect(inan,pl);
    
    for ii=1:length(ind)
        tmp(ind(ii),3) = ro2(w(pl(ind(ii)),2));  % use the known data points
        tmp(ind(ii),2) = ro1(w(pl(ind(ii)),2));  % use the known data points
        tmp(ind(ii),1) = 0.5*(tx(w(pl(ind(ii))-1,1)) + tx(w(pl(ind(ii))+1,1)));  % interpolate between the depth points
    end
    
    tw = fliplr(tmp(:,1)');
    rwd = fliplr(tmp(:,2)');
    rwi = fliplr(tmp(:,3)');

    % get the extra points at the end  -- only if unique
    append_t = fliplr(tr(w(1:ind2(end),2)));
    append_r1 = fliplr(ro1(w(1:ind2(end),2)));
    append_r2 = fliplr(ro2(w(1:ind2(end),2)));
    
    % get the extra start points that fall before the target time series
    start_t = fliplr(tr(w(pl(end)+2:end,2)));
    start_r1 = fliplr(ro1(w(pl(end)+2:end,2)));
    start_r2 = fliplr(ro2(w(pl(end)+2:end,2)));
    
    % compute the warping scaling factor
    %w_scale = (tx(w(1,1))-tx(w(pl(end),1)))/(tr(w(1,2))-tr(w(pl(end),2)));
    w_scale = (tx(w(1,1))-tx(w(pl(end),1)))/(tr(w(pl(1),2))-tr(w(pl(end),2)));
    
    
    % shift the extra time at the ends, then scale it
    if ~isempty(start_t)
        start_t = w_scale*(start_t - start_t(end)) + tw(1);
    end
    
    if ~isempty(append_t)
        append_t = w_scale*(append_t - append_t(1)) + tw(end);
    end
    
else % include averaging over the first few time points
    %tw=fliplr(tr(w(:,1))); % use the rows of w to get the matching time
    tw=fliplr(tx(w(:,1))); % use the rows of w to get the matching time
    rwd=fliplr(ro1(w(:,2))); % use the cols of w to get the corresponding values from original time series (not zero mean)
    rwi=fliplr(ro2(w(:,2))); % use the cols of w to get the corresponding values from original time series (not zero mean)
end

%initialize matrix
pl=find(diff(tw)~=0);
dj=max(diff([0 pl length(tw)]));
tw2=NaN(dj,length(pl)+1);
rwd2=NaN(dj,length(pl)+1);
rwi2=NaN(dj,length(pl)+1);

%initial values
tw2(1,1)=tw(1);
rwd2(1,1)=rwd(1);
rwi2(1,1)=rwi(1);
ct2=1;
ct3=1;

%loop through all entries
for ct=2:length(tw);
    if tw(ct)~=tw(ct-1),
        ct2=ct2+1;
        ct3=1;
    else
        ct3=ct3+1;
    end;
    tw2(ct3,ct2)=tw(ct);
    rwd2(ct3,ct2)=rwd(ct);
    rwi2(ct3,ct2)=rwi(ct);
end;

%average matrix into a monotonic vector
%tw3=nanmean (tw2);  % causing a problem that the mean of same time is not numerically the same as the extract time
tw3=tw2(1,:);

% average everything but the last point
if dj~=1
    rwd3=nanmean(rwd2);
    rwd3(end) = rwd2(1,end);
    rwi3=nanmean(rwi2);
    rwi3(end) = rwi2(1,end);
else
    rwd3=rwd2;
    rwi3=rwi2;
end

% get the sequences not averaged
r_alld = rwd2(~isnan(rwd2(:)));
r_alli = rwi2(~isnan(rwi2(:))); 
t_all = tw2(~isnan(tw2(:))); 
 
% check dimensions of r_all, should be a row vector
if size(r_alld,1)~=1
    r_alld = r_alld';
    r_alli = r_alli';
    t_all = t_all';
end

%% Compute the correlation coefficient 

% make sure that length of t_out is the same as tx for computing xc
if length(tw3)~=length(tx) 
    [~,ind_tw,ind_tx] = intersect(tw3,tx);
    
    t_out =tw3(ind_tw);
    rid = rwd3(ind_tw);
    rii = rwi3(ind_tw);
    
    tx_out = tx(ind_tx);
    xo1=xo(:,1);
    xo2=xo(:,2);
    xo_out1 = xo1(ind_tx);
    xo_out2 = xo2(ind_tx);
    xo_out = [xo_out1 xo_out2];
   
else
    t_out = tw3;
    rid = rwd3;
    rii = rwi3;
    
    tx_out = tx;
    xo_out = xo;
end

% make sure that length of t_out is the same as tx
% if length(tw3)~=length(tx)
%     [~,ind_tw,ind_tx] = intersect(tw3,tx);
%     
%     t_out = NaN(size(tx));
%     ri = NaN(size(tx));
%     
%     t_out(ind) = tw3;
%     ri(ind) = rw3;
% else
%     t_out = tw3;
%     ri = rw3;
% end

% compute correlation co-efficient
% xo_out1=xo_out(:,1);
% 
% if (length(rid)/length(r1))<0.1
%     xc1=NaN;
% else
%     xc1=xcPH(rid,xo_out1,1);
% end
% 
% if size(xo_out,2)>1
%     xo_out2=xo_out(:,2);
% 
%     if (length(rii)/length(r2))<0.1
%         xc2=NaN;
%     else
%         xc2=xcPH(rii,xo_out2,1);
%     end
%     
% else
%     xc2=NaN;
% end

%overlapd = length(rid);
overlap = length(rii);

%% compute the output sequence

% % get the full sequence
if 0 % if you want to use the averaged sequence
    rid = rwd3;  
    rii = rwi3; 
    t_out = tw3;
else % if you have points with the same depth
    rid = r_alld; 
    rii = r_alli;
    t_out = t_all; 
end

% add extra values to the start of the section
[~, ind]=unique(start_t);
if length(ind)>1
    sr1 = start_r1(ind(1:end-1));
    rid = [sr1', rid];
    sr2 = start_r2(ind(1:end-1));
    rii = [sr2' rii];
    st = start_t(ind(1:end-1));
    t_out = [st' t_out];
end

[~, ind]=unique(append_t);
if length(ind)>1
    ar1 = append_r1(ind(2:end));
    rid = [rid ar1'];
    ar2 = append_r2(ind(2:end));
    rii = [rii ar2'];
    at = append_t(ind(2:end));
    t_out = [t_out at'];
end

rid = [start_r1' rid append_r1'];
rii = [start_r2' rii append_r2'];
t_out = [start_t' t_out append_t'];

targetdata=[tx,x1.',x2.'];
candidatedata=[t_out.',rid.',rii.'];

[output1 output2] = psvr2(targetdata, candidatedata);
xc = output2;

%% plot everything
if qplot,
    colormap hot; h=colormap; colormap(flipud(h));
    num=5;
    for ct=1:num-1;
        sp(ct,:)=(2:num)+(ct-1)*num;
        sp2(ct)=(ct-1)*num+1;
    end;
    subplot(num,num,sp(:));
    cla; hold on;
    imagesc(d); colorbar; caxis([0 60]) % plot the squared distance matrix, d
    axis tight; axis off;
    h=plot(1:length(x),1:length(x),'k--'); set(h,'linewidth',2);  % plot a line along the diagonal
    %h=plot(1:length(r),1:length(r),'k--'); set(h,'linewidth',2);  % plot a line along the diagonal
    pl=find( (w(:,1)>1) & (w(:,2)>1));   % plot the minimum distance path
    h=plot(w(pl,2),w(pl,1),'k'); set(h,'linewidth',2);
    
    set(gca,'ydir','normal');
    set(gca,'yaxisloc','right');
    set(gca,'xaxisloc','top');
    
    subplot(num,num,sp2'); cla; hold on;
    h=plot(xo,tx,'*r'); set(h,'linewidth',1);  % plot the "known" time series xo
    h=plot(rii,t_out,'*k'); set(h,'linewidth',1);  % plot the new warped time series
    h=plot(rid,t_out,'*k'); set(h,'linewidth',1);  % plot the new warped time series
%    h=plot(ri*std(xo)/nanstd(ri),tx,'r'); set(h,'linewidth',1);  % plot the new warped time series, ri, scaled to match x
    axis tight;
    h=suptitle(['cross-correlation: ',num2str(xc,2),', g: ',num2str(g,2),', edge: ',num2str(edge,2)]); %font(h,12);
    h=xlabel('d\delta^{13}C'); %font(h,12);
    h=ylabel('depth (m)'); %font(h,12);
    set(gca,'xaxisloc','top');
    %font(gca,12);
    
    subplot(num,num,num*(num-1)+2:num*num)
    h=plot(tr,ro1,'*k'); set(h,'linewidth',1); % plot the original non-warped time series, r0
    h=plot(tr,ro2,'*k'); set(h,'linewidth',1); % plot the original non-warped time series, r0
    axis tight;
    h=ylabel('d\delta^{13}C'); %font(h,12);
    h=xlabel('time (m)'); %font(h,12);
    %font(gca,12);
    set(gca,'yaxisloc','right');
    drawnow;
end;

