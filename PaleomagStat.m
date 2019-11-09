%% PALEOMAGSTAT

% Cedric J. Hagen, Oregon State University, 2-13-19

% A function for determining the statistical significance of resulting
% correlation solutions from dtw_paleomag.m. Run time can be substantial,
% depending on your chosen nitt value (number of synthetics). 

%% Parameters
g = 1.01 ; %Test 0.98, 0.99, 1.00, and 1.01
edge = 0.20;
nitt = 10000;  % number of synthetic time series generated

%% Load in Data
psv2269 = xlsread('Paleomag_Data/NNA_2269&2322_4dtw.xlsx', 1); %Target
psv2322 = xlsread('Paleomag_Data/NNA_2269&2322_4dtw.xlsx', 2); %Candidate

deep2269 = psv2269(:,1);
dec2269 = psv2269(:,2);
inc2269 = psv2269(:,3);
int2269 = inc2269*0+1;

deep2322 = psv2322(:,1);
dec2322 = psv2322(:,2);
inc2322 = psv2322(:,3);
int2322 = inc2322*0+1;

%% Call xyz to separate x y z components of inc and dec
[x2269 y2269 z2269] = xyz(dec2269, inc2269, int2269);
[x2322 y2322 z2322] = xyz(dec2322, inc2322, int2322);

%% Make synthetic x y z, combine with decinc for synthetic inc and dec
if 0
    %normalize
    x2269=(x2269-mean(x2269))/(std(x2269));
    y2269=(y2269-mean(y2269))/(std(y2269));
    z2269=(z2269-mean(z2269))/(std(z2269));

    x2322=(x2322-mean(x2322))/(std(x2322));
    y2322=(y2322-mean(y2322))/(std(y2322));
    z2322=(z2322-mean(z2322))/(std(z2322));

    %sort the depth values
    [deep2269,ind2269] = sort(deep2269);
    x2269 = x2269(ind2269);
    y2269 = y2269(ind2269);
    z2269 = z2269(ind2269);

    [deep2322,ind2322] = sort(deep2322);
    x2322 = x2322(ind2322);
    y2322 = y2322(ind2322);
    z2322 = z2322(ind2322);

    %only keep the unique values
    [deep2269,ind2269,~] = unique(deep2269);
    x2269 = x2269(ind2269);
    y2269 = y2269(ind2269);
    z2269 = z2269(ind2269);

    [deep2322,ind2322,~] = unique(deep2322);
    x2322 = x2322(ind2322);
    y2322 = y2322(ind2322);
    z2322 = z2322(ind2322);

    ntime_2269 = length(x2269);
    ntime_2322 = length(x2322);

    %detrend data
    x2269_detrend = detrend(x2269);
    x2269_trend = x2269 - x2269_detrend;
    y2269_detrend = detrend(y2269);
    y2269_trend = y2269 - y2269_detrend;
    z2269_detrend = detrend(z2269);
    z2269_trend = z2269 - z2269_detrend;

    x2322_detrend = detrend(x2322);
    x2322_trend = x2322 - x2322_detrend;
    y2322_detrend = detrend(y2322);
    y2322_trend = y2322 - y2322_detrend;
    z2322_detrend = detrend(z2322);
    z2322_trend = z2322 - z2322_detrend;

    %interpolate the detrended time series onto an evenly spaced grid
    tt_2269 = deep2269(1):mean(diff(deep2269)):deep2269(end);
    tt_2322 = deep2322(1):mean(diff(deep2322)):deep2322(end);

    yy_x2269 = interp1(deep2269,x2269_detrend,tt_2269,'pchip');
    yy_y2269 = interp1(deep2269,y2269_detrend,tt_2269,'pchip');
    yy_z2269 = interp1(deep2269,z2269_detrend,tt_2269,'pchip');

    yy_x2322 = interp1(deep2322,x2322_detrend,tt_2322,'pchip');
    yy_y2322 = interp1(deep2322,y2322_detrend,tt_2322,'pchip');
    yy_z2322 = interp1(deep2322,z2322_detrend,tt_2322,'pchip');

    %estimate the autocorrelation of the interpolated sample
    [xc_x2269,lags_x2269] = xcorr(yy_x2269,'coeff');
    [xc_y2269,lags_y2269] = xcorr(yy_y2269,'coeff');
    [xc_z2269,lags_z2269] = xcorr(yy_z2269,'coeff');

    [xc_x2322,lags_x2322] = xcorr(yy_x2322,'coeff');
    [xc_y2322,lags_y2322] = xcorr(yy_y2322,'coeff');
    [xc_z2322,lags_z2322] = xcorr(yy_z2322,'coeff');

    xc1_x2269 = xc_x2269(ntime_2269:end);
    xc1_y2269 = xc_y2269(ntime_2269:end);
    xc1_z2269 = xc_z2269(ntime_2269:end);

    xc1_x2322 = xc_x2322(ntime_2322:end);
    xc1_y2322 = xc_y2322(ntime_2322:end);
    xc1_z2322 = xc_z2322(ntime_2322:end);

    %perform the cholesky decomposition to make synthetics
    S_2269 = zeros(ntime_2269);
    S_2322 = zeros(ntime_2322);

    S_x2269 = S_2269 + diag(xc1_x2269(1)*ones(ntime_2269,1),0);
    S_y2269 = S_2269 + diag(xc1_y2269(1)*ones(ntime_2269,1),0);
    S_z2269 = S_2269 + diag(xc1_z2269(1)*ones(ntime_2269,1),0);

    S_x2322 = S_2322 + diag(xc1_x2322(1)*ones(ntime_2322,1),0);
    S_y2322 = S_2322 + diag(xc1_y2322(1)*ones(ntime_2322,1),0);
    S_z2322 = S_2322 + diag(xc1_z2322(1)*ones(ntime_2322,1),0);

    for ii=2:ntime_2269
        S_x2269 = S_x2269 + diag(xc1_x2269(ii)*ones(ntime_2269+1-ii,1),ii-1);
        S_x2269 = S_x2269 + diag(xc1_x2269(ii)*ones(ntime_2269+1-ii,1),-(ii-1));
        S_y2269 = S_y2269 + diag(xc1_y2269(ii)*ones(ntime_2269+1-ii,1),ii-1);
        S_y2269 = S_y2269 + diag(xc1_y2269(ii)*ones(ntime_2269+1-ii,1),-(ii-1));
        S_z2269 = S_z2269 + diag(xc1_z2269(ii)*ones(ntime_2269+1-ii,1),ii-1);
        S_z2269 = S_z2269 + diag(xc1_z2269(ii)*ones(ntime_2269+1-ii,1),-(ii-1));
    end

    for ii=2:ntime_2322
        S_x2322 = S_x2322 + diag(xc1_x2322(ii)*ones(ntime_2322+1-ii,1),ii-1);
        S_x2322 = S_x2322 + diag(xc1_x2322(ii)*ones(ntime_2322+1-ii,1),-(ii-1));
        S_y2322 = S_y2322 + diag(xc1_y2322(ii)*ones(ntime_2322+1-ii,1),ii-1);
        S_y2322 = S_y2322 + diag(xc1_y2322(ii)*ones(ntime_2322+1-ii,1),-(ii-1));
        S_z2322 = S_z2322 + diag(xc1_z2322(ii)*ones(ntime_2322+1-ii,1),ii-1);
        S_z2322 = S_z2322 + diag(xc1_z2322(ii)*ones(ntime_2322+1-ii,1),-(ii-1));
    end

    %cholesky decomposition
    [L_x2269,p_x2269] = chol(S_x2269,'lower');
    [L_y2269,p_y2269] = chol(S_y2269,'lower');
    [L_z2269,p_z2269] = chol(S_z2269,'lower');

    [L_x2322,p_x2322] = chol(S_x2322,'lower');
    [L_y2322,p_y2322] = chol(S_y2322,'lower');
    [L_z2322,p_z2322] = chol(S_z2322,'lower');

    samples_x2269 = NaN(ntime_2269,nitt);
    samples_y2269 = NaN(ntime_2269,nitt);
    samples_z2269 = NaN(ntime_2269,nitt);
    samples_x2322 = NaN(ntime_2322,nitt);
    samples_y2322 = NaN(ntime_2322,nitt);
    samples_z2322 = NaN(ntime_2322,nitt);

    samples_dec2269 = NaN(ntime_2269,nitt);
    samples_inc2269 = NaN(ntime_2269,nitt);
    samples_int2269 = NaN(ntime_2269,nitt);

    samples_dec2322 = NaN(ntime_2322,nitt);
    samples_inc2322 = NaN(ntime_2322,nitt);
    samples_int2322 = NaN(ntime_2322,nitt);

    for ii=1:nitt
        % random sample from sample distribution
        R_2269 = randn([ntime_2269 1]);
        R_2322 = randn([ntime_2322 1]);
    
        % multiple the cholesky decomposition by the random numbers
        samples_x2269(:,ii) = L_x2269*R_2269;
        samples_y2269(:,ii) = L_y2269*R_2269;
        samples_z2269(:,ii) = L_z2269*R_2269;
    
        samples_x2322(:,ii) = L_x2322*R_2322;
        samples_y2322(:,ii) = L_y2322*R_2322;
        samples_z2322(:,ii) = L_z2322*R_2322;
    
        for qq=1:length(L_z2269*R_2269)
            [samples_dec2269(:,ii), samples_inc2269(:,ii), samples_int2269(:,ii)] = decinc(L_x2269*R_2269, L_y2269*R_2269, L_z2269*R_2269);
        end
    
        for qq=1:length(L_z2322*R_2322)
            [samples_dec2322(:,ii), samples_inc2322(:,ii), samples_int2322(:,ii)] = decinc(L_x2322*R_2322, L_y2322*R_2322, L_z2322*R_2322);
        end
    
        % compute the autocorrelation of the synthetics
        [xc2_x2269,lags1_x2269] = xcorr(samples_x2269(:,ii),'coeff');
        [xc2_y2269,lags1_y2269] = xcorr(samples_y2269(:,ii),'coeff');
        [xc2_z2269,lags1_z2269] = xcorr(samples_z2269(:,ii),'coeff');

        [xc2_x2322,lags1_x2322] = xcorr(samples_x2322(:,ii),'coeff');
        [xc2_y2322,lags1_y2322] = xcorr(samples_y2322(:,ii),'coeff');
        [xc2_z2322,lags1_z2322] = xcorr(samples_z2322(:,ii),'coeff');
    end

    %save the synthetics
    save('NNA_x2269_synthetics_10000.mat','samples_x2269','tt_2269')
    save('NNA_y2269_synthetics_10000.mat','samples_y2269','tt_2269')
    save('NNA_z2269_synthetics_10000.mat','samples_z2269','tt_2269')

    save('NNA_x2322_synthetics_10000.mat','samples_x2322','tt_2322')
    save('NNA_y2322_synthetics_10000.mat','samples_y2322','tt_2322')
    save('NNA_z2322_synthetics_10000.mat','samples_z2322','tt_2322')

    save('NNA_dec2269_synthetics_10000.mat','samples_dec2269','tt_2269')
    save('NNA_inc2269_synthetics_10000.mat','samples_inc2269','tt_2269')
    save('NNA_int2269_synthetics_10000.mat','samples_int2269','tt_2269')

    save('NNA_dec2322_synthetics_10000.mat','samples_dec2322','tt_2322')
    save('NNA_inc2322_synthetics_10000.mat','samples_inc2322','tt_2322')
    save('NNA_int2322_synthetics_10000.mat','samples_int2322','tt_2322')
end

%% Determine correlation with psvr2

%normalize target
dx = deep2269;
dec2269=(dec2269-mean(dec2269))/(std(dec2269));
inc2269=(inc2269-mean(inc2269))/(std(inc2269));
tdi=[dec2269,inc2269];

%load the synthetic samples
load('NNA_dec2322_synthetics_10000.mat')
load('NNA_inc2322_synthetics_10000.mat')

%load in the estimate correlation coefficients
tmp = load('NNA_xc.mat');
n_overlap = tmp.n_overlap;

%normalize candidate
dr3 = deep2322;
dec2322=(dec2322-mean(dec2322))/(std(dec2322));
inc2322=(inc2322-mean(inc2322))/(std(inc2322));
cdi=[dec2322,inc2322];

%run dtw
xc_mean = NaN(length(g),length(edge));
xc_std = NaN(length(g),length(edge));
xc_95CI = NaN(length(g),length(edge));
xc_all = NaN(length(g),length(edge),nitt);
overlap_all = NaN(length(g),length(edge),nitt);

tt_2322 = tt_2322.';

cnt = 1;
for ii=1:length(g)
    for jj=1:length(edge)        
        xc_tmp = NaN(nitt,1);
        g_val = g(ii);
        edge_val = edge(jj);
        fprintf('\t Completion: ');
        showTimeToCompletion; startTime = tic;
        prog = parfor_progress(nitt);
        parfor kk=1:nitt
            rdec = samples_dec2322(:,kk);
            rinc = samples_inc2322(:,kk);
            r11=[rdec,rinc];
            [rid,rii,xc,d_out,Dist,w_scale,overlap]=dtw_paleomag(tdi,r11,dx,tt_2322,g_val,edge_val,[],[],[],0,1,0);
            overlap_all(ii,jj,kk) = overlap;
            xc_tmp(kk) = xc;
            pause(1);
            prog = parfor_progress;
            showTimeToCompletion(prog/100, [], [], startTime);
%             figure, hold on
%                plot(tdi(:,1),dx,'*','Color',[0.5 0.5 0.5])
%                plot(tdi(:,2),dx,'o','Color',[0.5 0.5 0.5])
%                plot(rid,d_out,'*b')
%                plot(rii,d_out,'ob')
%                ylim([-500 2500])
        end
        parfor_progress(0);
        xc_all(ii,jj,:) = xc_tmp;
        
        cnt = cnt+1;
    end
end

delete(gcp('nocreate'))

outputfile = ['Correlation_test/NNA/correlation_test_NNA_g',num2str(g),'.mat'];
save(outputfile,'xc_all','g','edge','overlap_all','nitt')

%% Use determine p-val
%count the number of synthetics that it finds alignment for
cnt = NaN(length(g),length(edge));
for ii=1:length(g)
    for jj=1:length(edge)
        cnt(ii,jj) = sum(~isnan(xc_all(ii,jj,:)));
    end
end

map = lbmap(10,'RedBlue');
map2 = lbmap(20,'RedBlue');

disp([num2str(cnt),' alignments out of ',num2str(nitt),' synthetics'])

figure  
   imagesc(cnt)
   colormap(flipud(map))
   colorbar
   title('Number of Alignments')
   set(gca,'YTick',[1:1:25],'YTickLabel',g)
   xlabel('edge parameter')
   ylabel('g value')

  print('-dpng','-r1100','Correlation_test/NNA/num_alignments.png')

%count the percent of alignments that are more than what we obtain
tmp = load('NNA_xc.mat');
g_ii = find(abs(tmp.g-g)==min(abs(tmp.g-g)));
xc = nanmean(tmp.NNA_xc(g_ii,:));

if 1
    % load in the estimate correlation coefficients
    tmp = load('NNA_xc.mat');

    percent = NaN(length(g),length(edge));
    for ii=1:length(g)
        for jj=1:length(edge)
            g_ii = find(tmp.g==g(ii));
            edge_jj = find(tmp.edge==edge(jj));
            
            if ~isnan(xc)
                ind = find((xc_all(ii,jj,:)>xc));
                percent(ii,jj) = sum(xc_all(ii,jj,:)>xc)/cnt(ii,jj);

            end
        end
    end

    disp(['p value = ',num2str(percent)])
    
    figure  
       imagesc(percent)
       colormap(flipud(map))
       colorbar
        caxis([0.5 1])
       title(['p-value: ',num2str(percent)])
       set(gca,'XTick',[1:2:15],'XTickLabel',edge(1):0.02:edge(end))
       set(gca,'YTick',[1:1:25],'YTickLabel',g)
       xlabel('edge parameter')
       ylabel('g value')
       
      outputimage = ['Correlation_test/NNA/p-value_g',num2str(g),'.png'];
      print('-dpng','-r1100',outputimage)
end

% plot the number of overlapping points in the original warped section
if 1
     figure  
        imagesc(tmp.n_overlap)
        colormap(flipud(map2))
        colorbar
        caxis([0 500])
        title('overlapping points')
        set(gca,'XTick',[1:2:15],'XTickLabel',edge(1):0.02:edge(end))
        set(gca,'YTick',[1:1:25],'YTickLabel',g)
        xlabel('edge parameter')
        ylabel('g value')

       outputimage = ['Correlation_test/NNA/n_overlap_g',num2str(g),'.png'];
       print('-dpng','-r1100',outputimage)
end

%plot the results
if 1
    figure
        tmp = squeeze(xc_all);
        counts = hist(tmp,20);
        hist(tmp,20)
        hold on
        %figure
        plot(xc*ones(max(counts)+1,1),0:max(counts),'r-','LineWidth',2)
        plot(prctile(squeeze(xc_all),95)*ones(max(counts)+1,1),0:max(counts),'k--','LineWidth',2)
        title(['g: ',num2str(g),', corr: ',num2str(xc),', p-val: ',num2str(percent)])
        outputimage = ['Correlation_test/NNA/hist_NNA_g',num2str(g),'.png'];
        print('-dpng','-r600',outputimage)
end

disp(g)
disp(xc)
disp(prctile(squeeze(xc_all),95))
