function DTWPSV_Stat(psv_T, psv_C, name_T, name_C, g, edge, nitt, speed, gensyn)

% PALEOMAGSTAT

% Cedric J. Hagen, Oregon State University, 2-13-19

% A function for determining the statistical significance of resulting
% correlation solutions from dtw_paleomag.m. Run time can be substantial,
% depending on your chosen nitt value (number of synthetics). 

% BR edits 02-25-19
% -renamed variables to reflect target versus candidate PSV records
% -added comments to clarify code
% -erased call to xyz and decinc functions and added code inline
% -Changed method of normalization and detrending of xyz coponents
% -center data using vector mean and normalized variance to xyz st. dev.
% -changed data paths
% -made into a function

%% Parameters
% choose g value to test
g = g ; %(suggested = 0.98, 0.99, 1.00, and 1.01)
edge = edge; %(suggested = 0.01 to 0.15 by 0.01)
nitt = nitt;  % number of synthetic time series generated (suggested = 10000)

%% Load in Data
%load data files into psv_T and psv_C vectors
%name target and candidate PSV records
%T = Target; C = Candidate
psv_T = psv_T; %Target
psv_C = psv_C; %Candidate
name_T = name_T; %Name of Target PSV Curve
name_C = name_C; %Name of Candidate PSV Curve

%Break target and candidate data into depth, dec, and inc column vectors.  
%May require some editing by user.
deep_T = psv_T(:,1);
dec_T = psv_T(:,2);
inc_T = psv_T(:,3);
int_T = inc_T*0+1;

deep_C = psv_C(:,1);
dec_C = psv_C(:,2);
inc_C = psv_C(:,3);
int_C = inc_C*0+1;

%% Calculate x y z components of inc and dec
if speed ~= 1
x_T = int_T.*cosd(inc_T).*cosd(dec_T);
y_T = int_T.*cosd(inc_T).*sind(dec_T);
z_T = int_T.*sind(inc_T);
end

x_C = int_C.*cosd(inc_C).*cosd(dec_C);
y_C = int_C.*cosd(inc_C).*sind(dec_C);
z_C = int_C.*sind(inc_C);

%% Make synthetic x y z, combine with decinc for synthetic inc and dec
if gensyn == 1 %set to 0 if synthetics are already generated to save time

    %sort the depth values
    if speed ~= 1
    [deep_T,ind_T] = sort(deep_T);
    x_T = x_T(ind_T);
    y_T = y_T(ind_T);
    z_T = z_T(ind_T);
    end

    [deep_C,ind_C] = sort(deep_C);
    x_C = x_C(ind_C);
    y_C = y_C(ind_C);
    z_C = z_C(ind_C);

    %only keep the unique values
    if speed ~= 1
    [deep_T,ind_T,~] = unique(deep_T);
    x_T = x_T(ind_T);
    y_T = y_T(ind_T);
    z_T = z_T(ind_T);
    end

    [deep_C,ind_C,~] = unique(deep_C);
    x_C = x_C(ind_C);
    y_C = y_C(ind_C);
    z_C = z_C(ind_C);
    
    if speed ~= 1
    ntime_T = length(x_T);
    end
    ntime_C = length(x_C);
    
    
    %interpolate the time series onto an evenly spaced grid
    if speed ~= 1
    tt_T = deep_T(1):mean(diff(deep_T)):deep_T(end);
    end
    tt_C = deep_C(1):mean(diff(deep_C)):deep_C(end);
    
    if speed ~= 1
    yy_x_T = interp1(deep_T,x_T,tt_T,'pchip');
    yy_y_T = interp1(deep_T,y_T,tt_T,'pchip');
    yy_z_T = interp1(deep_T,z_T,tt_T,'pchip');
    end

    yy_x_C = interp1(deep_C,x_C,tt_C,'pchip');
    yy_y_C = interp1(deep_C,y_C,tt_C,'pchip');
    yy_z_C = interp1(deep_C,z_C,tt_C,'pchip');
    
    
    
    %Center, normalize and detrend data
    if speed ~= 1
    x_T_bar = mean(yy_x_T);
    x_T_std = std(yy_x_T);
    yy_x_T = detrend((yy_x_T - x_T_bar)/x_T_std);
    y_T_bar = mean(yy_y_T);
    y_T_std = std(yy_y_T);
    yy_y_T = detrend((yy_y_T - y_T_bar)/y_T_std);
    z_T_bar = mean(yy_z_T);
    z_T_std = std(yy_z_T);
    yy_z_T = detrend((yy_z_T - z_T_bar)/z_T_std);
    end
    
    x_C_bar = mean(yy_x_C);
    x_C_std = std(yy_x_C);
    yy_x_C = detrend((yy_x_C - x_C_bar)/x_C_std);
    y_C_bar = mean(yy_y_C);
    y_C_std = std(yy_y_C);
    yy_y_C = detrend((yy_y_C - y_C_bar)/y_C_std);
    z_C_bar = mean(yy_z_C);
    z_C_std = std(yy_z_C);
    yy_z_C = detrend((yy_z_C - z_C_bar)/z_C_std);

    %estimate the autocorrelation of the interpolated sample
    if speed ~= 1
    [xc_x_T,lags_x_T] = xcorr(yy_x_T,'coeff');
    [xc_y_T,lags_y_T] = xcorr(yy_y_T,'coeff');
    [xc_z_T,lags_z_T] = xcorr(yy_z_T,'coeff');
    end

    [xc_x_C,lags_xC] = xcorr(yy_x_C,'coeff');
    [xc_y_C,lags_y_C] = xcorr(yy_y_C,'coeff');
    [xc_z_C,lags_z_C] = xcorr(yy_z_C,'coeff');

    if speed ~= 1
    xc1_x_T = xc_x_T(ntime_T:end);
    xc1_y_T = xc_y_T(ntime_T:end);
    xc1_z_T = xc_z_T(ntime_T:end);
    end

    xc1_x_C = xc_x_C(ntime_C:end);
    xc1_y_C = xc_y_C(ntime_C:end);
    xc1_z_C = xc_z_C(ntime_C:end);

    %perform the cholesky decomposition to make synthetics
    if speed ~= 1
    S_T = zeros(ntime_T);
    end
    S_C = zeros(ntime_C);

    if speed ~= 1
    S_x_T = S_T + diag(xc1_x_T(1)*ones(ntime_T,1),0);
    S_y_T = S_T + diag(xc1_y_T(1)*ones(ntime_T,1),0);
    S_z_T = S_T + diag(xc1_z_T(1)*ones(ntime_T,1),0);
    end

    S_x_C = S_C + diag(xc1_x_C(1)*ones(ntime_C,1),0);
    S_y_C = S_C + diag(xc1_y_C(1)*ones(ntime_C,1),0);
    S_z_C = S_C + diag(xc1_z_C(1)*ones(ntime_C,1),0);

    if speed ~= 1
    for ii=2:ntime_T
        S_x_T = S_x_T + diag(xc1_x_T(ii)*ones(ntime_T+1-ii,1),ii-1);
        S_x_T = S_x_T + diag(xc1_x_T(ii)*ones(ntime_T+1-ii,1),-(ii-1));
        S_y_T = S_y_T + diag(xc1_y_T(ii)*ones(ntime_T+1-ii,1),ii-1);
        S_y_T = S_y_T + diag(xc1_y_T(ii)*ones(ntime_T+1-ii,1),-(ii-1));
        S_z_T = S_z_T + diag(xc1_z_T(ii)*ones(ntime_T+1-ii,1),ii-1);
        S_z_T = S_z_T + diag(xc1_z_T(ii)*ones(ntime_T+1-ii,1),-(ii-1));
    end
    end

    for ii=2:ntime_C
        S_x_C = S_x_C + diag(xc1_x_C(ii)*ones(ntime_C+1-ii,1),ii-1);
        S_x_C = S_x_C + diag(xc1_x_C(ii)*ones(ntime_C+1-ii,1),-(ii-1));
        S_y_C = S_y_C + diag(xc1_y_C(ii)*ones(ntime_C+1-ii,1),ii-1);
        S_y_C = S_y_C + diag(xc1_y_C(ii)*ones(ntime_C+1-ii,1),-(ii-1));
        S_z_C = S_z_C + diag(xc1_z_C(ii)*ones(ntime_C+1-ii,1),ii-1);
        S_z_C = S_z_C + diag(xc1_z_C(ii)*ones(ntime_C+1-ii,1),-(ii-1));
    end

    %cholesky decomposition
    if speed ~= 1
    [L_x_T,p_x_T] = chol(S_x_T,'lower');
    [L_y_T,p_y_T] = chol(S_y_T,'lower');
    [L_z_T,p_z_T] = chol(S_z_T,'lower');
    end

    [L_x_C,p_x_C] = chol(S_x_C,'lower');
    [L_y_C,p_y_C] = chol(S_y_C,'lower');
    [L_z_C,p_z_C] = chol(S_z_C,'lower');

    if speed ~= 1
    samples_x_T = NaN(ntime_T,nitt);
    samples_y_T = NaN(ntime_T,nitt);
    samples_z_T = NaN(ntime_T,nitt);
    end
    samples_x_C = NaN(ntime_C,nitt);
    samples_y_C = NaN(ntime_C,nitt);
    samples_z_C = NaN(ntime_C,nitt);

    if speed ~= 1
    samples_dec_T = NaN(ntime_T,nitt);
    samples_inc_T = NaN(ntime_T,nitt);
    samples_int_T = NaN(ntime_T,nitt);
    end

    samples_dec_C = NaN(ntime_C,nitt);
    samples_inc_C = NaN(ntime_C,nitt);
    samples_int_C = NaN(ntime_C,nitt);

    for ii=1:nitt
        % random sample from sample distribution
        if speed ~= 1
        R_T = randn([ntime_T 1]);
        end
        R_C = randn([ntime_C 1]);
  
        % multiple the cholesky decomposition by the random numbers
        % weight by standard deviation
        % add back vector mean
        if speed ~= 1
        samples_x_T(:,ii) = (L_x_T*R_T)/std(L_x_T*R_T) * x_T_std + x_T_bar;
        samples_y_T(:,ii) = (L_y_T*R_T)/std(L_y_T*R_T) * y_T_std + y_T_bar;
        samples_z_T(:,ii) = (L_z_T*R_T)/std(L_z_T*R_T) * z_T_std + z_T_bar;
        end
    
        samples_x_C(:,ii) = (L_x_C*R_C)/std(L_x_C*R_C) * x_C_std + x_C_bar;
        samples_y_C(:,ii) = (L_y_C*R_C)/std(L_y_C*R_C) * y_C_std + y_C_bar;
        samples_z_C(:,ii) = (L_z_C*R_C)/std(L_z_C*R_C) * z_C_std + z_C_bar;
        
        %Calculate Inc and Dec
        if speed ~= 1
        samples_int_T(:,ii) = sqrt(samples_x_T(:,ii).^2 ...
            +samples_y_T(:,ii).^2+samples_z_T(:,ii).^2);
        samples_dec_T(:,ii) = atan2d(samples_y_T(:,ii), samples_x_T(:,ii));
        samples_inc_T(:,ii) = asind(samples_z_T(:,ii)./samples_int_T(:, ii));
        end
        
        samples_int_C(:,ii) = sqrt(samples_x_C(:,ii).^2 ...
            +samples_y_C(:,ii).^2+samples_z_C(:,ii).^2);
        samples_dec_C(:,ii) = atan2d(samples_y_C(:,ii), samples_x_C(:,ii));
        samples_inc_C(:,ii) = asind(samples_z_C(:,ii)./samples_int_C(:, ii));
     
        % compute the autocorrelation of the synthetics
        if speed ~= 1
        [xc2_x_T,lags1_x_T] = xcorr(samples_x_T(:,ii),'coeff');
        [xc2_y_T,lags1_y_T] = xcorr(samples_y_T(:,ii),'coeff');
        [xc2_z_T,lags1_z_T] = xcorr(samples_z_T(:,ii),'coeff');
        end

        [xc2_x_C,lags1_x_C] = xcorr(samples_x_C(:,ii),'coeff');
        [xc2_y_C,lags1_y_C] = xcorr(samples_y_C(:,ii),'coeff');
        [xc2_z_C,lags1_z_C] = xcorr(samples_z_C(:,ii),'coeff');
    end

    %save the synthetics
    outputfolder_C = ['PSV_Data/' name_C '/Synthetics/'];
    if ~exist(outputfolder_C, 'dir')
        mkdir(outputfolder_C)
    end
    if speed ~= 1
    outputfolder_T = ['PSV_Data/' name_T '/Synthetics/'];
    if ~exist(outputfolder_T, 'dir')
        mkdir(outputfolder_T)
    end
    end
    
    if speed ~= 1
    save([outputfolder_T 'DTWPSV_x_' name_T '_synthetics_' num2str(nitt) '.mat'],...
        'samples_x_T', 'tt_T')
    save([outputfolder_T 'DTWPSV_y_' name_T '_synthetics_' num2str(nitt) '.mat'],...
        'samples_y_T', 'tt_T')
    save([outputfolder_T 'DTWPSV_z_' name_T '_synthetics_' num2str(nitt) '.mat'],...
        'samples_z_T', 'tt_T')
    end
    
    save([outputfolder_C 'DTWPSV_x_' name_C '_synthetics_' num2str(nitt) '.mat'],...
        'samples_x_C', 'tt_C')
    save([outputfolder_C 'DTWPSV_y_' name_C '_synthetics_' num2str(nitt) '.mat'],...
        'samples_y_C', 'tt_C')
    save([outputfolder_C 'DTWPSV_z_' name_C '_synthetics_' num2str(nitt) '.mat'],...
        'samples_z_C', 'tt_C')
    
    if speed ~= 1
    save([outputfolder_T 'DTWPSV_dec_' name_T '_synthetics_' num2str(nitt) '.mat'],...
        'samples_dec_T', 'tt_T')
    save([outputfolder_T 'DTWPSV_inc_' name_T '_synthetics_' num2str(nitt) '.mat'],...
        'samples_inc_T', 'tt_T')
    save([outputfolder_T 'DTWPSV_int_' name_T '_synthetics_' num2str(nitt) '.mat'],...
        'samples_int_T', 'tt_T')
    end
    
    save([outputfolder_C 'DTWPSV_dec_' name_C '_synthetics_' num2str(nitt) '.mat'],...
        'samples_dec_C', 'tt_C')
    save([outputfolder_C 'DTWPSV_inc_' name_C '_synthetics_' num2str(nitt) '.mat'],...
        'samples_inc_C', 'tt_C')
    save([outputfolder_C 'DTWPSV_int_' name_C '_synthetics_' num2str(nitt) '.mat'],...
        'samples_int_C', 'tt_C')
end

%% Determine correlation with psvr2

%target depths
dx = deep_T;
tdi=[dec_T,inc_T];
outputfolder_C = ['PSV_Data/' name_C '/Synthetics/'];
outputfolder_T = ['PSV_Data/' name_T '/Synthetics/'];

%load the synthetic samples for candidate curve
load([outputfolder_C 'DTWPSV_dec_' name_C '_synthetics_' num2str(nitt) '.mat']);
load([outputfolder_C 'DTWPSV_inc_' name_C '_synthetics_' num2str(nitt) '.mat']);

%load in the estimate correlation coefficients
% tmp = load('NNA_xc.mat');
outputdatafolder = ['PSV_Data/' name_C '/' name_T '-' name_C '/Output_Data/'];
tmp = load([outputdatafolder 'DTWPSV_' name_T '-' name_C '_xc.mat']);
n_overlap = tmp.n_overlap;

% Candidate Depths
dr3 = deep_C;
cdi=[dec_C,inc_C];

%run dtw
xc_mean = NaN(length(g),length(edge));
xc_std = NaN(length(g),length(edge));
xc_95CI = NaN(length(g),length(edge));
xc_all = NaN(length(g),length(edge),nitt);
overlap_all = NaN(length(g),length(edge),nitt);

tt_C = tt_C.';

cnt = 1;
for ii=1:length(g)
    for jj=1:length(edge)        
        xc_tmp = NaN(nitt,1);
        g_val = g(ii);
        edge_val = edge(jj);
%         fprintf('\t Completion: ');
%         showTimeToCompletion; startTime = tic;
%         prog = parfor_progress(nitt);
        parfor kk=1:nitt
            rdec = samples_dec_C(:,kk);
            rinc = samples_inc_C(:,kk);
            r11=[rdec,rinc];
            [rid,rii,xc,d_out,Dist,w_scale,overlap]=dtw_paleomag(tdi,r11,dx,tt_C,g_val,edge_val,0,1,0);
            overlap_all(ii,jj,kk) = overlap;
            xc_tmp(kk) = xc;
%             pause(1);
%             prog = parfor_progress;
%             showTimeToCompletion(prog/100, [], [], startTime);
%             figure, hold on
%                plot(tdi(:,1),dx,'*','Color',[0.5 0.5 0.5])
%                plot(tdi(:,2),dx,'o','Color',[0.5 0.5 0.5])
%                plot(rid,d_out,'*b')
%                plot(rii,d_out,'ob')
%                ylim([-500 2500])
        end
%         parfor_progress(0);
        xc_all(ii,jj,:) = xc_tmp;
        
        cnt = cnt+1;
    end
end

delete(gcp('nocreate'))

outputcorrfolder = ['PSV_Data/' name_C '/' name_T '-' name_C '/Correlation_Test/'];
if ~exist(outputcorrfolder, 'dir')
   mkdir(outputcorrfolder)
end

outputfile = [outputcorrfolder 'correlation_test_' name_T '-' name_C '_g',num2str(g), '_' num2str(nitt) '.mat'];
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
% map = 'jet';
% map2 = 'jet';

disp([num2str(cnt),' alignments out of ',num2str(nitt),' synthetics'])

figure  
   imagesc(cnt)
   colormap(flipud(map))
   colorbar
   title('Number of Alignments')
   set(gca,'YTick',[1:1:25],'YTickLabel',g)
   xlabel('edge parameter')
   ylabel('g value')

  print('-dpng','-r1100', [outputcorrfolder 'num_alignments.png'])
  close
%%
%count the percent of alignments that are more than what we obtain
tmp = load([outputdatafolder 'DTWPSV_' name_T '-' name_C '_xc.mat']);
g_ii = find(abs(tmp.g-g)==min(abs(tmp.g-g)));
xc = nanmean(tmp.DTWPSV_C_xc(g_ii,:));

if 1
    % load in the estimate correlation coefficients
    tmp = load([outputdatafolder 'DTWPSV_' name_T '-' name_C '_xc.mat']);

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
       
      outputimage = [outputcorrfolder 'p-value_g',num2str(g),'.png'];
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

       outputimage = [outputcorrfolder 'n_overlap_g',num2str(g),'.png'];
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
        outputimage = [outputcorrfolder 'hist_NNA_g',num2str(g),'.png'];
        print('-dpng','-r600',outputimage)
end

disp(g)
disp(xc)
disp(prctile(squeeze(xc_all),95))
