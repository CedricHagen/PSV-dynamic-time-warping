function DTWPSV_main(psv_T, psv_C, name_T, name_C, g, edge)

%% organize data
% target sequence
i = ~isnan(psv_T(:, 1));
dx_T = psv_T(i,1);
di_T = psv_T(i, [2 3]);

% candidate sequence
i = ~isnan(psv_C(:, 1));
dx_C = psv_C(i,1);
di_C = psv_C(i, [2 3]);

%% Plots initial data

DTWPSV_C_xc = NaN(length(g),length(edge));
ri_all = NaN(length(g),length(edge),length(di_T)+length(di_C));
dout_all = NaN(length(g),length(edge),length(di_T)+length(di_C));
n_overlap = NaN(length(g),length(edge)); % count the number of points that align

tx_fix=[];
tr_fix=[];
tx_fix_sd = [];
tr_fix_sd = [];
sd=[];
acc_wei = [];

% make directory if it doesn't exist already
outputfolder = ['PSV_Data/' name_C '/' name_T '-' name_C '/Output_Images/'];
if ~exist(outputfolder, 'dir')
    mkdir(outputfolder)
end

figure(1); clf; %Plots original data
    subplot(2, 1, 1)
    hold on 
    plot(dx_T,di_T(:, 1),'k'); 
    plot(dx_C,di_C(:, 1),'r');
    mm = [di_T(:, 1); di_C(:, 1)];
    ylim([min(mm) max(mm)])
    mm = [dx_T; dx_C];
    xlim([min(mm) max(mm)])
    ylabel('Dec (^o)')
    title([name_T '-' name_C '-Original'])
    
    subplot(2, 1, 2)
    hold on
    plot(dx_T,di_T(:, 2),'k'); 
    plot(dx_C,di_C(:, 2),'r');
    mm = [di_T(:, 2); di_C(:, 2)];
    ylim([min(mm) max(mm)])
    mm = [dx_T; dx_C];
    xlim([min(mm) max(mm)])
    ylabel('Inc (^o)')
    xlabel('Depth')
 
    saveas(figure(1), [outputfolder name_T '-' name_C '_original.png']);
    close
%% Plots distance matrix and time warped data
xc_vals=NaN(length(g),length(edge));
pointfinder_depth=[];
new_depth_vals=[];
    for n=1:(length(g));
        for m=1:(length(edge));
            [rid,rii,xc,d_out,Dist,w_scale,overlap]=dtw_paleomag(di_T,di_C,dx_T,dx_C,g(n),edge(m),0,0,0);
            xc_vals(n,m)=xc;
            DTWPSV_C_xc(n,m) = xc;
            rid_all(n,m,1:length(rid)) = rid;
            rii_all(n,m,1:length(rii)) = rii;
            dout_all(n,m,1:length(d_out)) = d_out;
            n_overlap(n,m) = overlap;

            % plot results
            figure(1); clf; 
            subplot(2, 1, 1)
            hold on
            plot(dx_T, di_T(:, 1), 'k')
            plot(d_out, rid, 'r')
            mm = [di_T(:, 2); rid'];
            ylim([min(mm) max(mm)])
            mm = [dx_T; d_out'];
            xlim([min(mm) max(mm)])
            ylabel('Dec (^o)')
            title([name_T '-' name_C ' cross-corr.: ',num2str(xc,2),', g: ',num2str(g(n)),', edge: ',num2str(edge(m))]);
            
            subplot(2, 1, 2)
            hold on
            plot(dx_T, di_T(:, 2), 'k')
            plot(d_out, rii, 'r')
            mm = [di_T(:, 2); rii'];
            ylim([min(mm) max(mm)])
            mm = [dx_T; d_out'];
            xlim([min(mm) max(mm)])
            ylabel('Inc (^o)')
            xlabel('Depth')
            
            saveas(figure(1),[outputfolder name_T '-' name_C '_g',num2str(g(n)),'_edge',num2str(edge(m)), '_accW', num2str(acc_wei) '.png']);
   
            
            % make depth-depth conversion
            d_outp = nan(size(d_out));
            for ii = 1:length(d_out)
                idx = find(di_C(:, 1) == rid(ii) & di_C(:, 2) == rii(ii));
                if ii == 1
                    d_outp(ii) = dx_C(idx(1));
                else
                    tdep = dx_C(idx);
                    idx2 = find(tdep >= d_outp(ii-1));
                    if ~isempty(idx2)
                        d_outp(ii) = tdep(idx2(1));
                    end
                end
            end
        
            uo = unique(d_outp);
            uo = uo(~isnan(uo));
            uw = nan(size(uo));
            for ii = 1:length(uo)
                idx = d_outp == uo(ii);
                uw(ii) = mean(d_out(idx));
            end
            
            uw1 = unique(uw);
            uw1 = uw1(~isnan(uw1));
            uo1 = nan(size(uw1));
            for ii = 1:length(uw1);
                idx = uw == uw1(ii);
                uo1(ii) = mean(uo(idx));
            end
                
         
            % create data output files
            outputdatafolder = ['PSV_Data/' name_C '/' name_T '-' name_C '/Output_Data/'];
            if ~exist(outputdatafolder, 'dir')
                mkdir(outputdatafolder)
            end
            
            outputdatavalues=[rii;rid;d_out;d_outp];
            outputdata = fopen([outputdatafolder name_T '-' name_C '_g',num2str(g(n)),'_edge',num2str(edge(m)), '_accW', num2str(acc_wei) '.txt'],'w');
            fprintf(outputdata,'%6s\t %6s\t %12s\t %12s\n','inc','dec','warped_depth', 'orrig_depth');
            fprintf(outputdata,'%6.2f\t %6.2f\t %12.2f\t %12.2f\n',outputdatavalues);
            fclose(outputdata);
            
            outputdatavaluesdepth = [uo1; uw1];
            outputdatadepths = fopen([outputdatafolder name_T '-' name_C '_depthconv_g',num2str(g(n)),'_edge',num2str(edge(m)), '_accW', num2str(acc_wei) '.txt'],'w');
            fprintf(outputdatadepths,'%12s\t %12s\n','orrig_depth', 'warped_depth');
            fprintf(outputdatadepths,'%12.2f\t %12.2f\n',outputdatavaluesdepth);
            fclose(outputdatadepths);
            close
           
        end
    end

map=lbmap(20,'RedBlue');
% map = 'jet';
figure(1);
imagesc(xc_vals);
colormap(flipud(map))
colorbar;
ylabel('g value');
xlabel('edge parameter');
set(gca,'XTick',[1:2:15],'XTickLabel',edge(1):0.02:edge(end));
set(gca,'YTick',[1:1:4],'YTickLabel',g);
saveas(figure(1),[outputfolder name_T '-' name_C  '_accW', num2str(acc_wei) '_corr.png']);

close

[M,I] = max(DTWPSV_C_xc(:));
[I_row, I_col] = ind2sub(size(DTWPSV_C_xc),I);
g_pref = g(I_row);
edge_pref = edge(I_col);

save([outputdatafolder 'DTWPSV_' name_T '-' name_C '_xc.mat'],'g','edge','DTWPSV_C_xc','g_pref','edge_pref','n_overlap','rid_all','rii_all','dout_all')
