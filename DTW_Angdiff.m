%% Load in data
psv2269 = xlsread('Paleomag_Data/NNA_2269&2322_4dtw.xlsx', 1);
psvU1305 = xlsread('Paleomag_Data/NNA_U1305_4dtw.xlsx', 1);

% target sequence
% i = ~isnan(psv(:, 31)) & psv(:, 23) < 21;
% dx = psv(i, 23);
% tdi = psv(i, 30:31);
% tdi(:, 2) = tdi(:, 2) - 8;
dx = psv2269(:,1);
tdi1 = psv2269(:,2);
tdi2 = psv2269(:,3);
tdi = [tdi1,tdi2];

% candidate sequence
% i = ~isnan(psv(:, 20));
% dr3 = psv(i, 12);
% cdi = psv(i, 19:20);
dr3 = psvU1305(:,1);
cdi1 = psvU1305(:,2);
cdi2 = psvU1305(:,3);
cdi = [cdi1,cdi2];


%% Plots initial data
g = [0.98:0.01:1.01];
edge = [0.01:0.01:0.15];

NNA_U1305_xc = NaN(length(g),length(edge));
ri_all = NaN(length(g),length(edge),length(tdi)+length(cdi));
dout_all = NaN(length(g),length(edge),length(tdi)+length(cdi));
n_overlap = NaN(length(g),length(edge)); % count the number of points that align

tx_fix=[];
tr_fix=[];
sd=[];

figure(1); clf; hold on; %Plots original data
    plot(dx,tdi,'k'); 
    plot(dr3,cdi,'r');
    saveas(figure(1),'Paleomag_Data/Output_Images/2269-U1305/2269-U1305_original.png');
    close
%% Plots distance matrix and time warped data
xc_vals=NaN(length(g),length(edge));
pointfinder_depth=[];
new_depth_vals=[];
    for n=1:(length(g));
        for m=1:(length(edge));
            [rid,rii,xc,d_out,Dist,w_scale,overlap]=dtw_paleomag(tdi,cdi,dx,dr3,g(n),edge(m),tr_fix,tx_fix,sd,0,0,0);
            xc_vals(n,m)=xc;
            NNA_U1305_xc(n,m) = xc;
            rid_all(n,m,1:length(rid)) = rid;
            rii_all(n,m,1:length(rii)) = rii;
            dout_all(n,m,1:length(d_out)) = d_out;
            n_overlap(n,m) = overlap;
            figure; clf; hold on;

% This section allows for one to find how specific points have move
% according to the fits. Place original depth constraints for d_out below
% (1). pointfinder_depth then produces a list of depth values. Find the average depth for each non-unique fit, and
% then match that to the corresponding age for the candidate section. Thus,
% each non-unique solution will (may) produce a different value for the
% point of interest. 
%             a=397;
%             b=396;
%             pointfinder_deep=find(dr3>b & dr3<=a); % (1)
%             if isempty(pointfinder_deep)==1;
%                 pointfinder_depth=[pointfinder_depth,0];
%             else
%                 count=0;
%                 sum=0;
%                 for q=1:(length(pointfinder_deep));
%                     sum=sum+pointfinder_deep(q);
%                     count=count+1;
%                 end
%                 truenewpoint=[];
%                 average_pointfinder=sum/count;
%                 indy=average_pointfinder;
%                 rindy=round(indy);
%                 x_val=dr3(rindy);
%                 newpoint=find(ri==x_val);
%                 newdepth=d_out(newpoint);
%                 if length(newdepth)>1;
%                     newdepth=min(newdepth);
%                 end
%                 new_depth_vals=[new_depth_vals,newdepth];
%                 pointfinder_depth=[pointfinder_depth,average_pointfinder]; %(2)
%             end
            %new_depth_vals
            plot(tdi,dx,'k.',rid,d_out,'b.',rii,d_out,'b.')
            title(['cross-corr.: ',num2str(xc,2),', g: ',num2str(g(n)),', edge: ',num2str(edge(m))]);
            saveas(figure(1),['Paleomag_Data/Output_Images/2269-U1305/2269-U1305_g',num2str(g(n)),'_edge',num2str(edge(m)),'.png']);
%             d_age=[];
%             for j=1:length(d_out);
%                 d_out_pt=d_out(j);
%                 agepoint=find(depthkey==d_out_pt);
%                 d_age_pt=agekey(agepoint);
%                 if d_out_pt>max(depthkey);
%                     d_age_pt=-999;
%                 end
%                 if length(d_age_pt)>1;
%                     d_age_pt=mean(d_age_pt);
%                 end
%                 if length(agepoint)<1;
%                     d_age_pt=-999;
%                 end
%                 d_age=[d_age,d_age_pt];
%             end
%             outputdatavalues1=[ri;d_age];
            outputdatavalues=[rii;rid;d_out];
%            outputdatavalues2_i=[rii;d_out];
%             outputdata1 = fopen(['Paleomag_Data/Output_Data/Utah-Idaho/Age_Utah-Idaho_g',num2str(g(n)),'_edge',num2str(edge(m)),'.txt'],'w');
%             fprintf(outputdata1,'%6s %12s\n','d13c','age');
%             fprintf(outputdata1,'%6.2f %12.8f\n',outputdatavalues1);
%             fclose(outputdata1);
            outputdata = fopen(['Paleomag_Data/Output_Data/2269-U1305/Depth_2269-U1305_g',num2str(g(n)),'_edge',num2str(edge(m)),'.txt'],'w');
            fprintf(outputdata,'%6s %6s %12s\n','inc','dec','m');
            fprintf(outputdata,'%6.2f %6.2f %12.8f\n',outputdatavalues);
            fclose(outputdata);
%             outputdata2_i = fopen(['Paleomag_Data/Output_Data/2269-2322/Depth_2269-2322_g',num2str(g(n)),'_edge',num2str(edge(m)),'.txt'],'w');
%             fprintf(outputdata2_i,'%6s %12s\n','d13c','m');
%             fprintf(outputdata2_i,'%6.2f %12.8f\n',outputdatavalues2_i);
%             fclose(outputdata2_i);
            close
        end
    end

map=lbmap(20,'RedBlue');
figure;
imagesc(xc_vals);
colormap(flipud(map))
colorbar;
ylabel('g value');
xlabel('edge parameter');
set(gca,'XTick',[1:2:15],'XTickLabel',edge(1):0.02:edge(end));
set(gca,'YTick',[1:1:4],'YTickLabel',g);
saveas(figure(1),'Paleomag_Data/Output_Images/2269-U1305/2269-U1305_corr.png');
close

[M,I] = max(NNA_U1305_xc(:));
[I_row, I_col] = ind2sub(size(NNA_U1305_xc),I);
g_pref = g(I_row);
edge_pref = edge(I_col);

save('NNA_U1305_xc.mat','g','edge','NNA_U1305_xc','g_pref','edge_pref','n_overlap','rid_all','rii_all','dout_all')