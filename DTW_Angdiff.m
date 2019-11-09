%% Load in data
psv2269 = xlsread('Paleomag_Data/NNA_2269&2322_4dtw.xlsx', 1); %Change to your file location
psvU1305 = xlsread('Paleomag_Data/NNA_U1305_4dtw.xlsx', 1); %Change to your file location

% target sequence
dx = psv2269(:,1);
tdi1 = psv2269(:,2);
tdi2 = psv2269(:,3);
tdi = [tdi1,tdi2];

% candidate sequence
dr3 = psvU1305(:,1);
cdi1 = psvU1305(:,2);
cdi2 = psvU1305(:,3);
cdi = [cdi1,cdi2];


%% Plots initial data
g = [0.98:0.01:1.01];
edge = [0.01:0.01:0.15];

NNA_U1305_xc = NaN(length(g),length(edge)); %Rename according to your data
ri_all = NaN(length(g),length(edge),length(tdi)+length(cdi));
dout_all = NaN(length(g),length(edge),length(tdi)+length(cdi));
n_overlap = NaN(length(g),length(edge)); % count the number of points that align

tx_fix=[];
tr_fix=[];
sd=[];

figure(1); clf; hold on; %Plots original data
    plot(dx,tdi,'k'); 
    plot(dr3,cdi,'r');
    saveas(figure(1),'Paleomag_Data/Output_Images/2269-U1305/2269-U1305_original.png'); %Change according to your desired directory
    close
    
    
%% Plots distance matrix and time warped data
xc_vals=NaN(length(g),length(edge));
pointfinder_depth=[];
new_depth_vals=[];
    for n=1:(length(g));
        for m=1:(length(edge));
            [rid,rii,xc,d_out,Dist,w_scale,overlap]=dtw_paleomag(tdi,cdi,dx,dr3,g(n),edge(m),tr_fix,tx_fix,sd,0,0,0);
            xc_vals(n,m)=xc;
            NNA_U1305_xc(n,m) = xc; %Rename according to your data
            rid_all(n,m,1:length(rid)) = rid;
            rii_all(n,m,1:length(rii)) = rii;
            dout_all(n,m,1:length(d_out)) = d_out;
            n_overlap(n,m) = overlap;
            figure; clf; hold on;

            plot(tdi,dx,'k.',rid,d_out,'b.',rii,d_out,'b.')
            title(['cross-corr.: ',num2str(xc,2),', g: ',num2str(g(n)),', edge: ',num2str(edge(m))]);
            saveas(figure(1),['Paleomag_Data/Output_Images/2269-U1305/2269-U1305_g',num2str(g(n)),'_edge',num2str(edge(m)),'.png']); %Change according to your desired directory

            outputdatavalues=[rii;rid;d_out];

            outputdata = fopen(['Paleomag_Data/Output_Data/2269-U1305/Depth_2269-U1305_g',num2str(g(n)),'_edge',num2str(edge(m)),'.txt'],'w'); %Change according to your desired directory
            fprintf(outputdata,'%6s %6s %12s\n','inc','dec','m');
            fprintf(outputdata,'%6.2f %6.2f %12.8f\n',outputdatavalues);
            fclose(outputdata);

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
saveas(figure(1),'Paleomag_Data/Output_Images/2269-U1305/2269-U1305_corr.png'); %Change according to your desired directory
close

[M,I] = max(NNA_U1305_xc(:)); %Rename according to your data
[I_row, I_col] = ind2sub(size(NNA_U1305_xc),I); %Rename according to your data
g_pref = g(I_row);
edge_pref = edge(I_col);

save('NNA_U1305_xc.mat','g','edge','NNA_U1305_xc','g_pref','edge_pref','n_overlap','rid_all','rii_all','dout_all') %Rename according to your data