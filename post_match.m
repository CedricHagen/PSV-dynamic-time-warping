%Load in original candidate data
psv2322 = xlsread('Paleomag_Data/NNA_2269&2322_4dtw.xlsx', 2);
dr3 = psv2322(:,1);
origdec = psv2322(:,2);
originc = psv2322(:,3);

cnt=1;
progressbar
for i=0.98:0.01:1.01
    pause(0.01)
    progressbar(cnt/4)
    for j=0.01:0.01:0.15
        
        %Load in warped candidate resulting data
        warped = importdata(['Paleomag_Data/Output_Data/2269-2322/Depth_2269-2322_g',num2str(i),'_edge',num2str(j),'.txt']);
        inc = warped.data(:,1);
        dec = warped.data(:,2);
        newdepth = warped.data(:,3);
        
        origdepth=[];
        
        for n=1:length(newdepth)
            a=find(round(origdec,1)==dec(n));
            b=find(round(originc,1)==inc(n));
            c=[];
            d=[];
            for m=1:length(a)
                c=[c,find(b==a(m))];
            end
            for p=1:length(b)
                d=[d,find(a==b(p))];
            end
            if c>0 & d>0;
                if length(c)>1
                    c=c(1);
                end
                origdepth=[origdepth,dr3(b(c))];
            else
                origdepth=[origdepth,-9999];
                disp(['Warning! Problem at: ', num2str(n)])
            end
        end
        
        for w=1:length(origdepth)
            if w>1
                if origdepth(w)>=origdepth(w-1)
                elseif w < max(length(origdepth))
                    newind=find(dr3==origdepth(w-1));
                    if newind == length(dr3)
                        origdepth(w)=dr3(newind);
                    else
                        origdepth(w)=dr3(newind+1);
                    end
                end
            end
        end
        
        for t=1:length(origdepth)
            if t>1
                if origdepth(t)>=origdepth(t-1)
                else
                    disp(['t = ',num2str(t),', g = ',num2str(i),', edge = ',num2str(j)])
                end
            end
        end
        origdepth=origdepth';
        outputdatavalues=[inc';dec';newdepth';origdepth'];
        
        outputdata = fopen(['Paleomag_Data/Output_Data/2269-2322/Depth_2269-2322_g',num2str(i),'_edge',num2str(j),'_postmatch.txt'],'w');
        fprintf(outputdata,'%6s %6s %6s %6s\n','inc','dec','m', 'orig m');
        fprintf(outputdata,'%6.2f %6.2f %6.2f %6.2f\n',outputdatavalues);
        fclose(outputdata);
    end
    cnt=cnt+1;
end