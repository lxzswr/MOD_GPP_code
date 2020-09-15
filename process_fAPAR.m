% retrieve the coordinates of pixels for each tower.


% PART 1:
% input files are 3km * 3km ~500m pixels, therefore, 13 rows * 13 cols
% but I only know the coordinates of left bottom pixels.

% therefore, need to get a matrix that has coordinates
% pixel ID, x (lb,lu,ru,rb), y (lb,lu,ru,rb)

% therefore, we can plot out each pixels as polygons

trig1=0; % a trig to run/unrun this routine

if (trig1)
    
    % directory of download fAPAR file from MODIStool R code
    dl_fAPAR_dir = '/Volumes/Elements/project14_AFF/MODIS_cutouts/downloadedFpar/';
    out_dir = '/Volumes/Elements/project14_AFF/MODIS_cutouts/sites_sin_coords/';
    
    % get the list of sites to reformat
    d = dir(dl_fAPAR_dir);
    isub = [d(:).isdir]; % returns logical vector
    sites = {d(isub).name}';
    % remove the non-sites '.' and '..'
    sites(ismember(sites,{'.','..'})) = [];
    
    
    for i=1:length(sites) % for each site
        
        cSite=sites(i);
        disp(cSite);

        home=strcat(dl_fAPAR_dir, cSite);

        filename1=strcat(home,'/*.csv');

        filename=dir(filename1{:});

        % load the first yearly data file, extract data, and concatenate
            tmp=strcat(home,'/',filename(1).name);
            inData=importdata(tmp{:});

            % extract the data from text
            info=inData.textdata(1:15,1); 
        
        pixel_coords=[]; % initiate the pixel_coords
        % it should include
        % pixel ID, x (lb,lu,ru,rb), y (lb,lu,ru,rb)
    
        % read in the basic information
        tmp=char(info(3,1));
        cell_size = str2num(tmp(13:end)); % cell size, about 463 m


        % left bottom corner of the image
        tmp=char(info(1,1));
        basic_x = str2num(tmp(14:end));
        tmp=char(info(2,1));
        basic_y =str2num(tmp(14:end));
        
        tmp=char(info(4,1));
        n_row=str2num(tmp(10:end));
        tmp=char(info(5,1));
        n_col=str2num(tmp(10:end));
        
        tmp=char(info(8,1));
        scale=str2num(tmp(10:end));
        
        tmp=char(info(9,1));
        site_lat=str2num(tmp(13:end));
        tmp=char(info(10,1));
        site_long=str2num(tmp(14:end));
        
        for r=1:n_row % for each row 
            
            for c=1:n_col % for each col
                
                pixel_ID = c + (r-1)*n_col;
                
                pixel_coords(pixel_ID,1)=pixel_ID;
                
                % x coordinates
                pixel_coords(pixel_ID,2) = basic_x + (c-1)*cell_size;
                pixel_coords(pixel_ID,3) = basic_x + (c-1)*cell_size;
                pixel_coords(pixel_ID,4) = basic_x + (c)*cell_size;
                pixel_coords(pixel_ID,5) = basic_x + (c)*cell_size;
                
                
                % y coordinates
                pixel_coords(pixel_ID,6) = basic_y + (n_row-r)*cell_size;
                pixel_coords(pixel_ID,7) = basic_y + (n_row-r+1)*cell_size;
                pixel_coords(pixel_ID,8) = basic_y + (n_row-r)*cell_size;
                pixel_coords(pixel_ID,9) = basic_y + (n_row-r+1)*cell_size;
                
            end
            
        end
    
        save([out_dir,char(cSite),'_sin.mat'],'pixel_coords');
    
    end
    
    
    
end

% test
% figure;
% scatter(pixel_coords(:,2),pixel_coords(:,6));



% PART 2: use R code to convert these coordinates to lat and long...
% name of the R routine: convertfromsin2ll.r
% location: folder sites_ll
% output file: pixel_ll
% pixel ID, lat (lb,lu,ru,rb), longs (lb,lu,ru,rb)

% the following function is just a test of the coordinates converted
% adding the site coordinates to each file

trig2=0;

if (trig2)
    
    % read every sites coordinates
    filename2=strcat('/Volumes/Elements/project14_AFF/MODIS_cutouts/','All_BASE-site_list.csv');
    
    M=readtable(filename2);
    
    name_list=table2array(M(:,1));
    latlong_list=table2array(M(:,3:4));
    
    ll_dir = '/Volumes/Elements/project14_AFF/MODIS_cutouts/sites_ll/';
    ll_files=dir([ll_dir,'*mat']);
    
    figure;
    
    for i=1:length(ll_files)
        
        site_name=ll_files(i).name(1:6);
        ind=strcmp(name_list,site_name);
        site_latlong=latlong_list(ind,:);
        
        %load corresponding mat file
        load(strcat(ll_dir,ll_files(i).name));
        
        
        subplot(5,10,i);
        hold on;
        scatter(pixel_ll(:,6),pixel_ll(:,2));
        scatter(site_latlong(2),site_latlong(1));
        
        
    end

end




% PART 3: use matlab to calculate the weight of each pixel for a tower
% (monthly)
% output format:
% ID, weight for 1st month, 2nd month,............

trig3=0;

warning off;

%figure;

if (trig3)
    
    ll_dir = '/Volumes/Elements/project14_AFF/MODIS_cutouts/sites_ll/';
    foot_dir = '/Volumes/Elements/project14_AFF/AMF-footprint-analysis/footprint-climatology/';
    out_dir = '/Volumes/Elements/project14_AFF/MODIS_cutouts/sites_weight/';
    
    ll_files=dir([ll_dir,'*mat']);
    
    for i=1:length(ll_files) % do it site by site, some site may not have footprint information
        
        %initiate output file
        ft_weight=[];
        
        site_name=ll_files(i).name(1:6);
        load(strcat(ll_dir,ll_files(i).name)); % load pixel_ll file for the site
        
        ft_files=dir([foot_dir,'*',char(site_name),'*.csv']);
        
        for j=1:length(ft_files) % for every month data
            
            M=readtable(strcat(foot_dir,ft_files(j).name));
            tmp=table2array(M);
            
            % get the polygon for footprint, use ft size and model
            ind=tmp(:,3)==0.8 & tmp(:,4)==1;
            
            if sum(ind)==0 % if ft model 1 does not work, use the 2nd model
                ind=tmp(:,3)==0.8 & tmp(:,4)==2;
            end
            
            ft_poly=tmp(ind,1:2);
            
            res_area=[]; % store the intersected area
            
            for m=1:length(pixelll) % for every pixel
                
                poly1=polyshape(ft_poly(:,1),ft_poly(:,2)); %long and lat
                
                % the x coordinates must follow the sequence of its
                % rectangle...
                poly2=polyshape(pixelll(m,[8,6,7,9])',pixelll(m,[4,2,3,5])'); % long and lat
                
                [polyout,~,~] = intersect(poly1,poly2);
                                
                res_area(m,1)=area(polyout);
                
%                 pg=plot(poly2);
%                 pg.FaceColor = [0.8,0.8,0.8];
%                 hold on;
            end
            
%             pg2=plot(poly1);
%             pg2.FaceColor = [0.8,0,0];
            % get the weight for every pixel for this month
            ft_weight(:,j)=res_area(:,1)./sum(res_area);
            
        end
        
        save([out_dir,char(site_name),'_weight.mat'],'ft_weight');
        
        disp(i);
        
    end
    
    
end



% PART 4: based on PART 3 results, to get the monthly APAR for every sites
% output: site-specific fAPAR file
% rows: months, cols: ft_fAPAR, 1*1,2*2,3*3,4*4,5*5

trig4=0;

if (trig4)

    fAPAR_dir='/Volumes/Elements/project14_AFF/MODIS_cutouts/formattedFpar/';
    fAPAR_files=dir([fAPAR_dir,'*.mat']);
    
    weight_dir = '/Volumes/Elements/project14_AFF/MODIS_cutouts/sites_weight/';
    
    out_dir = '/Volumes/Elements/project14_AFF/MODIS_cutouts/weight_fAPAR/';
    
    % this create other weight options
    org = zeros(13,13);
    
        % 1*1
        org(7,7)=1;
        weight_1=reshape(org,[],1);
        
        % 3*3
        org(6:8,6:8)=1/9;
        weight_3=reshape(org,[],1);
        
        % 5*5
        org(5:9,5:9)=1/25;
        weight_5=reshape(org,[],1);

    
    %  the end of create other weight options

    % for every site
    for i=1:length(fAPAR_files)
        
        % initiate the output file
        % doy, ft_fAPAR, 1*1, 3*3, 5*5
        ft_fAPAR=[];
        
        site_name=fAPAR_files(i).name(1:6); 

        % read in 8-day fAPAR value
        load([fAPAR_dir,fAPAR_files(i).name]); %FparData
        
        % interpolate 8 day fAPAR to monthly
        month_d=[1,32,60,91,121,152,182,213,244,274,305,335];
        m_fAPAR=nan(12,171);
        m_fAPAR(:,1)=FparData(1,1); %year
        m_fAPAR(:,2)=month_d';
        
        for dd=1:169 % do it for every pixel
            
            m_fAPAR(:,dd+2) = interp1(FparData(:,2),FparData(:,dd+2),m_fAPAR(:,2));
            
        end
 

        % read in fAPAR weight
        load(strcat(weight_dir,char(site_name),'_weight.mat')); %ft_weight

        
        if ~isempty(ft_weight)
            
            for mm=1:12 % for every month
            
                ft_fAPAR(mm,1)=nanmean(m_fAPAR(mm,3:end).*ft_weight(:,mm)',2);
            
            end



            % get other weighted fAPAR 

            w1_fAPAR=nanmean(m_fAPAR(:,3:end).*weight_1',2);
            w3_fAPAR=nanmean(m_fAPAR(:,3:end).*weight_3',2);
            w5_fAPAR=nanmean(m_fAPAR(:,3:end).*weight_5',2);

            ft_fAPAR(:,2)=w1_fAPAR;
            ft_fAPAR(:,3)=w3_fAPAR;
            ft_fAPAR(:,4)=w5_fAPAR;

           % interpolate them to months

    %         figure;
    %         hold on;
    %         plot(ft_fAPAR,'k');
    %         plot(w1_fAPAR);
    %         plot(w3_fAPAR);
    %         plot(w5_fAPAR);

            % save for each site

            save([out_dir,char(site_name),'_fAPAR.mat'],'ft_fAPAR');
            
        end
        
        
    end

end









% PART 5, calculate MODIS GPP using different fAPAR.
% output: site_specific GPP
% month, GPP using different fAPAR.

trig5=0;

if (trig5)
    
    meteo_dir='/Volumes/Elements/project14_AFF/AMF-footprint-analysis/met-data/';
    
    fAPAR_dir= '/Volumes/Elements/project14_AFF/MODIS_cutouts/weight_fAPAR/';
    
    out_dir = '/Volumes/Elements/project14_AFF/MODIS_cutouts/estGPP/';
    
    % modis GPP parameters
    % source : http://www.ntsg.umt.edu/files/modis/MOD17UsersGuide2015_v3.pdf
        parameter_table=nan(6,11);

        %pft_number 
        parameter_table(1,:)=[1,2,3,4,5,6,7,8,9,10,12];
        pft_names={'ENF';'EBF';'DNF';'DBF';...
        'MF';'CSH';'OSH';'WSA';'SAV';'GRA';...
        'CRO'};
    
        %max LUE kgC/m2/d/MJ
        parameter_table(2,:)=[0.000962,0.001268,0.001086...
            ,0.001165,0.001051,0.001281,0.000841...
            ,0.001239,0.001206,0.000860,0.001044];
        %Tmin_min C
        parameter_table(3,:)=[-8,-8,-8,-6,-7,-8,-8,-8,-8,-8,-8];
        %Tmin_max
        parameter_table(4,:)=[8.31,9.09,10.44,9.94,9.50,8.61,8.80,11.39,11.39,12.02,12.02];
        %VPD_min Pa, need to convert to kpa
        parameter_table(5,:)=[650,800,650,650,650,650,650,650,650,650,650]./1000;
        %VPD_max Pa to kPa
        parameter_table(6,:)=[4600,3100,2300,1650,2400,4700,4800,3200,3100,5300,4300]./1000;
    
    
          % get site name and PFT first.
          filename2=strcat('/Volumes/Elements/project14_AFF/MODIS_cutouts/','All_BASE-site_list.csv');
    
          M=readtable(filename2);

          name_list=table2array(M(:,1:2));   
        
        
        
     file_list=dir([fAPAR_dir,'*mat']);
     
     for i = 1:length(file_list) % for each site
         
         site_name = file_list(i).name(1:6);
         
         % load meteo data
         meteo_list=dir([meteo_dir,'*cat1_day_data.csv']);
         M=readtable([meteo_dir,meteo_list(1).name]);
         
         met_data=table2array(M);
         
         % load fAPAR data
         load([fAPAR_dir,file_list(i).name]); %ft_fAPAR
         
         
        % load PFT
        p_ind=strcmp(name_list(:,1),site_name);
        pft_site=name_list(p_ind,2);
        
        p2_ind=strcmp(pft_names,pft_site);
        
        if sum(p2_ind)>0 % if exist that PFT
            
            parameters=parameter_table(:,p2_ind);
            
        else
            
            parameters=parameter_table(:,end); % otherwise regard it as CRO
        end
         
        
        %calculate GPP for each fAPAR
        
        for ff = 1:4
            
            cd ('/Volumes/Elements/project14_AFF/code/');
            Tscale = ftair(met_data(:,3),parameters(3),parameters(4));
            Vscale = fvpd(met_data(:,5)./10,parameters(5),parameters(6)); % meteo VPD from hpa to kpa
            PAR = met_data(:,6).*0.5.* 3600 *24./1000000; % from SW W/m2 to PAR MJ/m2/d
            
            % fAPAR, SW, LUE, T and VPD
            % default unit is Kg/d
            % change to g/month
            estGPP(:,ff)=ft_fAPAR(:,ff).* PAR .* parameters(2) .* Tscale .*Vscale .*30. *1000;
           
            cumGPP = cumsum(estGPP,1);
            
            
        end
        
        
        save([out_dir,char(site_name),'_GPP.mat'],'estGPP','cumGPP');
         
     end
    
end











% PART 6, overview figure for each site
% include: 1, size of June footprint compared to pixel size
%          2. fAPAR dynamic of different area of interests
%          3. accumulated GPP for different fAPAR


trig6=1;


if (trig6)
    
    
    estGPP_dir='/Volumes/Elements/project14_AFF/MODIS_cutouts/estGPP/';
    
    ll_dir = '/Volumes/Elements/project14_AFF/MODIS_cutouts/sites_ll/'; % directory for surrounding pixels
    foot_dir = '/Volumes/Elements/project14_AFF/AMF-footprint-analysis/footprint-climatology/'; % footprint directory
    
    fAPAR_dir= '/Volumes/Elements/project14_AFF/MODIS_cutouts/weight_fAPAR/';
    
    file_list=dir([estGPP_dir,'*.mat']);
    
    for i=1:length(file_list)
    
        site_name=file_list(i).name(1:6);
        
        f1=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 25, 10], ...
        'OuterPosition', [2, 2, 25, 10]);
        
    
        % footprint of June
        cd('/Volumes/Elements/project6/code4Remi/functions');
        subaxis(1,3,1,'SpacingVert',0.05,'SpacingHoriz',0.1,'MR',0.05,'MarginTop',.05,'MarginBottom',.15); 
        
        load([ll_dir,char(site_name),'.mat']); % pixelll
        
        
        ft_files=dir([foot_dir,'*',char(site_name),'*.csv']);
        
            
            M=readtable(strcat(foot_dir,ft_files(6).name)); %only use June data
            tmp=table2array(M);
            
            % get the polygon for footprint, use ft size and model
            ind=tmp(:,3)==0.8 & tmp(:,4)==1;
            
            if sum(ind)==0 % if ft model 1 does not work, use the 2nd model
                ind=tmp(:,3)==0.8 & tmp(:,4)==2;
            end
            
            ft_poly=tmp(ind,1:2);
            
            for m=1:length(pixelll) % for every pixel
                
                poly1=polyshape(ft_poly(:,1),ft_poly(:,2)); %long and lat

                poly2=polyshape(pixelll(m,[8,6,7,9])',pixelll(m,[4,2,3,5])'); % long and lat

                % plot pixels
                pg=plot(poly2);
                pg.FaceColor = [0.8,0.8,0.8];
                pg.EdgeAlpha = 0.2;
                hold on;
            end
            
            % plot foot print
            pg2=plot(poly1);
            pg2.FaceColor = [1,0,0];

        xlabel('longitude','FontSize',12);
        ylabel('latitude','FontSize',12);
        
        
        
        % monthly fAPAR
        cd('/Volumes/Elements/project6/code4Remi/functions');
        subaxis(1,3,2,'SpacingVert',0.05,'SpacingHoriz',0.1,'MR',0.05,'MarginTop',.05,'MarginBottom',.15); 
        
        
        load([fAPAR_dir,char(site_name),'_fAPAR.mat']); % ft_fAPAR
        
        pp1=plot(1:12,ft_fAPAR);

        ylabel('fAPAR','FontSize',12);
        xlabel('month','FontSize',12); 
        
        ylim([0 1]);
        
        legend(pp1,{'footprint','1*1','3*3','5*5'},'Location','northwest');
        
        % cumulated GPP
        cd('/Volumes/Elements/project6/code4Remi/functions');
        subaxis(1,3,3,'SpacingVert',0.05,'SpacingHoriz',0.1,'MR',0.05,'MarginTop',.05,'MarginBottom',.15); 
        
        load([estGPP_dir,char(site_name),'_GPP.mat']); % ft_fAPAR
        
        pp2=plot(1:12,cumGPP);
        ylabel('cumulative GPP (g/m^2)','FontSize',12);
        xlabel('month','FontSize',12);     
        
        legend(pp2,{'footprint','1*1','3*3','5*5'},'Location','northwest');
        
        
        % out put figure
        set(f1,'PaperPositionMode','auto');
        print(strcat('/Volumes/Elements/project14_AFF/Figures/',char(site_name),'_overview'),'-djpeg','-r600');
        
        close(f1);
    
    end
    
    
    
end

