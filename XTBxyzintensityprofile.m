% Extract intensity profile between two spots in 3D
%
%    <CustomTools>
%      <Menu name="Bioimaging XT">
%       <Submenu name="3D">
%        <Item name="I-profile between 2 spots" icon="Matlab" tooltip="Intensity profile between 2 spots in 3D">
%          <Command>MatlabXT::XTBxyzintensityprofile(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%    </CustomTools>
%
% Copyright (c) Feb 2017, Bioimaging Core Facility
% v1. Nicolas Liaudet




function [ output_args ] = XTBxyzintensityprofile(aImarisApplicationID)
%XTNSEG Summary of this function goes here
%   Extract intensity profile between two spots in 3D.
%    <CustomTools>
%      <Menu>
%        <Item name="I-profile between 2 spots" icon="Matlab" tooltip="Intensity profile between 2 spots in 3D">
%          <Command>MatlabXT::XTBxyzintensityprofile(%i)</Command>
%        </Item>
%      </Menu>
%    </CustomTools>




%% connect to Imaris interface
disp('Imaris connection...')
tic
if ~isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
    javaaddpath ImarisLib.jar
    vImarisLib = ImarisLib;
    if ischar(aImarisApplicationID)
        aImarisApplicationID = round(str2double(aImarisApplicationID));
    end
    vImarisApplication = vImarisLib.GetApplication(aImarisApplicationID);
else
    vImarisApplication = aImarisApplicationID;
end
vDataSet = vImarisApplication.GetDataSet;
 %vImarisDataSet = vImarisApplication.GetDataSet.Clone;
% vImarisDataSet = vImarisApplication.GetDataSet;
toc



%% Get DATA geometry
DataSizeXYZTC = [vDataSet.GetSizeX, vDataSet.GetSizeY,...
    vDataSet.GetSizeZ, vDataSet.GetSizeT,...
    vDataSet.GetSizeC];

DataPosXYZ = [vDataSet.GetExtendMinX, vDataSet.GetExtendMaxX;...
              vDataSet.GetExtendMinY, vDataSet.GetExtendMaxY;...
              vDataSet.GetExtendMinZ, vDataSet.GetExtendMaxZ];
          
DataResXYZ = diff(DataPosXYZ,1,2)./DataSizeXYZTC(1:3)';

%% get the spots

vSpots = vImarisApplication.GetFactory.ToSpots(vImarisApplication.GetSurpassSelection);
if isempty(vSpots)
    errordlg('Please select a "Spots" object!','Selection Error');        
    return
end

Chname = cell(1,DataSizeXYZTC(5));
for idxC = 1:DataSizeXYZTC(5) 
    Chname{idxC} = char(vDataSet.GetChannelName(idxC-1));
end


XYZ        = vSpots.GetPositionsXYZ;
Tidx       = vSpots.GetIndicesT;


button = questdlg('Define the origine according to:','Give Origin','Selection','Intensity','Intensity');

switch button
    case 'Intensity'
        %------------------- SPOT REF ACCORDGING TO MAX INTENSITY VAL ------------
        %get spots A and B acording to their intensity
        if rem(length(Tidx),2) ~= 0
             errordlg('Some spots are not always there!','Dataset Error');
            return
        end
        choice = chooseSpotRef(Chname);
        stats = vSpots.GetStatistics;
        statnames = cell(stats.mNames);
        idx = ismember(statnames,'Intensity Max');
        aFactors = transpose(cell(stats.mFactors));
        aFactors = aFactors(idx,2);
        values = stats.mValues;
%         ids = stats.mIds;
        values = values(idx);
%         ids = ids(idx);
        idx = ismember(aFactors,num2str(choice));
%         ids = ids(idx);
        values = values(idx);
        [~,idx] = sort(values,'Descend');
                        
        SpotAidx = idx(1);
        SpotBidx = idx(2);
        
        
        
    case 'Selection'
        %------------------- SPOT REF ACCORDGING TO USER SELECTION IN IMARIS---------
        %Reassign the A or B labelling according to the highligthed spot in the current scene
        
        vIds = vSpots.GetIds;
        SpotAidx =vIds(1);
        SpotBidx =vIds(2);
        SelectedSpots = vIds(vSpots.GetSelectedIndices+1);
               
        
        if ~isempty(SelectedSpots)
            if ~(SpotAidx == SelectedSpots)
                
                tmp1 = SpotAidx;
                tmp2 = SpotBidx;
                
                SpotBidx = tmp1;
                SpotAidx = tmp2;
            end
        end
        %find whom is whom here
        SpotAidx = find(vIds==SpotAidx);
        SpotBidx = find(vIds==SpotBidx);
        
end


TidxA = Tidx(SpotAidx);
% TidxB = Tidx(SpotBidx);

[~,pTA] = sort(TidxA);
[~,pTB] = sort(TidxA);

SPOTA_XYZ_um  = XYZ(SpotAidx,:);
SPOTB_XYZ_um  = XYZ(SpotBidx,:);

SPOTA_XYZ_um = SPOTA_XYZ_um(pTA,:);
SPOTB_XYZ_um = SPOTB_XYZ_um(pTB,:);

SPOTA_XYZ = SPOTA_XYZ_um-repmat(DataPosXYZ(:,1)',size(SPOTA_XYZ_um,1),1);
SPOTB_XYZ = SPOTB_XYZ_um-repmat(DataPosXYZ(:,1)',size(SPOTB_XYZ_um,1),1);
SPOTA_XYZ = SPOTA_XYZ./repmat(DataResXYZ',size(SPOTA_XYZ_um,1),1); 
SPOTB_XYZ = SPOTB_XYZ./repmat(DataResXYZ',size(SPOTB_XYZ_um,1),1);
SPOTA_XYZ = round(SPOTA_XYZ);
SPOTB_XYZ = round(SPOTB_XYZ);

SPOTB_XYZ_um  = XYZ(SpotBidx,:);

t = cell(1,length(TidxA));
for idxt = 1:length(TidxA)
    t(idxt) = vSpots.GetTimePoint(TidxA(idxt));
end




%%

[selection,ok] = listdlg('PromptString','Select channel to be analyzed:',...
                             'SelectionMode','multiple',...
                             'ListString',Chname);

if ~ok
    errordlg('No channel selected...','Selection Error');        
    return
end
                         

Intensity.chname = Chname(selection);
Intensity.chid   = selection;

Intensity.absolute.d      = cell(1,length(TidxA));%distances in absolute units
Intensity.relative.d      = cell(1,length(TidxA));%distances in relative units 100% is SpotA to SpotB
for k=1:length(selection)    
    mch = ['mCh' num2str(k)];
    sch = ['sCh' num2str(k)];
    Intensity.absolute.(mch) = cell(1,length(TidxA));%average intensity at a given distance in Ch
    Intensity.absolute.(sch) = cell(1,length(TidxA));%std intensity at a given distance in Ch    
    Intensity.relative.(mch) = cell(1,length(TidxA));%average intensity at a given rel. distance in Ch
    Intensity.relative.(sch) = cell(1,length(TidxA));%std intensity at a given rel. distance in Ch
end


    
prompt = {'Enter the radius (\mum):';...
          'Enter the distance step along the axis from A to B (\mum)';...
          'Enter the distance step along the axis (A to B is 100%)  (%)'};
dlg_title  = 'Integration cylinder';
num_lines  = 1;

ds     = 0.3;% in um
dperc = 0.03;% in%
defaultans = {'1',num2str(ds) ,num2str(100*dperc)};

options.Interpreter = 'Tex';
answer = inputdlg(prompt,dlg_title,num_lines,defaultans,options);

rcyl = str2double(answer{1});
if isempty(rcyl)||isnan(rcyl)||rcyl<=0
    errordlg('Please give a radius which is a number >0','Input Error');        
    return
end

ds = str2double(answer{2});
if isempty(ds)||isnan(ds)||ds<=0
    errordlg('Please give a distance step which is a number >0','Input Error');        
    return
end

dperc = str2double(answer{3})/100;
if isempty(dperc)||isnan(dperc)||dperc<=0
    errordlg('Please give a relative distance step radius which is a number >0 and <100','Input Error');        
    return
end

Intensity.rcyl  = rcyl;
Intensity.ds    = ds;
Intensity.dperc = dperc;



DataSizeXYZTC(5) = DataSizeXYZTC(5)+2;
      
vDataSet.SetSizeC(DataSizeXYZTC(5));
  
        
        
vDataSet.SetChannelName(DataSizeXYZTC(5)-2,['Link'])
vDataSet.SetChannelRange(DataSizeXYZTC(5)-2,0,1)

vDataSet.SetChannelName(DataSizeXYZTC(5)-1,['Geodesic distance'])
vDataSet.SetChannelRange(DataSizeXYZTC(5)-1,0,max((sum((SPOTB_XYZ-SPOTA_XYZ).^2,2)).^0.5))


hwbar = waitbar(0, ['Processing...']);


for idxT=1:DataSizeXYZTC(4)
    hwbar = waitbar(idxT/(DataSizeXYZTC(4)), hwbar,...
        ['Processing z-stack ',...
        num2str(idxT) '/' num2str(DataSizeXYZTC(4))]);
    
    if strcmp(vDataSet.GetType,'eTypeUInt8')
        tmp = zeros([DataSizeXYZTC(1:3)],'uint8');
    elseif strcmp(vDataSet.GetType,'eTypeUInt16')
        tmp = zeros([DataSizeXYZTC(1:3)],'uint16');
    elseif strcmp(vDataSet.GetType,'eTypeFloat')
        tmp = zeros([DataSizeXYZTC(1:3)],'single');
    end
    
    a_xyz = SPOTA_XYZ(idxT,:);
    b_xyz = SPOTB_XYZ(idxT,:);
   
    %make the cylinder axis    
    [link_x,link_y,link_z] = bresenham_line3d(a_xyz,b_xyz, 0);               
    link_ind = sub2ind(DataSizeXYZTC([1 2 3]),link_x,link_y,link_z);
    tmp(link_ind) = 1;        
    vDataSet.SetDataVolumeAs1DArrayShorts(tmp(:),DataSizeXYZTC(5)-2 , idxT-1 )
    
    
    %keep voxels closer than xxx um away from the axis
    D = bwdist(logical(tmp));    
    D = uint16(D<= rcyl/DataResXYZ(1)) ;            
    idx_D = find(D~=0);    
    [Dx,Dy,Dz] = ind2sub(DataSizeXYZTC([2 1 3]),idx_D);
    Dxyz = [Dx,Dy,Dz];
    
    
    %kill voxel behind A    
    n_AB = b_xyz-a_xyz;
    D_a_xyz = Dxyz-repmat(a_xyz,[size(Dxyz,1) 1 ]);    
    Pscal_abDa = dot(repmat(n_AB,[size(D_a_xyz,1) 1 ]),D_a_xyz,2);
    idx_rmv_A = Pscal_abDa<0;    
    D(idx_D(idx_rmv_A)) = 0;


    %kill voxel behind B
    n_BA = -n_AB;
    D_b_xyz = Dxyz-repmat(b_xyz,[size(Dxyz,1) 1 ]);    
    Pscal_baDb = dot(repmat(n_BA,[size(D_b_xyz,1) 1 ]),D_b_xyz,2);
    idx_rmv_B = Pscal_baDb<0;
    D(idx_D(idx_rmv_B)) = 0;
    
    
   
    
    %find "A face"
    Bc = strel(ones(3,3,3));
    gD = D-imerode(D,Bc);
    
    tmp = false(size(D));
    dm = floor(rcyl/DataResXYZ(1));
    
    tmp(a_xyz(1),a_xyz(2),a_xyz(3)) = true;
    
    sw = (2*dm-1)/2;
    ses2 = ceil(2*dm/2);
    [y,x,z] = meshgrid(-sw:sw, -sw:sw, -sw:sw);
    m = sqrt(x.^2 + y.^2 + z.^2);
    bc = (m <= m(ses2, ses2, 2*dm));
    Bc = strel('arbitrary', bc);    
    
    tmp = imdilate(tmp,Bc);
        
    msk = gD~=0;
    msk = msk&tmp;
    
    
    %distance along the cylinder
    D = bwdistgeodesic(D~=0,msk,'quasi-euclidean');
    
    
    vDataSet.SetDataVolumeAs1DArrayShorts(uint16(D(:)),DataSizeXYZTC(5)-1 , idxT-1 )
    
    %group voxels according to the distance with ds um steps
    D = D*DataResXYZ(1);%in um...   
    abs_edges = [0:ds:floor(max(D(:)))];
    
    %group voxels according to percentage of the max ditance
    rel_edges = linspace(0,max(D(:)),round(1/dperc+1));
    
    Intensity.absolute.d = {abs_edges(1:end-1)+ds/2};
    Intensity.relative.d = {100*[dperc/2:dperc:1-dperc/2]};
    
    
    
    [~,D_Bin_abs] = histcountsmex(D(:),abs_edges);
    [~,D_Bin_rel] = histcountsmex(D(:),rel_edges);
    
    
    for idxC = 1:length(Intensity.chid)        
        
        mch = ['mCh' num2str(idxC)];
        sch = ['sCh' num2str(idxC)];
                        
        if strcmp(vDataSet.GetType,'eTypeUInt8')
            Ch = zeros([DataSizeXYZTC(1:3)],'uint8');
            Ch(:) = typecast(vDataSet.GetDataVolumeAs1DArrayShorts(...
                Intensity.chid(idxC)-1,idxT-1), 'uint8');            
        elseif strcmp(vDataSet.GetType,'eTypeUInt16')
            Ch = zeros([DataSizeXYZTC(1:3)],'uint16');
            Ch(:) = typecast(vDataSet.GetDataVolumeAs1DArrayShorts(...
                Intensity.chid(idxC)-1,idxT-1), 'uint16');      
        elseif strcmp(vDataSet.GetType,'eTypeFloat')
            Ch = zeros([DataSizeXYZTC(1:3)],'single');
            Ch(:) = typecast(vDataSet.GetDataVolumeAs1DArrayShorts(...
                Intensity.chid(idxC)-1,idxT-1), 'single');      
        end
        
        m = zeros(1,length(abs_edges)-1);
        s = zeros(1,length(abs_edges)-1);       
        for didx = 1:length(abs_edges)-1                        
            Binidx = D_Bin_abs==didx; 
            valch = Ch(Binidx);                        
            m(didx) = mean(valch);
            s(didx) = std(single(valch));            
        end
        Intensity.absolute.(mch) = {m};
        Intensity.absolute.(sch) = {s};
                
        m = zeros(1,length(rel_edges)-1);
        s = zeros(1,length(rel_edges)-1);        
        for didx = 1:length(rel_edges)-1            
            Binidx = D_Bin_rel==didx; 
            valch = Ch(Binidx);                        
            m(didx) = mean(valch);
            s(didx) = std(single(valch));             
        end
        
        m = smooth(m,3)';
        
        Intensity.relative.(mch) = {m};
        Intensity.relative.(sch) = {s};
        
    end
       
end




close(hwbar)





FilePathName = char(vImarisApplication.GetCurrentFileName);



makefigure(SPOTA_XYZ_um,SPOTB_XYZ_um,Intensity,FilePathName)







DATA.FilePathName = FilePathName;
DATA.Intensity    = Intensity;
% 
DATA.SPOTA_XYZ_um = SPOTA_XYZ_um;
DATA.SPOTB_XYZ_um = SPOTB_XYZ_um;
 
save([FilePathName(1:end-4) '.mat'],'DATA')



    function makefigure(SPOTA_XYZ_um,SPOTB_XYZ_um,Intensity,FilePathName)
        
        [pathstr,filename,fileext] = fileparts(FilePathName);
        
        
        % non normalized distance
        figure('PaperPosition',[0.635 0.635 28.431 19.731],'Color',[1 1 1],...
            'PaperSize',[29.7 21],'units','normalized','outerposition',[0 0 1 1],'Name',[filename,fileext])
        
      
        ax = zeros(1,length(Intensity.chid)*2+2);
        
        for idx = 1:length(ax)
            
            if idx == 1 
                
                ax(idx) = subplot(2,length(Intensity.chid)+1,[1 length(Intensity.chid)+2]);
            elseif idx == length(Intensity.chid)+2 
            
            else
            
                ax(idx) = subplot(2,length(Intensity.chid)+1,idx);
            end
        end
        
      
        %--------------------------------------------------------------------------
        hp(1) = line(SPOTA_XYZ_um(:,1),SPOTA_XYZ_um(:,2),SPOTA_XYZ_um(:,3),'Parent',ax(1),...
            'Color',[1 0 0],'Marker','o','MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[1 0 0],...
            'LineStyle','none');
        hp(2) = line(SPOTB_XYZ_um(:,1),SPOTB_XYZ_um(:,2),SPOTB_XYZ_um(:,3),'Parent',ax(1),...
            'Color',[0 1 0],'Marker','o','MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[0 1 0],...
            'LineStyle','none');
        
        
        xlabel(ax(1),'x(\mum)')
        ylabel(ax(1),'y(\mum)')
        zlabel(ax(1),'z(\mum)')
        
        title(ax(1),['Spots configuration, d_{AB}=' num2str((sum((SPOTA_XYZ_um-SPOTB_XYZ_um).^2,2)).^0.5) '\mum'])
        legend(ax(1),hp,{'Spot A','Spot B'},'Location','southoutside')
        set(ax(1),'View',[-37.5 30],'Box','on','BoxStyle','full',...
            'XGrid','on','YGrid','on','ZGrid','on','DataAspectRatio',[1 1 1])
        
      
        %--------------------------------------------------------------------------
        for idxC = 1:length(Intensity.chid)
            mch = ['mCh' num2str(idxC)];
            sch = ['sCh' num2str(idxC)];
            d = Intensity.absolute.d{:};
            m = Intensity.absolute.(mch){:};
            s = Intensity.absolute.(sch){:};
            
            h(1) = line(d,m,'Parent',ax(idxC+1),'MarkerSize',4,...
                'Marker','o','MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[lines(1)]);
            h(2) = patch([d fliplr(d)],[m+s fliplr(m-s)],1,...
                'FaceColor',lines(1)+[0.2 0.2 0.2],'EdgeColor','none',...
                'FaceAlpha',0.3,'Parent',ax(idxC+1));
            
            xlabel(ax(idxC+1),'d(\mum)')
            ylabel(ax(idxC+1),'I(-)')
            title(ax(idxC+1),['Channel ' Intensity.chname{idxC}])
            legend(h,'Mean','\pm std')
            set(ax(idxC+1),'XGrid','on','YGrid','on','Box','on')
            
            mch = ['mCh' num2str(idxC)];
            sch = ['sCh' num2str(idxC)];
            d = Intensity.relative.d{:};
            m = Intensity.relative.(mch){:};
            s = Intensity.relative.(sch){:};
            
            h(1) = line(d,m,'Parent',ax(idxC+length(Intensity.chid)+2),'MarkerSize',4,...
                'Marker','o','MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[lines(1)]);
            h(2) = patch([d fliplr(d)],[m+s fliplr(m-s)],1,...
                'FaceColor',lines(1)+[0.2 0.2 0.2],'EdgeColor','none',...
                'FaceAlpha',0.3,'Parent',ax(idxC+length(Intensity.chid)+2));
            
            xlabel(ax(idxC+length(Intensity.chid)+2),'d(%)')
            ylabel(ax(idxC+length(Intensity.chid)+2),'I(-)')
            title(ax(idxC+length(Intensity.chid)+2),['Channel ' Intensity.chname{idxC}])
            legend(h,'Mean','\pm std')
            set(ax(idxC+length(Intensity.chid)+2),'XGrid','on','YGrid','on','Box','on')
            
        end
           
 
        
        
        print('-dpng','-r300',[FilePathName(1:end-4) '.png'],'-loose')
        
    end





end


function choice = chooseSpotRef(Chnames)

scrsz = get(0,'ScreenSize');


d = dialog('Position',[(scrsz(3)-250)/2 (scrsz(4)-150)/2 250 150],'Name','Spot Reference');
txt = uicontrol('Parent',d,...
    'Style','text',...
    'Position',[20 80 210 40],...
    'String','Select a channel');

popup = uicontrol('Parent',d,...
    'Style','popup',...
    'Position',[75 70 100 25],...
    'String',Chnames,...
    'Callback',@popup_callback);

btn = uicontrol('Parent',d,...
    'Position',[89 20 70 25],...
    'String','Ok',...
    'Callback','delete(gcf)');

choice = 1;
movegui(d,'center')

uiwait(d);


    function popup_callback(popup,callbackdata)
        idx = popup.Value;
        popup_items = popup.String;
        choice = get(popup,'Value');
        %           popup_items = get(popup,'String');
        %           choice = char(popup_items(idx,:));
    end
end