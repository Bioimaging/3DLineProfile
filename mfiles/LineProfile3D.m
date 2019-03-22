function LineProfile3D(DATA)
%LineProfile3D Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
%
%   v1.0 07-Feb-2019 NL

CMAP_CONDITION = DATA.CMAP_CONDITION;
LBL_CONDITION  = DATA.LBL_CONDITION;
LBL_FACTOR     = DATA.LBL_FACTOR;
LineProf       = DATA.LineProf;
filepath       = DATA.filepath;

hfig = figure;
hfig.Name = '3D Line Profile';


for idxCond = 1:size(LBL_CONDITION,1)
    idxKeep1 = ismember(LineProf.Factor1,LBL_CONDITION{idxCond});
    if length(LBL_FACTOR) == 1
        idxKeep2 = true(size(idxKeep1));
    elseif length(LBL_FACTOR) == 2
        idxKeep2 = ismember(LineProf.Factor2,LBL_CONDITION{idxCond,2});
    end
    drel   = LineProf{idxKeep1&idxKeep2,'drel'};
    HURP_I = LineProf{idxKeep1&idxKeep2,'HURP_I'};
    TUBL_I = LineProf{idxKeep1&idxKeep2,'TUBL_I'};
    
    d = drel{1};
    I_HURP = cat(1,HURP_I{:});
    med_I_HURP = median(I_HURP);
    p25_I_HURP = prctile(I_HURP,25);
    p75_I_HURP = prctile(I_HURP,75);
    
    I_TUBL = cat(1,TUBL_I{:});
    med_I_TUBL = median(I_TUBL);
    p25_I_TUBL = prctile(I_TUBL,25);
    p75_I_TUBL = prctile(I_TUBL,75);
    
    %_____________________________________________________
    ax(idxCond) = subplot(3,size(LBL_CONDITION,1),idxCond);
    for idxF = 1:length(HURP_I)
        line(drel{idxF},HURP_I{idxF},'Parent',ax(idxCond),'Color',CMAP_CONDITION(idxCond,:))
    end
    ylabel('HURP intensity profile (-)','Parent',ax(idxCond))
    xlabel('Relative distance (%)','Parent',ax(idxCond))
    %         title([LBL_FACTOR{:} ': ' LBL_CONDITION{idxCond}],'Parent',ax(idxCond))
    %_____________________________________________________
    ax(idxCond+size(LBL_CONDITION,1)) = subplot(3,size(LBL_CONDITION,1),idxCond+size(LBL_CONDITION,1));
    for idxF = 1:length(TUBL_I)
        line(drel{idxF},TUBL_I{idxF},'Parent',ax(idxCond+size(LBL_CONDITION,1)),'Color',CMAP_CONDITION(idxCond,:))
    end
    ylabel('Tubulin intensity profile (-)','Parent',ax(idxCond+size(LBL_CONDITION,1)))
    xlabel('Relative distance (%)','Parent',ax(idxCond+size(LBL_CONDITION,1)))
    %         title([LBL_FACTOR{:} ': ' LBL_CONDITION{idxCond}],'Parent',ax(idxCond+size(LBL_CONDITION,1)))
    %_____________________________________________________
    ax(idxCond+2*size(LBL_CONDITION,1)) = subplot(3,size(LBL_CONDITION,1),idxCond+2*size(LBL_CONDITION,1));
    hold on
    hp(1) = patch([d, fliplr(d)]',[p25_I_TUBL, fliplr(p75_I_TUBL)]',1,...
        'FaceColor',CMAP_CONDITION(idxCond,:),'EdgeColor',CMAP_CONDITION(idxCond,:),...
        'FaceAlpha',0.25,'Parent',ax(idxCond+2*size(LBL_CONDITION,1)));
    line(d,med_I_TUBL,'Color',CMAP_CONDITION(idxCond,:),'LineWidth',2,'Parent',ax(idxCond+2*size(LBL_CONDITION,1)),...
        'Marker','none')
    hp(2) = patch([d, fliplr(d)]',[p25_I_HURP, fliplr(p75_I_HURP)]',1,...
        'FaceColor',CMAP_CONDITION(idxCond,:),'EdgeColor',CMAP_CONDITION(idxCond,:),...
        'FaceAlpha',0.5,'Parent',ax(idxCond+2*size(LBL_CONDITION,1)));
    line(d,med_I_HURP,'Color',CMAP_CONDITION(idxCond,:),'LineWidth',2,'Parent',ax(idxCond+2*size(LBL_CONDITION,1)),'Marker','none')
    hold off    
    title('Q1 median Q3','Parent',ax(idxCond+2*size(LBL_CONDITION,1)))
    ylabel('Median profiles (-)','Parent',ax(idxCond+2*size(LBL_CONDITION,1)))
    xlabel('Relative distance (%)','Parent',ax(idxCond+2*size(LBL_CONDITION,1)))
    legend(hp,{'Tubulin','HURP'})
    if length(LBL_FACTOR) == 1
        title([LBL_FACTOR{:} ': ' LBL_CONDITION{idxCond}],'Parent',ax(idxCond))        
    elseif length(LBL_FACTOR) == 2
        title([LBL_FACTOR{1} ': ' LBL_CONDITION{idxCond,1} ' - ',...
            LBL_FACTOR{2} ': ' LBL_CONDITION{idxCond,2}],'Parent',ax(idxCond))
    end
    
end
set(ax,'box','on')
for idxAx=1:3
    linkaxes(ax(1:size(LBL_CONDITION,1)),'y')
    linkaxes(ax(size(LBL_CONDITION,1)+1:2*size(LBL_CONDITION,1)),'y')
    linkaxes(ax(2*size(LBL_CONDITION,1)+1:3*size(LBL_CONDITION,1)),'y')
end
linkaxes(ax,'x')
saveas(hfig,[filepath filesep hfig.Name '.fig'])



end


