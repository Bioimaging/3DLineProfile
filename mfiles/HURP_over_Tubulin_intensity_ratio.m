function HURP_over_Tubulin_intensity_ratio(DATA)
%HURP_over_Tubulin_intensity_ratio Summary of this function goes here
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
hfig.Name = 'HURP over Tubulin intensity ratio';

ax = subplot(1,1,1);

set(ax,'XTick',[1:size(LBL_CONDITION,1)],'Box','on','XLIM',[1-0.5 size(LBL_CONDITION,1)+0.5])

txt_xlbl = cell(size(LBL_CONDITION,1),1);
hr = cell(size(LBL_CONDITION,1),1);
for idxC = 1:size(LBL_CONDITION,1)
    idxKeep1 = ismember(LineProf.Factor1,LBL_CONDITION{idxC});
    if length(LBL_FACTOR) == 1
        idxKeep2 = true(size(idxKeep1));
        txt_xlbl{idxC} = LBL_CONDITION{idxC};
    elseif length(LBL_FACTOR) == 2
        idxKeep2 = ismember(LineProf.Factor2,LBL_CONDITION{idxC,2});
        txt_xlbl{idxC} = [LBL_CONDITION{idxC,1} ' - ' LBL_CONDITION{idxC,2}];
    end
    HURP_spindle = LineProf{idxKeep1&idxKeep2,'HURP_spindle'};
    HURP_spindle_ratio = HURP_spindle(:,2)./HURP_spindle(:,1);
    
    TUBL_spindle = LineProf{idxKeep1&idxKeep2,'TUBL_spindle'};
    TUBL_spindle_ratio = TUBL_spindle(:,2)./TUBL_spindle(:,1);
    
    HT_ratio = HURP_spindle_ratio./TUBL_spindle_ratio;
    hr{idxC} = HT_ratio;
    
    
    line(idxC+0.2*(rand(size(HT_ratio))-0.5),HT_ratio,...
        'LineStyle','none','Marker','o',...
        'MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',CMAP_CONDITION(idxC,:),...
        'Parent',ax)
    line(idxC+[-0.5*0.1 0.5*0.1],median(HT_ratio)*ones(2,1),...
        'Color',[0 0 0],'LineWidth',2,'Parent',ax)
    line(idxC*ones(2,1),[prctile(HT_ratio,25) prctile(HT_ratio,75)],...
        'Color',[0 0 0],'LineWidth',1,'Parent',ax)
    line(idxC+[-0.5*0.05 0.5*0.05],prctile(HT_ratio,25)*ones(2,1),...
        'Color',[0 0 0],'LineWidth',2,'Parent',ax)
    line(idxC+[-0.5*0.05 0.5*0.05],prctile(HT_ratio,75)*ones(2,1),...
        'Color',[0 0 0],'LineWidth',2,'Parent',ax)     
  
end
title(ax,'HURP over Tubulin half-spindle ratios, Q1 median Q3')
line(get(ax,'XLIM'),[1 1],'Color',[0.2 0.2 0.2],'LineStyle',':','Parent',ax)
set(ax,'XTickLabel',txt_xlbl)
ylabel('HURP over Tubulin half-spindle ratio (-)')


if size(LBL_CONDITION,1) == 2
    
    [ p,test_txt ] = MakeunPairedTest1Factor(hr{1},hr{2},0.05);
    title(ax,{'HURP over Tubulin half-spindle ratios, Q1 median Q3';...
        [test_txt ' , p=' strrep([num2str(p,'%0.3e') ],'e','\cdot10^{') '}']} )
    
    
saveas(hfig,[filepath filesep hfig.Name '.fig'])
    
elseif  size(LBL_CONDITION,1) == 4
    
     
    r = groot;
    r.ShowHiddenHandles = 'on';
   
saveas(hfig,[filepath filesep hfig.Name '.fig'])
    
    HTr = (LineProf.HURP_spindle(:,2)./LineProf.HURP_spindle(:,1))./(LineProf.TUBL_spindle(:,2)./LineProf.TUBL_spindle(:,1));
    [p,tbl,stats,terms]  = anovan(HTr,...
        {LineProf.Factor1,LineProf.Factor2},'model','interaction','varnames',LBL_FACTOR);        
    ht = findobj(gcf,'Type','Uicontrol','-and','Tag','Heading');
    ht.String = 'Analysis of HURP over Tubulin half-spindle ratio Variance';
    figure
    [c,m,h,gnames] = multcompare(stats,'Dimension',[1 2]);
    set(gcf,'Name',['HURP over Tubulin half-spindle ratio 2-way post-hoc ' get(gcf,'Name')])
    writetable(table([gnames(c(:,1)) gnames(c(:,2)) num2cell(c(:,end))]),[filepath filesep  'HURP over Tubulin ratio half-spindle.xlsx'])
    
    
end

end

