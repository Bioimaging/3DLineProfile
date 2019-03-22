function HURP_Iratio_vs_HalfSpindleratio_corr(DATA)
%HURP_Iratio_vs_HalfSpindleratio_corr Summary of this function goes here
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

if isempty(LineProf.ratioManualD) || (length(LineProf.ratioManualD)~= height(LineProf))
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    disp('Some manual ratio distance are missing: cannot proceed to the correlation analysis!!!!')
    return
end

hfig = figure;
hfig.Name = 'HURP intensity ratio vs half-spindle ratio correlation';

ax = gobjects(size(LBL_CONDITION,1));

    

for idxCond = 1:size(LBL_CONDITION,1)
    idxCF1 = ismember(LineProf.Factor1,LBL_CONDITION(idxCond,1));
    if length(LBL_FACTOR) == 2
        idxCF2 = ismember(LineProf.Factor2,LBL_CONDITION(idxCond,2));
    else
        idxCF2 = true(size(idxCF1));
    end
    HURP_spindle = LineProf{idxCF1&idxCF2,'HURP_spindle'};
    HURP_spindle_ratio = HURP_spindle(:,2)./HURP_spindle(:,1);
    ratioManualD = LineProf{idxCF1&idxCF2,'ratioManualD'};
    
    ax(idxCond) = subplot(1,size(LBL_CONDITION,1),idxCond);
    line(ratioManualD,HURP_spindle_ratio,...
        'LineStyle','none','Marker','o',...
        'MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',CMAP_CONDITION(idxCond,:),'Parent',ax(idxCond))
    
    [r,p] = corr(ratioManualD,HURP_spindle_ratio,'type','Spearman');
    
    if length(LBL_FACTOR) == 1
        title({[LBL_CONDITION{idxCond}]  ;...
            ['\rho_{Spearman} = ' num2str(r) '(p = ' num2str(p) ')']},'Parent',ax(idxCond))
    else
        title({[LBL_CONDITION{idxCond,1} ' - ' LBL_CONDITION{idxCond,2}]  ;...
            ['\rho_{Spearman} =  ' num2str(r) '(p = ' num2str(p) ')']},'Parent',ax(idxCond))        
    end
    
    ylabel('HURP half-spindle ratio','Parent',ax(idxCond))
    xlabel('Distance ratio','Parent',ax(idxCond))
    
    set(ax(idxCond),'box','on',...
        'XLim',[min([ratioManualD;0.5]) max([ratioManualD;3])],...
        'YLim',[min([HURP_spindle_ratio;0.5]) max([HURP_spindle_ratio;3])])
end

linkaxes
for idxA =1:length(ax)
    line(get(ax(idxA),'XLim'),ones(2,1),'Color',[0 0 0],'Parent',ax(idxA))
    line(ones(2,1),get(ax(idxA),'YLim'),'Color',[0 0 0],'Parent',ax(idxA))
end

saveas(hfig,[filepath filesep hfig.Name '.fig'])


end

