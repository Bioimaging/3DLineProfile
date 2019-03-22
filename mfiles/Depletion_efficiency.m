function Depletion_efficiency(DATA)
%DEPLETION_EFFICIENCY Summary of this function goes here
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






if length(LBL_CONDITION) == 1
    return
end

choice = choosedialog(unique(LBL_CONDITION),'Depletion reference');

hfig = figure;
hfig.Name = 'Depletion efficiency';


if length(LBL_FACTOR) == 2
    
    ax(1) = subplot(1,2,1);
    ax(2) = subplot(1,2,2);
    
    tmpidx = any(ismember(LBL_CONDITION,choice));
    idxvarFactor = find(tmpidx);
    idxrefFactor = find(~tmpidx);
    refCondition = unique(LBL_CONDITION(:,idxrefFactor));
    varCondition = unique(LBL_CONDITION(:,idxvarFactor));
    %choice condition as the reference, ie the first in the list...
    tmpidx = single(ismember(varCondition,choice));
    tmpidx(tmpidx==0) = 2;
    varCondition = varCondition(tmpidx);
    
    for idxRefC = 1:length(refCondition)
        idxR = ismember(LineProf.(['Factor' num2str(idxrefFactor)]),refCondition{idxRefC});
        h_1 = cell(length(varCondition),1);
        h_2 = cell(length(varCondition),1);
        h_1_m = zeros(length(varCondition),1);
        h_2_m = zeros(length(varCondition),1);
        
        for idxVarC = 1:length(varCondition)
            idxV = ismember(LineProf.(['Factor' num2str(idxvarFactor)]),varCondition{idxVarC});
            a = ismember(LBL_CONDITION,varCondition{idxVarC});
            b = ismember(LBL_CONDITION,refCondition{idxRefC});
            idxcmap = all((a|b),2);
            
            HURP_spindle = LineProf{idxR&idxV,'HURP_spindle'};
            
            h_1{idxVarC} = HURP_spindle(:,1);
            h_2{idxVarC} = HURP_spindle(:,2);
            
            h_1_m(idxVarC) = median(HURP_spindle(:,1));
            h_2_m(idxVarC) = median(HURP_spindle(:,2));
            
            line(idxVarC+0.2*(rand(size(HURP_spindle(:,1)))-0.5),HURP_spindle(:,1),...
                'LineStyle','none','Marker','o',...
                'MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',CMAP_CONDITION(idxcmap,:),'Parent',ax(idxRefC))
            line(idxVarC+[-0.5*0.2 0.5*0.2],median(HURP_spindle(:,1))*ones(2,1),...
                'Color',[0 0 0],'LineWidth',2,'Parent',ax(idxRefC))
            line(idxVarC*ones(2,1),[prctile(HURP_spindle(:,1),25) prctile(HURP_spindle(:,1),75)],...
                'Color',[0 0 0],'LineWidth',2,'Parent',ax(idxRefC))
            
            line(idxVarC+2+0.2*(rand(size(HURP_spindle(:,2)))-0.5),HURP_spindle(:,2),...
                'LineStyle','none','Marker','o',...
                'MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',CMAP_CONDITION(idxcmap,:),'Parent',ax(idxRefC))
            line(idxVarC+2+[-0.5*0.2 0.5*0.2],median(HURP_spindle(:,2))*ones(2,1),...
                'Color',[0 0 0],'LineWidth',2,'Parent',ax(idxRefC))
            line(idxVarC+2*ones(2,1),[prctile(HURP_spindle(:,2),25) prctile(HURP_spindle(:,2),75)],...
                'Color',[0 0 0],'LineWidth',2,'Parent',ax(idxRefC))
            
            txt_xlbl{idxVarC}   = ['1st ' varCondition{idxVarC}];
            txt_xlbl{idxVarC+2} = ['2nd ' varCondition{idxVarC}];
            
        end
        [p1, test_txt1] = MakeunPairedTest1Factor(h_1{1},h_1{2},0.05);
        [p2, test_txt2] = MakeunPairedTest1Factor(h_2{1},h_2{2},0.05);
        ax(idxRefC).Title.String = {[refCondition{idxRefC} ', \epsilon_{1^{st}} = ' num2str(100*(h_1_m(1)-h_1_m(2))/h_1_m(1) ,'%02.2f') '%',...
            ', \epsilon_{2^{nd}} = ' num2str(100*(h_2_m(1)-h_2_m(2))/h_2_m(1) ,'%02.2f') '%'];...
            ['1^{st} ' test_txt1 ' , p=' num2str(p1) '; 2^{nd} ' test_txt2 ' , p=' num2str(p2)]};
        ylabel('HURP half-spindle','Parent',ax(idxRefC))
    end
    
    
    
    set(ax,'XTick',[1:length(txt_xlbl)],'XTickLabel',txt_xlbl)
    
    
    ylabel('HURP half-spindle','Parent',ax(2))
else
    ax = axes;
    %choice condition as the reference, ie the first in the list...
    varCondition = unique(LBL_CONDITION);
    tmpidx = single(ismember(varCondition,choice));
    tmpidx(tmpidx==0) = 2;
    varCondition = varCondition(tmpidx);
    h_1 = cell(length(varCondition),1);
    h_2 = cell(length(varCondition),1);
    h_1_m = zeros(length(varCondition),1);
    h_2_m = zeros(length(varCondition),1);
    for idxVarC = 1:length(varCondition)
        idxV = ismember(LineProf.Factor1,varCondition{idxVarC});
        
        idxcmap = idxVarC;
        
        HURP_spindle = LineProf{idxV,'HURP_spindle'};
        
        h_1{idxVarC} = HURP_spindle(:,1);
        h_2{idxVarC} = HURP_spindle(:,2);
        
        h_1_m(idxVarC) = median(HURP_spindle(:,1));
        h_2_m(idxVarC) = median(HURP_spindle(:,2));
        
        line(idxVarC+0.2*(rand(size(HURP_spindle(:,1)))-0.5),HURP_spindle(:,1),...
            'LineStyle','none','Marker','o',...
            'MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',CMAP_CONDITION(idxcmap,:),'Parent',ax)
        line(idxVarC+[-0.5*0.2 0.5*0.2],median(HURP_spindle(:,1))*ones(2,1),...
            'Color',[0 0 0],'LineWidth',2,'Parent',ax)
        line(idxVarC*ones(2,1),[prctile(HURP_spindle(:,1),25) prctile(HURP_spindle(:,1),75)],...
            'Color',[0 0 0],'LineWidth',2,'Parent',ax)
        
        line(idxVarC+2+0.2*(rand(size(HURP_spindle(:,2)))-0.5),HURP_spindle(:,2),...
            'LineStyle','none','Marker','o',...
            'MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',CMAP_CONDITION(idxcmap,:),'Parent',ax)
        line(idxVarC+2+[-0.5*0.2 0.5*0.2],median(HURP_spindle(:,2))*ones(2,1),...
            'Color',[0 0 0],'LineWidth',2,'Parent',ax)
        line(idxVarC+2*ones(2,1),[prctile(HURP_spindle(:,2),25) prctile(HURP_spindle(:,2),75)],...
            'Color',[0 0 0],'LineWidth',2,'Parent',ax)
        
        txt_xlbl{idxVarC}   = ['1st ' varCondition{idxVarC}];
        txt_xlbl{idxVarC+2} = ['2nd ' varCondition{idxVarC}];
        
    end
    [p1, test_txt1] = MakeunPairedTest1Factor(h_1{1},h_1{2},0.05);
    [p2, test_txt2] = MakeunPairedTest1Factor(h_2{1},h_2{2},0.05);
    ax.Title.String = {['\epsilon_{1^{st}} = ' num2str(100*(h_1_m(1)-h_1_m(2))/h_1_m(1) ,'%02.2f') '%',...
        ', \epsilon_{2^{nd}} = ' num2str(100*(h_2_m(1)-h_2_m(2))/h_2_m(1) ,'%02.2f') '%'];...
        ['1^{st} ' test_txt1 ' , p=' num2str(p1) '; 2^{nd} ' test_txt2 ' , p=' num2str(p2)]};
    ylabel('HURP half-spindle','Parent',ax)
    
end


linkaxes

saveas(hfig,[filepath filesep hfig.Name '.fig'])

end

