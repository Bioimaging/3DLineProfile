function LineProfile3D_Comparison(DATA)
%LineProfile3D_Comparison Summary of this function goes here
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

all_lbl_condition = unique(LBL_CONDITION);
if length(all_lbl_condition) > 1
    hfig = figure;
    hfig.Name = '3D Line Profile Comparison';
end



if length(all_lbl_condition) == 2
    ax = subplot(1,1,1);
    hpatch = gobjects(2,1);
    hold(ax,'on');
    ax2 = axes('Position',get(ax,'Position'),...
            'XAxisLocation','top',...
            'YTick',{},...
            'YTickLabel',{},...
            'YAxisLocation','right',...
            'Color','none');   
        
        drel = LineProf.drel{1}+0.5;
        drel_lbl = cellfun(@(x) ['d_' num2str(x)] ,num2cell(drel),'UniformOutput',false)';
        drel = table(drel','VariableNames',{'RelativeDist'});
        
        RMA_table = array2table(cat(1,LineProf{:,'HURP_I'}{:}));
        RMA_table = [LineProf.Factor1 RMA_table];        
        RMA_table.Properties.VariableNames =[ LBL_FACTOR ; drel_lbl];
                
        
        rm   = fitrm(RMA_table,[drel_lbl{1} '-' drel_lbl{end} '~ ' LBL_FACTOR{:}],...
            'WithinDesign',drel,'WithinModel','orthogonalcontrasts');
        ranovatbl = ranova(rm);
        mauchlyT  = mauchly(rm);
        epsT      = epsilon(rm);
        D_total   = multcompare(rm,LBL_FACTOR{:},'By','RelativeDist');
        for idxCond = 1:length(LBL_CONDITION)
            idx = ismember(RMA_table{:,1},LBL_CONDITION{idxCond});
            %      m = mean(Rmatable{idx,[2:end]},1);
            %      s = std(Rmatable{idx,[2:end]},[],1);
            %      sem = s./length(s)^0.5;
            
            m = median(RMA_table{idx,[2:end]},1);
            %      s = std(Rmatable{idx,[2:end]},[],1);
            quart1 = prctile(RMA_table{idx,[2:end]},25,1);
            quart3 = prctile(RMA_table{idx,[2:end]},75,1);
            
            d = drel.RelativeDist-0.5;

            line(d,m,...
                'Color',CMAP_CONDITION(idxCond,:),'LineStyle','-',...
                'LineWidth',2,'Marker','o','MarkerSize',5,'Parent',ax,...
                'MarkerFaceColor',CMAP_CONDITION(idxCond,:),...
                'MarkerEdgeColor',[0.2 0.2 0.2]);
            
            hpatch(idxCond) = patch([d' fliplr(d')],[quart1 fliplr(quart3) ],...
                CMAP_CONDITION(idxCond,:),...
                'FaceAlpha',0.3,'EdgeColor',CMAP_CONDITION(idxCond,:),'Parent',ax);                        
        end
        v = get(ax,'YLim');        
        set(ax,'YLim',[v(1) v(2)*1.1]);
        pval = D_total.pValue([1:2:end]);
        txt = cell(length(d),1);
        for idxD =1:length(d)
            if pval(idxD)<0.001
                lbl = '***';
            elseif pval(idxD)>=0.001&&pval(idxD)<0.01
                lbl = '**';
            elseif pval(idxD)>=0.01&&pval(idxD)<0.05
                lbl = '*';
            elseif pval(idxD)>=0.05
                lbl = 'ns';
            end
            txt{idxD} = ['p=' strrep([num2str(pval(idxD),'%0.3e')],'e','\cdot10^{') '}' ];
                text(d(idxD),v(2),{[lbl ]},...
                'HorizontalAlignment','center','Parent',ax,'FontSize',10,...
                'Rotation',0)
           
        end                
        
        ylabel('Parent',ax,'HURP Intensity (-)')
        xlabel('Parent',ax,'Relative distance(%)')
        
        set(ax2,'XTick',d,'XTickLabel',txt,'XLim',[0 100],...
            'XTickLabelRotation',30,'FontSize',8)
        
        title(ax2,{[LBL_CONDITION{1} ' vs ' LBL_CONDITION{2} ', Q1 median Q3'];...
            ['Repeated measures ANOVA with a Greenhouse-Geisser correction: F(',...
            num2str(ranovatbl.DF(1)*epsT.GreenhouseGeisser,'%0.3f') ',',...
            num2str(ranovatbl.DF(3)*epsT.GreenhouseGeisser,'%0.3f') ')=' num2str(ranovatbl.F(1),'%0.3f'),...
            '; p=' strrep([num2str(ranovatbl.pValueGG(1),'%0.3e') ],'e','\cdot10^{') '}']} )
        ax2.Title.FontSize = 8;
        
%         title('Q1 median Q3','Parent',ax(idxCond+2*size(LBL_CONDITION,1)))
%         ylabel('Median profiles (-)','Parent',ax(idxCond+2*size(LBL_CONDITION,1)))
%         xlabel('Relative distance (%)','Parent',ax(idxCond+2*size(LBL_CONDITION,1)))
        legend(hpatch,LBL_CONDITION)
        hold off
        
        

elseif length(all_lbl_condition) == 4    
    ax = gobjects(2,2);
    for idx =1:4
        ax(idx) = subplot(2,2,idx);
        hold(ax(idx),'on');
%         set(ax(idx),'Box','on')
    end
    ax2 = gobjects(2,2);
    hpatch = gobjects(4,2);
    for idx  = 1:4
        ax2(idx) = axes('Position',get(ax(idx),'Position'),...
            'XAxisLocation','top',...
            'YTick',{},...
            'YTickLabel',{},...
            'YAxisLocation','right',...
            'Color','none');
    end
    
    
    for idxCond = 1:length(all_lbl_condition)
        [~,idxFactor] = find(ismember(LBL_CONDITION,all_lbl_condition{idxCond}));
        idxFactor = unique(idxFactor);
        idx_FixedCond = ismember(LineProf.(['Factor' num2str(idxFactor)]),all_lbl_condition{idxCond});
        CONDITION_OTHERFACTOR = unique(LBL_CONDITION(:,mod(idxFactor,2)+1));
        
        drel = LineProf.drel{1}+0.5;
        drel_lbl = cellfun(@(x) ['d_' num2str(x)] ,num2cell(drel),'UniformOutput',false)';
        drel = table(drel','VariableNames',{'RelativeDist'});
        
        RMA_table = array2table(cat(1,LineProf{idx_FixedCond,'HURP_I'}{:}));
        RMA_table = [LineProf(idx_FixedCond,['Factor' num2str(mod(idxFactor,2)+1)]) RMA_table];
%         RMA_table = LineProf(idx_FixedCond,{['Factor' num2str(mod(idxFactor,2)+1)], 'HURP_I'})
        RMA_table.Properties.VariableNames =[ LBL_FACTOR(mod(idxFactor,2)+1) ; drel_lbl];
        
        
        
        rm   = fitrm(RMA_table,[drel_lbl{1} '-' drel_lbl{end} '~ ' LBL_FACTOR{mod(idxFactor,2)+1}],...
            'WithinDesign',drel,'WithinModel','orthogonalcontrasts');
        ranovatbl = ranova(rm);
        mauchlyT  = mauchly(rm);
        epsT      = epsilon(rm);
        D_total = multcompare(rm,LBL_FACTOR{mod(idxFactor,2)+1},'By','RelativeDist');
        
       
        for idxOtherCond = 1:length(CONDITION_OTHERFACTOR)
            idx = ismember(RMA_table{:,1},CONDITION_OTHERFACTOR{idxOtherCond});
            %      m = mean(Rmatable{idx,[2:end]},1);
            %      s = std(Rmatable{idx,[2:end]},[],1);
            %      sem = s./length(s)^0.5;
            
            m = median(RMA_table{idx,[2:end]},1);
            %      s = std(Rmatable{idx,[2:end]},[],1);
            quart1 = prctile(RMA_table{idx,[2:end]},25,1);
            quart3 = prctile(RMA_table{idx,[2:end]},75,1);
            
            d = drel.RelativeDist-0.5;
            
            a = ismember(LBL_CONDITION,all_lbl_condition{idxCond});
            b = ismember(LBL_CONDITION,CONDITION_OTHERFACTOR{idxOtherCond});
            idxtmp = all((a|b),2);
            
            line(d,m,...
                'Color',CMAP_CONDITION(idxtmp,:),'LineStyle','-',...
                'LineWidth',2,'Marker','o','MarkerSize',5,'Parent',ax(idxCond),...
                'MarkerFaceColor',CMAP_CONDITION(idxtmp,:),...
                'MarkerEdgeColor',[0.2 0.2 0.2]);
            
            hpatch(idxCond,idxOtherCond) = patch([d' fliplr(d')],[quart1 fliplr(quart3) ],...
                CMAP_CONDITION(idxtmp,:),...
                'FaceAlpha',0.3,'EdgeColor',CMAP_CONDITION(idxtmp,:),'Parent',ax(idxCond));
            
        end
       
        v = get(ax(idxCond),'YLim');        
        set(ax(idxCond),'YLim',[v(1) v(2)*1.1]);
        pval = D_total.pValue([1:2:end]);
        txt = cell(length(d),1);
        for idxD =1:length(d)
            if pval(idxD)<0.001
                lbl = '***';
            elseif pval(idxD)>=0.001&&pval(idxD)<0.01
                lbl = '**';
            elseif pval(idxD)>=0.01&&pval(idxD)<0.05
                lbl = '*';
            elseif pval(idxD)>=0.05
                lbl = 'ns';
            end
            txt{idxD} = ['p=' strrep([num2str(pval(idxD),'%0.3e')],'e','\cdot10^{') '}' ];
                text(d(idxD),v(2),{[lbl ]},...
                'HorizontalAlignment','center','Parent',ax(idxCond),'FontSize',10,...
                'Rotation',0)
           
        end                
        
        ylabel('Parent',ax(idxCond),'HURP Intensity (-)')
        xlabel('Parent',ax(idxCond),'Relative distance(%)')
        
        set(ax2(idxCond),'XTick',d,'XTickLabel',txt,'XLim',[0 100],...
            'XTickLabelRotation',30,'FontSize',8)
        
        title(ax2(idxCond),{[LBL_FACTOR{idxFactor} ': ' all_lbl_condition{idxCond},...
            ' - ' CONDITION_OTHERFACTOR{1} ' vs ' CONDITION_OTHERFACTOR{2} ', Q1 median Q3'];...
            ['Repeated measures ANOVA with a Greenhouse-Geisser correction: F(',...
            num2str(ranovatbl.DF(1)*epsT.GreenhouseGeisser,'%0.3f') ',',...
            num2str(ranovatbl.DF(3)*epsT.GreenhouseGeisser,'%0.3f') ')=' num2str(ranovatbl.F(1),'%0.3f'),...
            '; p=' strrep([num2str(ranovatbl.pValueGG(1),'%0.3e') ],'e','\cdot10^{') '}']} )
        ax2(idxCond).Title.FontSize = 8;
        
%         title('Q1 median Q3','Parent',ax(idxCond+2*size(LBL_CONDITION,1)))
%         ylabel('Median profiles (-)','Parent',ax(idxCond+2*size(LBL_CONDITION,1)))
%         xlabel('Relative distance (%)','Parent',ax(idxCond+2*size(LBL_CONDITION,1)))
        legend(hpatch(idxCond,:),CONDITION_OTHERFACTOR)

    end
    linkaxes(ax,'xy')
    hold off
end
if length(all_lbl_condition) > 1
    
    saveas(hfig,[filepath filesep hfig.Name '.fig'])
end

end

