function Intensity_vs_stripe_corr(DATA)
%Intensity_vs_stripe_corr Summary of this function goes here
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
hfig.Name = 'Intensity vs stripe correlation';

clear hscat
if length(LBL_FACTOR) == 1
    ax = gobjects(3,1);
    hscat = gobjects(3,size(LBL_CONDITION,1));
    
    for idx = 1:3
        ax(idx) = subplot(1,3,idx);
        hold(ax(idx),'on');
        set(ax(idx),'Box','on')
    end    
    
    rho1  = zeros(size(LBL_CONDITION,1),1);
    rho2  = zeros(size(LBL_CONDITION,1),1);
    rho3  = zeros(size(LBL_CONDITION,1),1);
    pval1 = zeros(size(LBL_CONDITION,1),1);
    pavl2 = zeros(size(LBL_CONDITION,1),1);
    pval3 = zeros(size(LBL_CONDITION,1),1);
    for idxCond = 1:size(LBL_CONDITION,1)
        
        idxKeep1 = ismember(LineProf.Factor1,LBL_CONDITION{idxCond});
%         if length(LBL_FACTOR) == 1
%             idxKeep2 = true(size(idxKeep1));
%         elseif length(LBL_FACTOR) == 2
%             idxKeep2 = ismember(LineProf.Factor2,LBL_CONDITION{idxCond,2});
%         end
        
        HalfSpindle_abs  = LineProf{idxKeep1,'HalfSpindle_abs'};
        StripeLength_abs = LineProf{idxKeep1,'StripeLength_abs'};
        StripeHURP       = LineProf{idxKeep1,'StripeHURP'};
            
        [rho1(idxCond),pval1(idxCond)] = corr(StripeLength_abs(:),HalfSpindle_abs(:),'Type','Spearman');
        [rho2(idxCond),pval2(idxCond)] = corr(StripeLength_abs(:),StripeHURP(:),'Type','Spearman');
        [rho3(idxCond),pval3(idxCond)] = corr(HalfSpindle_abs(:),StripeHURP(:),'Type','Spearman');        
        
        hscat(1,idxCond) =  scatter(StripeLength_abs(:),HalfSpindle_abs(:),...
            'Parent',ax(1),...
            'MarkerFaceColor',CMAP_CONDITION(idxCond,:),...
            'MarkerEdgeColor',[0.2 0.2 0.2]);
        hscat(2,idxCond) = scatter(StripeLength_abs(:),StripeHURP(:),...
            'Parent',ax(2),...
            'MarkerFaceColor',CMAP_CONDITION(idxCond,:),...
            'MarkerEdgeColor',[0.2 0.2 0.2]);
        hscat(3,idxCond) = scatter(HalfSpindle_abs(:),StripeHURP(:),...
            'Parent',ax(3),...
            'MarkerFaceColor',CMAP_CONDITION(idxCond,:),...
            'MarkerEdgeColor',[0.2 0.2 0.2]);
        hscat(1,idxCond).MarkerFaceAlpha = .5;
        hscat(1,idxCond).MarkerEdgeAlpha = .9;
        hscat(2,idxCond).MarkerFaceAlpha = .5;
        hscat(2,idxCond).MarkerEdgeAlpha = .9;
        hscat(3,idxCond).MarkerFaceAlpha = .5;
        hscat(3,idxCond).MarkerEdgeAlpha = .9;
            
        
    end
    txt = [LBL_CONDITION{1} ': \rho_{Spearman}=' num2str(rho1(1),'%3.3f'),...
        ' (p=' strrep([num2str(pval1(1),'%3.3e') ')'],'e','\cdot10^{') '})'];
    if size(LBL_CONDITION,1) == 2
        txt = [txt  ' - ',...
        LBL_CONDITION{2} ': \rho_{Spearman}=' num2str(rho1(2),'%3.3f'),...
        ' (p=' strrep([num2str(pval1(2),'%3.3e') ')'],'e','\cdot10^{') '})'];
    end    
    title('Parent',ax(1),txt);
   
     txt = [LBL_CONDITION{1} ': \rho_{Spearman}=' num2str(rho2(1),'%3.3f'),...
        ' (p=' strrep([num2str(pval2(1),'%3.3e') ')'],'e','\cdot10^{') '})'];
    if size(LBL_CONDITION,1) == 2
        txt = [txt  ' - ',...
        LBL_CONDITION{2} ': \rho_{Spearman}=' num2str(rho2(2),'%3.3f'),...
        ' (p=' strrep([num2str(pval2(2),'%3.3e') ')'],'e','\cdot10^{') '})'];
    end    
    title('Parent',ax(2),txt);
    
    
     txt = [LBL_CONDITION{1} ': \rho_{Spearman}=' num2str(rho3(1),'%3.3f'),...
        ' (p=' strrep([num2str(pval3(1),'%3.3e') ')'],'e','\cdot10^{') '})'];
    if size(LBL_CONDITION,1) == 2
        txt = [txt  ' - ',...
        LBL_CONDITION{2} ': \rho_{Spearman}=' num2str(rho3(2),'%3.3f'),...
        ' (p=' strrep([num2str(pval3(2),'%3.3e') ')'],'e','\cdot10^{') '})'];
    end    
    title('Parent',ax(3),txt);
    
    for idx = 1:3
        legend(ax(idx),hscat(idx,:),LBL_CONDITION)
    end
    
    
    xlabel('Parent',ax(1),'Stripe Length (\mum)')
    ylabel('Parent',ax(1),'Half-spindle Length (\mum)')
    
    xlabel('Parent',ax(2),'Stripe Length (\mum)')
    ylabel('Parent',ax(2),'HURP Intensity (-)')
    
    xlabel('Parent',ax(3),'Half-spindle Length (\mum)')
    ylabel('Parent',ax(3),'HURP Intensity (-)')
    hold off
        
        
elseif length(LBL_FACTOR) == 2
    
    ax = gobjects(12,1);
    hscat = gobjects(12,2);
    for idx =1:12
        ax(idx) = subplot(4,3,idx);
        hold(ax(idx),'on');
        set(ax(idx),'Box','on')
    end
    
    
    all_lbl_condition = unique(LBL_CONDITION);
    rho1  = zeros(12,2);
    rho2  = zeros(12,2);
    rho3  = zeros(12,2);
    pval1 = zeros(12,2);
    pval2 = zeros(12,2);
    pval3 = zeros(12,2);
    
    for idxCond = 1:length(all_lbl_condition)
        [~,idxFactor] = find(ismember(LBL_CONDITION,all_lbl_condition{idxCond}));
        idxFactor = unique(idxFactor);
        idx_FixedCond = ismember(LineProf.(['Factor' num2str(idxFactor)]),all_lbl_condition{idxCond});
        CONDITION_OTHERFACTOR = unique(LBL_CONDITION(:,mod(idxFactor,2)+1));
        for idxOtherCond = 1:length(CONDITION_OTHERFACTOR)
           
            idxTestedCond = ismember(LineProf.(['Factor' num2str(mod(idxFactor,2)+1)]),CONDITION_OTHERFACTOR{idxOtherCond});
            
            HalfSpindle_abs  = LineProf{idx_FixedCond&idxTestedCond,'HalfSpindle_abs'};
            StripeLength_abs = LineProf{idx_FixedCond&idxTestedCond,'StripeLength_abs'};
            StripeHURP       = LineProf{idx_FixedCond&idxTestedCond,'StripeHURP'};
             
            
            idxRef1 = (idxCond-1)*3+1;
            idxRef2 = idxRef1+1;
            idxRef3 = idxRef2+1;
            
            [rho1(idxRef1,idxOtherCond),pval1(idxRef1,idxOtherCond)] = corr(StripeLength_abs(:),HalfSpindle_abs(:),'Type','Spearman');
            [rho2(idxRef2,idxOtherCond),pval2(idxRef2,idxOtherCond)] = corr(StripeLength_abs(:),StripeHURP(:),'Type','Spearman');
            [rho3(idxRef3,idxOtherCond),pval3(idxRef3,idxOtherCond)] = corr(HalfSpindle_abs(:),StripeHURP(:),'Type','Spearman');
            
            
            a = ismember(LBL_CONDITION,all_lbl_condition{idxCond});
            b = ismember(LBL_CONDITION,CONDITION_OTHERFACTOR{idxOtherCond});
            idxtmp = all((a|b),2);
            
            
            hscat(idxRef1,idxOtherCond) =  scatter(StripeLength_abs(:),HalfSpindle_abs(:),...
                'Parent',ax(idxRef1),...
                'MarkerFaceColor',CMAP_CONDITION(idxtmp,:),...
                'MarkerEdgeColor',[0.2 0.2 0.2]);
            hscat(idxRef2,idxOtherCond) = scatter(StripeLength_abs(:),StripeHURP(:),...
                'Parent',ax(idxRef2),...
                'MarkerFaceColor',CMAP_CONDITION(idxtmp,:),...
                'MarkerEdgeColor',[0.2 0.2 0.2]);
            hscat(idxRef3,idxOtherCond) = scatter(HalfSpindle_abs(:),StripeHURP(:),...
                'Parent',ax(idxRef3),...
                'MarkerFaceColor',CMAP_CONDITION(idxtmp,:),...
                'MarkerEdgeColor',[0.2 0.2 0.2]);
            
            hscat(idxRef1,idxOtherCond).MarkerFaceAlpha = .5;
            hscat(idxRef1,idxOtherCond).MarkerEdgeAlpha = .9;
            hscat(idxRef2,idxOtherCond).MarkerFaceAlpha = .5;
            hscat(idxRef2,idxOtherCond).MarkerEdgeAlpha = .9;
            hscat(idxRef3,idxOtherCond).MarkerFaceAlpha = .5;
            hscat(idxRef3,idxOtherCond).MarkerEdgeAlpha = .9;
           
             
        end

        title('Parent',ax(idxRef1),...
            {[LBL_FACTOR{idxFactor} ': ' all_lbl_condition{idxCond}];...
            [CONDITION_OTHERFACTOR{1} ': \rho_{Spearman}=' num2str(rho1(idxRef1,1),'%3.3f'),...
            ' (p=' strrep([num2str(pval1(idxRef1,1),'%3.3e') ')'],'e','\cdot10^{') '}) - ',...
            CONDITION_OTHERFACTOR{2} ': \rho_{Spearman}=' num2str(rho1(idxRef1,2),'%3.3f'),...
            ' (p=' strrep([num2str(pval1(idxRef1,2),'%3.3e') ')'],'e','\cdot10^{') '})']});
        title('Parent',ax(idxRef2),...
            {[LBL_FACTOR{idxFactor} ': ' all_lbl_condition{idxCond} ];...
            [CONDITION_OTHERFACTOR{1} ': \rho_{Spearman}=' num2str(rho2(idxRef2,1),'%3.3f'),...
            ' (p=' strrep([num2str(pval2(idxRef2,1),'%3.3e') ')'],'e','\cdot10^{') '}) - ',...
            CONDITION_OTHERFACTOR{2} ': \rho_{Spearman}=' num2str(rho2(idxRef2,2),'%3.3f'),...
            ' (p=' strrep([num2str(pval2(idxRef2,2),'%3.3e') ')'],'e','\cdot10^{') '})']});
        title('Parent',ax(idxRef3),...
            {[LBL_FACTOR{idxFactor} ': ' all_lbl_condition{idxCond} ];...
            [CONDITION_OTHERFACTOR{1} ': \rho_{Spearman}=' num2str(rho3(idxRef3,1),'%3.3f'),...
            ' (p=' strrep([num2str(pval3(idxRef3,1),'%3.3e') ')'],'e','\cdot10^{') '}) - ',...
            CONDITION_OTHERFACTOR{2} ': \rho_{Spearman}=' num2str(rho3(idxRef3,2),'%3.3f'),...
            ' (p=' strrep([num2str(pval3(idxRef3,2),'%3.3e') ')'],'e','\cdot10^{') '})']});
        
        legend(ax(idxRef1),[hscat(idxRef1,1) hscat(idxRef1,2)],CONDITION_OTHERFACTOR)
        
        legend(ax(idxRef2),[hscat(idxRef2,1) hscat(idxRef2,2)],CONDITION_OTHERFACTOR)
        
        legend(ax(idxRef3),[hscat(idxRef3,1) hscat(idxRef3,2)],CONDITION_OTHERFACTOR)
        
        xlabel('Parent',ax(idxRef1),'Stripe Length (\mum)')
        ylabel('Parent',ax(idxRef1),'Half-spindle Length (\mum)')
        
        xlabel('Parent',ax(idxRef2),'Stripe Length (\mum)')
        ylabel('Parent',ax(idxRef2),'HURP Intensity (-)')
        
        xlabel('Parent',ax(idxRef3),'Half-spindle Length (\mum)')
        ylabel('Parent',ax(idxRef3),'HURP Intensity (-)')
       
        
        
    end
    linkaxes(ax([1 2 4 5 7 8 10 11]),'x')
    linkaxes(ax([2 3 5 6 8 9 11 12]),'y')
    linkaxes(ax([3 6 9 12]),'x')
    linkaxes(ax([1 4 7 10]),'y')
    hold off
end

saveas(hfig,[filepath filesep hfig.Name '.fig'])

end

