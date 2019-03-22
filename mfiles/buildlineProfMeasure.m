function  DATA = buildlineProfMeasure(LBL_DAPI,LBL_HURP,LBL_TUBULIN,LastFolderPath)
%BUILDLINEPROFMEASURE Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
%
%   v1.0 07-Feb-2019 NL

 prompt   = {'Enter DAPI wavelength:',...
            'Enter HURP wavelength:',...
            'Enter TUBULIN wavelength:',...
            'Enter factor number (1 or 2):'};
        title    = 'Input';
        dims     = [1 35];
        definput = {LBL_DAPI{:},LBL_HURP{:},LBL_TUBULIN{:},'2'};
        answer   = inputdlg(prompt,title,dims,definput);
        
        if isempty(answer)
            return
        end
        
        LBL_DAPI    = answer(1);%435
        LBL_HURP    = answer(2);%594
        LBL_TUBULIN = answer(3);%676
        NB_FACTOR   = str2double(answer{4});
        if NB_FACTOR~=1&&NB_FACTOR~=2
            opts = struct('WindowStyle','modal',...
                'Interpreter','none');
            f = errordlg('This is either 1 or 2 factor(s)!!!',...
                'Parameters Error', opts);
            return
        end
        
        LBL_FACTOR     = cell(1,NB_FACTOR);
        LBL_CONDITION  = cell(1,NB_FACTOR);
        for idxfactor = 1:NB_FACTOR
            prompt   = {['Enter the name of factor ' num2str(idxfactor) ':'],...
                'Enter the two condition separated by semicolon (;)'};
            title    = 'Input';
            dims     = [1 35];
            definput = {'Symmetry','DMSO;Taxol';'Depletion','siLA;siHURP'};
            answer   = inputdlg(prompt,title,dims,definput(idxfactor,:));
            LBL_FACTOR(idxfactor) = answer(1);
            LBL_CONDITION{idxfactor} = regexp(answer{2},';','Split');
        end
        
        tmp = unique((cellfun(@(x) length(x),LBL_CONDITION)));
        if length(tmp)~=1||tmp>2
            opts = struct('WindowStyle','modal',...
                'Interpreter','none');
            f = errordlg('You have to enter 2 conditions per factor !!!',...
                'Parameters Error', opts);
            return
        end
        for idxF = 1:length(LBL_FACTOR)
            if ~isvarname(LBL_FACTOR{idxF})
                opts = struct('WindowStyle','modal',...
                    'Interpreter','none');
                f = errordlg({['<< '  LBL_FACTOR{idxF} ' >> is not a valid name'];...
                    'A valid variable name starts with a letter, followed by letters, digits, or underscores, no space!'},...
                    'Naming Error', opts);
                return
            end
        end
        
        
        LBL_CONDITION = cat(1,LBL_CONDITION{:})';
        if NB_FACTOR == 2
            new_lbl_condition = cell(4,2);
            new_lbl_condition(1,:) = [LBL_CONDITION(1,1) LBL_CONDITION(1,2)];
            new_lbl_condition(2,:) = [LBL_CONDITION(1,1) LBL_CONDITION(2,2)];
            new_lbl_condition(3,:) = [LBL_CONDITION(2,1) LBL_CONDITION(1,2)];
            new_lbl_condition(4,:) = [LBL_CONDITION(2,1) LBL_CONDITION(2,2)];
            LBL_CONDITION = new_lbl_condition;
        end
        
        
        CMAP_CONDITION = lines(size(LBL_CONDITION,1));
        if NB_FACTOR == 2
            for idxcondition = 1:size(LBL_CONDITION,1)
                tittxt = [LBL_CONDITION{idxcondition,1} ' - ' LBL_CONDITION{idxcondition,2}];
                CMAP_CONDITION(idxcondition,:) = ChooseColor(tittxt,255*CMAP_CONDITION(idxcondition,:));
            end
        else
            for idxcondition = 1:size(LBL_CONDITION,1)
                tittxt = [LBL_CONDITION{idxcondition,1}];
                CMAP_CONDITION(idxcondition,:) = ChooseColor(tittxt,255*CMAP_CONDITION(idxcondition,:));
            end
        end
        
        
        
        PATH_CONDITION = cell(size(LBL_CONDITION,1),1);
        for idxcondition = 1:size(LBL_CONDITION,1)
            if NB_FACTOR == 1
                matfolder = uigetdir(LastFolderPath, ['Please Select -> ',...
                    LBL_FACTOR{1},...
                    ': ',...
                    LBL_CONDITION{idxcondition},...
                    ' folder']);
                if isnumeric(matfolder)
                    opts = struct('WindowStyle','modal',...
                        'Interpreter','none');
                    f = errordlg('No folder was selected',...
                        'Parameters Error', opts);
                    return
                end
                PATH_CONDITION{idxcondition} = matfolder;
                LastFolderPath = fileparts(matfolder);
            elseif NB_FACTOR == 2
                matfolder = uigetdir(LastFolderPath, ['Please Select -> ',...
                    LBL_FACTOR{1},...
                    ': ',...
                    LBL_CONDITION{idxcondition,1},...
                    ' - ',...
                    LBL_FACTOR{2},...
                    ': ',...
                    LBL_CONDITION{idxcondition,2},...
                    ' folder']);
                if isnumeric(matfolder)
                    opts = struct('WindowStyle','modal',...
                        'Interpreter','none');
                    f = errordlg('No folder was selected',...
                        'Parameters Error', opts);
                    return
                end
                PATH_CONDITION{idxcondition} = matfolder;
                LastFolderPath = fileparts(matfolder);
            end
        end
        
        LineProf = cell(size(LBL_CONDITION,1),1);
        for idxcondition = 1:size(LBL_CONDITION,1)
            LineProf{idxcondition} = LineLoader(PATH_CONDITION{idxcondition},...
                LBL_CONDITION(idxcondition,:),...
                {LBL_DAPI,LBL_HURP,LBL_TUBULIN});
        end
        LineProf = cat(1,LineProf{:});
        
        
        
       
        %**************************************************************************
        % Get the HURP maxima and positions
        %**************************************************************************
        LineProf.HURP_Max     = zeros(height(LineProf),2);
        LineProf.HURP_Max_idx = zeros(height(LineProf),2);
        for idxF = 1:height(LineProf)
            drel = LineProf{idxF,'drel'}{:};
            idxCut = floor(length(drel)/2);
            HURP_I = LineProf{idxF,'HURP_I'}{:};
            % get the 2 max
            [HURP_Max1, HURP_Max1_idx] = max(HURP_I([1:idxCut]));
            [HURP_Max2, HURP_Max2_idx] = max(HURP_I([idxCut+1:end]));
            HURP_Max2_idx = HURP_Max2_idx+idxCut;
            
            LineProf{idxF,'HURP_Max'}     = [HURP_Max1 HURP_Max2];
            LineProf{idxF,'HURP_Max_idx'} = [HURP_Max1_idx HURP_Max2_idx];
        end
        
        %**************************************************************************
        % Get the Tubulin at HURP average maxima positions
        % for data conditions AND at its own average maxima positions
        %**************************************************************************
        choosen_condition = choosedialog(unique(LBL_CONDITION(:)),...
            'Select the condition to extract the Tubulin at HURP median maxima positions:');
        
                
        if NB_FACTOR == 2
            [~,idxFactor] = find(ismember(LBL_CONDITION,choosen_condition));
            idxFactor = unique(idxFactor);
            idx_ChoosenCond = ismember(LineProf.(['Factor' num2str(idxFactor)]),choosen_condition);
            CONDITION_OTHERFACTOR = unique(LBL_CONDITION(:,mod(idxFactor,2)+1));
            HURP_mean_Max_idx = zeros(2,2);
            for idxCond = 1:length(CONDITION_OTHERFACTOR)
                idx_otherFactorCond = ismember(LineProf.(['Factor' num2str(mod(idxFactor,2)+1)]),...
                    CONDITION_OTHERFACTOR{idxCond});
                I = LineProf{idx_otherFactorCond&idx_ChoosenCond,'HURP_I'};
                I_HURP  = cat(1,I{:});
                med_I_HURP = median(I_HURP);
                idxCut = floor(length(med_I_HURP)/2);
                [~,HURP_mean_Max1_idx] = max(med_I_HURP([1:idxCut]));
                [~,HURP_mean_Max2_idx] = max(med_I_HURP([idxCut+1:end]));
                HURP_mean_Max2_idx = HURP_mean_Max2_idx+idxCut;
                HURP_mean_Max_idx(idxCond,:) =[HURP_mean_Max1_idx HURP_mean_Max2_idx];
            end
            LineProf.TUBL_at_HURP_mean_Max  = zeros(height(LineProf),2);
            LineProf.TUBL_mean_Max_idx = zeros(height(LineProf),2);
            LineProf.TUBL_mean_Max = zeros(height(LineProf),2);
            
            for idxCond = 1:size(LBL_CONDITION,1)
                idxCF1 = ismember(LineProf.Factor1,LBL_CONDITION(idxCond,1));              
                idxCF2 = ismember(LineProf.Factor2,LBL_CONDITION(idxCond,2));
                idx = find(idxCF1&idxCF2);
                TUBL_I = LineProf{idx,'TUBL_I'};
                I_TUBL = cat(1,TUBL_I{:});
                
                med_I_TUBL = median(I_TUBL);
                idxCut = floor(length(med_I_TUBL)/2);
                [~,TUBL_mean_Max1_idx] = max(med_I_TUBL([1:idxCut]));
                [~,TUBL_mean_Max2_idx] = max(med_I_TUBL([idxCut+1:end]));
                TUBL_mean_Max2_idx = TUBL_mean_Max2_idx+idxCut;
                
                LineProf.TUBL_mean_Max_idx(idx,:) = repmat([TUBL_mean_Max1_idx TUBL_mean_Max2_idx],[length(idx) 1]);
                
                
                idxOtherCond = ismember(CONDITION_OTHERFACTOR,LBL_CONDITION(1,:));
                
                for idxF = 1:length(TUBL_I)
                    
                    
                    LineProf{idx(idxF),'TUBL_at_HURP_mean_Max'} = TUBL_I{idxF}(HURP_mean_Max_idx(idxOtherCond,:));
                    LineProf{idx(idxF),'TUBL_mean_Max'} = TUBL_I{idxF}(LineProf{idx(idxF),'TUBL_mean_Max_idx'});
                end
            end
            
            
            
            
%             CONDF2 = unique(LBL_CONDITION(:,idxFactor));
%             for idxCondF1 = 1:length(CONDITION_OTHERFACTOR)
%                 idxCond1 = ismember(LineProf.Factor1,CONDITION_OTHERFACTOR{idxCondF1});
%                 for idxCondF2 = 1:length(CONDF2)
%                     idxCond2 = ismember(LineProf.Factor2,CONDF2{idxCondF2});
%                     idx = find(idxCond1&idxCond2);
%                     TUBL_I = LineProf{idx,'TUBL_I'};
%                     
%                     I_TUBL = cat(1,TUBL_I{:});
%                     med_I_TUBL = median(I_TUBL);
%                     idxCut = floor(length(med_I_TUBL)/2);
%                     [~,TUBL_mean_Max1_idx] = max(med_I_TUBL([1:idxCut]));
%                     [~,TUBL_mean_Max2_idx] = max(med_I_TUBL([idxCut+1:end]));
%                     TUBL_mean_Max2_idx = TUBL_mean_Max2_idx+idxCut;
%                     
%                     LineProf.TUBL_mean_Max_idx(idx,:) = repmat([TUBL_mean_Max1_idx TUBL_mean_Max2_idx],[length(idx) 1]);
%                     
%                     for idxF = 1:length(TUBL_I)
%                         LineProf{idx(idxF),'TUBL_at_HURP_mean_Max'} = TUBL_I{idxF}(HURP_mean_Max_idx(idxCondF1,:));
%                         LineProf{idx(idxF),'TUBL_mean_Max'} = TUBL_I{idxF}(LineProf{idx(idxF),'TUBL_mean_Max_idx'});
%                     end
%                     
%                 end
%             end
            
        else
            idx_ChoosenCond = ismember(LineProf.Factor1,choosen_condition);
            HURP_mean_Max_idx = zeros(1,2);
            I = LineProf{idx_ChoosenCond,'HURP_I'};
            I_HURP  = cat(1,I{:});
            med_I_HURP = median(I_HURP);
            idxCut = floor(length(med_I_HURP)/2);
            [~,HURP_mean_Max1_idx] = max(med_I_HURP([1:idxCut]));
            [~,HURP_mean_Max2_idx] = max(med_I_HURP([idxCut+1:end]));
            HURP_mean_Max2_idx = HURP_mean_Max2_idx+idxCut;
            HURP_mean_Max_idx =[HURP_mean_Max1_idx HURP_mean_Max2_idx];
            
            LineProf.TUBL_at_HURP_mean_Max  = zeros(height(LineProf),2);
            LineProf.TUBL_mean_Max_idx = zeros(height(LineProf),2);
            LineProf.TUBL_mean_Max = zeros(height(LineProf),2);
            for idxCondF = 1:length(LBL_CONDITION)
%                 idxCond2 = ismember(LineProf.Factor2,CONDF2{idxCondF});
                idxCond1 = ismember(LineProf.Factor1,LBL_CONDITION{idxCondF});
                idx = find(idxCond1);
                TUBL_I = LineProf{idx,'TUBL_I'};
                I_TUBL = cat(1,TUBL_I{:});
                med_I_TUBL = median(I_TUBL);
                idxCut = floor(length(med_I_TUBL)/2);
                [~,TUBL_mean_Max1_idx] = max(med_I_TUBL([1:idxCut]));
                [~,TUBL_mean_Max2_idx] = max(med_I_TUBL([idxCut+1:end]));
                TUBL_mean_Max2_idx = TUBL_mean_Max2_idx+idxCut;
                LineProf.TUBL_mean_Max_idx(idx,:) = repmat([TUBL_mean_Max1_idx TUBL_mean_Max2_idx],[length(idx) 1]);
                for idxF = 1:length(TUBL_I)
                    LineProf{idx(idxF),'TUBL_at_HURP_mean_Max'} = TUBL_I{idxF}(HURP_mean_Max_idx);
                    LineProf{idx(idxF),'TUBL_mean_Max'} = TUBL_I{idxF}(LineProf{idx(idxF),'TUBL_mean_Max_idx'});
                end
            end
            
        end
        
        
        
        
        %**************************************************************************
        % Find the metaphase plate with DAPI and get the 2 half-spindles
        % HURP intensities and Tubulin intensities
        %**************************************************************************
        LineProf.DAPI_metaphase_idx = zeros(height(LineProf),1);
        LineProf.HURP_spindle       = zeros(height(LineProf),2);
        LineProf.TUBL_spindle       = zeros(height(LineProf),2);
        clear HURP_spindle TUBL_spindle
        for idxF = 1:height(LineProf)
            
            DAPI_I = LineProf{idxF,'DAPI_I'}{:};
            [~, DAPI_metaphase_idx] = max(DAPI_I);
            
            HURP_I = LineProf{idxF,'HURP_I'}{:};
            HURP_spindle(1) = mean(HURP_I(1:DAPI_metaphase_idx-1));
            HURP_spindle(2) = mean(HURP_I(DAPI_metaphase_idx+1:end));
            
            LineProf{idxF,'DAPI_metaphase_idx'} = DAPI_metaphase_idx;
            LineProf{idxF,'HURP_spindle'} = HURP_spindle;
            
            TUBL_I = LineProf{idxF,'TUBL_I'}{:};
            TUBL_spindle(1) = mean(TUBL_I(1:DAPI_metaphase_idx-1));
            TUBL_spindle(2) = mean(TUBL_I(DAPI_metaphase_idx+1:end));
            LineProf{idxF,'TUBL_spindle'} = TUBL_spindle;
            
        end
        
        %**************************************************************************
        % Find the HURP stripes of each half-spindle,
        HURP_Iloss = 50;% in percentage from HurpMax1 HurpMax2 and the minimum in
        % between
        %**************************************************************************
        
        
        LineProf.StripeHURP        = zeros(height(LineProf),2);
        LineProf.HalfSpindle_abs   = zeros(height(LineProf),2);
        LineProf.StripeLength_abs  = zeros(height(LineProf),2);
        LineProf.HalfSpindle_rel   = zeros(height(LineProf),2);
        LineProf.StripeLength_rel  = zeros(height(LineProf),2);
        
        for idxF = 1:height(LineProf)
            HURP_Max = LineProf{idxF,'HURP_Max'};
            HURP_Max_idx = LineProf{idxF,'HURP_Max_idx'};
            HURP_I = LineProf{idxF,'HURP_I'}{:};
            dabs = LineProf{idxF,'dabs'}{:};
            drel = LineProf{idxF,'drel'}{:};
            dabs = dabs*drel/100;
            
            [HURP_Min,HURP_Min_idx] = min(HURP_I([HURP_Max_idx(1):HURP_Max_idx(2)]));
            HURP_Min_idx = HURP_Min_idx+HURP_Max_idx(1)-1;
            %get the stripes positions
            HURP_Stripe1_I = (HURP_Max(1)-HURP_Min)*HURP_Iloss/100 + HURP_Min;
            HURP_Stripe2_I = (HURP_Max(2)-HURP_Min)*HURP_Iloss/100 + HURP_Min;
            
            
            [~,HURP_Stripe1_idx(1)] = min(abs(HURP_I(1:HURP_Max_idx(1))-HURP_Stripe1_I));
            [~,HURP_Stripe1_idx(2)] = min(abs(HURP_I(HURP_Max_idx(1):HURP_Min_idx)-HURP_Stripe1_I));
            
            [~,HURP_Stripe2_idx(1)] = min(abs(HURP_I(HURP_Min_idx:HURP_Max_idx(2))-HURP_Stripe2_I));
            [~,HURP_Stripe2_idx(2)] = min(abs(HURP_I(HURP_Max_idx(2):end)-HURP_Stripe2_I));
            
            HURP_Stripe1_idx(2) = HURP_Stripe1_idx(2)+ HURP_Max_idx(1) -1;
            HURP_Stripe2_idx(1) = HURP_Stripe2_idx(1)+ HURP_Min_idx -1;
            HURP_Stripe2_idx(2) = HURP_Stripe2_idx(2)+ HURP_Max_idx(2) -1;
            
            LineProf{idxF,'StripeLength_rel'}  = [drel(HURP_Stripe1_idx(2)) - drel(HURP_Stripe1_idx(1)), drel(HURP_Stripe2_idx(2)) - drel(HURP_Stripe2_idx(1))];
            LineProf{idxF,'HalfSpindle_rel'}   = [drel(HURP_Min_idx)-drel(1) drel(end)-drel(HURP_Min_idx)];
            
            LineProf{idxF,'StripeLength_abs'} = [dabs(HURP_Stripe1_idx(2)) - dabs(HURP_Stripe1_idx(1)), dabs(HURP_Stripe2_idx(2)) - dabs(HURP_Stripe2_idx(1))];
            LineProf{idxF,'HalfSpindle_abs'}  = [dabs(HURP_Min_idx)-dabs(1) dabs(end)-dabs(HURP_Min_idx)];
            LineProf{idxF,'StripeHURP'}       = [trapz(HURP_I(HURP_Stripe1_idx(1):HURP_Stripe1_idx(2))) trapz(HURP_I(HURP_Stripe2_idx(1):HURP_Stripe2_idx(2)))];
            
            
            %    line(drel,HURP_I)
            %    line([drel(HURP_Stripe1_idx(1)) drel(HURP_Stripe1_idx(2))],[HURP_I(HURP_Stripe1_idx(1)) HURP_I(HURP_Stripe1_idx(2))]);
            %    line([drel(HURP_Stripe2_idx(1)) drel(HURP_Stripe2_idx(2))],[HURP_I(HURP_Stripe2_idx(1)) HURP_I(HURP_Stripe2_idx(2))]);
            
            
            
        end
        
        idxKill = any(LineProf.StripeLength_abs == 0,2);
        if any(idxKill)
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            disp('The following cells will be removed because of a 0 stripe length:')
            disp(LineProf.FileNames(idxKill))
            LineProf(idxKill,:) = [];
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        end
        
 
%         
%         save([pwd filesep 'mfiles' filesep 'DefaultOptions.mat'],...
%             'LastFolderPath','LBL_DAPI','LBL_HURP','LBL_TUBULIN','NB_FACTOR',...
%             '-append');
        
             
          varname = {'LineProf',...
            'LBL_FACTOR',...
            'LBL_CONDITION',...
            'PATH_CONDITION',...
            'LBL_DAPI',...
            'LBL_HURP',...
            'LBL_TUBULIN',...
            'CMAP_CONDITION'};
        
        filepath = Myuisave(varname,'LineProfile');
        filepath = fileparts(filepath);    
        
        DATA.CMAP_CONDITION = CMAP_CONDITION;
        DATA.LBL_CONDITION  = LBL_CONDITION;
        DATA.LBL_DAPI       = LBL_DAPI;
        DATA.LBL_FACTOR     = LBL_FACTOR;
        DATA.LBL_HURP       = LBL_HURP;
        DATA.LBL_TUBULIN    = LBL_TUBULIN;
        DATA.LineProf       = LineProf;
        DATA.PATH_CONDITION = PATH_CONDITION;
        DATA.filepath = filepath;
        
end

