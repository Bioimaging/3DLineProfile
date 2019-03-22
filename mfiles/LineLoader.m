function LineProf = LineLoader(Folder,Condition,LBL_Wave)
%LINELOADER Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
%
%   v1.0 07-Feb-2019 NL

LBL_DAPI    = LBL_Wave{1};
LBL_HURP    = LBL_Wave{2};
LBL_TUBULIN = LBL_Wave{3};


if length(Condition)==2
    txtcondition = [Condition{1} ' ' Condition{2}];
else
    txtcondition = Condition{1};
end


FileNames = dir([Folder filesep '*.mat']);
FileNames = {FileNames.name}';

ratioManualD = zeros(length(FileNames),1);
drel         = cell(length(FileNames),1);
dabs         = cell(length(FileNames),1);
HURP_I       = cell(length(FileNames),1);
TUBL_I       = cell(length(FileNames),1);
DAPI_I       = cell(length(FileNames),1);


isSymFile = false;

if exist([Folder filesep 'symmetry.xlsx']) == 2   
    [xls_r ,xlsnames] = xlsread([Folder filesep 'symmetry.xlsx']);
    xls_r = xls_r(2:end,1);
    xlsnames = cellfun(@(x) [x(1:end-4) '.mat'], xlsnames(2:end,1),...
        'UniformOutput',false);
    isSymFile = true;
else
    opts = struct('WindowStyle','modal',...
        'Interpreter','none');
    h = warndlg(['No symmetry.xlsx file was found for the ' txtcondition,...
        ' condition!!!'],...
        'Missing file', opts);
    waitfor(h)
end


hwb = waitbar(0,'','Name',txtcondition);
set(findall(hwb), 'Units', 'Normalized')
set(hwb, 'Units', 'Pixels')
hwb.Position(3) = hwb.Position(3)+150;
hwb.Position(1) = hwb.Position(1)-50;
tmp = findobj(hwb,'Type','Axes');
tmp.Title.Interpreter = 'none';
tmp.Title.FontSize = 10;


last_time = clock;
init_time = last_time;
for idxF = 1:length(FileNames)    
    if etime(clock,last_time) >= 0.1
        progress = idxF/length(FileNames);
        elap = etime(clock,init_time);
        sec_remain = elap*(1/progress-1);
        r_mes = datestr(sec_remain/86400,'HH:MM:SS');
        waitbar(progress,hwb,...
            {FileNames{idxF};
            [num2str(idxF) '/' num2str(length(FileNames)),...
            ' - Remaining time: ' r_mes]})       
        last_time = clock;
    end
    
    if isSymFile        
        idx1 = ismember(xlsnames,[FileNames{idxF}]);
        idx2 = ismember(xlsnames,[FileNames{idxF}(1:end-4) 'D.mat']);%because of the "D" missing before the .mat...
        idx3 = ismember(xlsnames,[FileNames{idxF}(1:end-5) '.mat']);%because of the "D" that sometime has to be missing before the .mat...
        idx = idx1|idx2|idx3;
        if any(idx)
            ratioManualD(idxF) = xls_r(idx,1);
        end
        if isempty(ratioManualD(idxF))
            warndlg([FileNames{idxF} ' in ' txtcondition,...
                ' has no distance ratio!!!'])
        end
    end
  
    
    DATA = load([Folder filesep  FileNames{idxF}]);
    DATA = DATA.DATA;
    
    drel{idxF} = DATA.Intensity.relative.d{:};
    dabs{idxF} = (sum((DATA.SPOTB_XYZ_um-DATA.SPOTA_XYZ_um).^2))^0.5;%DATA.Intensity.absolute.d{:};
    
    %Get HURP SIGNAL
    tmp = find(cellfun(@(x) ~isempty(x),regexp(DATA.Intensity.chname,LBL_HURP)));
    mch = ['mCh' num2str(tmp)];
    hurp_i = DATA.Intensity.relative.(mch){:};
    mch = ['Ch' num2str(tmp)];
    HurpBackground = double(quantile(DATA.Intensity.cylinder.(mch),0.1));
    HURP_I{idxF} = hurp_i-HurpBackground;
    
    %Get TUBULIN SIGNAL
    tmp = find(cellfun(@(x) ~isempty(x),regexp(DATA.Intensity.chname,LBL_TUBULIN)));
    mch = ['mCh' num2str(tmp)];
    tubulin_i = DATA.Intensity.relative.(mch){:};
    mch = ['Ch' num2str(tmp)];
    TubulinBackground = double(quantile(DATA.Intensity.cylinder.(mch),0.1));
    TUBL_I{idxF} = tubulin_i-TubulinBackground;
    
    %Get DAPI SIGNAL
    tmp = find(cellfun(@(x) ~isempty(x),regexp(DATA.Intensity.chname,LBL_DAPI)));
    mch = ['mCh' num2str(tmp)];
    dapi_i = DATA.Intensity.relative.(mch){:};
    DAPI_I{idxF} = dapi_i;
    
end
delete(hwb)
 

if length(Condition)==2    
    Factor1 = repmat(Condition(1),[length(FileNames),1]);
    Factor2 = repmat(Condition(2),[length(FileNames),1]);
    LineProf = table(Factor1,Factor2,...
                    FileNames,...
                    drel,dabs,...
                    HURP_I,TUBL_I,DAPI_I,...
                    ratioManualD);
else
    Factor1 = repmat(Condition(1),[length(FileNames),1]);    
    LineProf = table(Factor1,...
                    FileNames,...
                    drel,dabs,...
                    HURP_I,TUBL_I,DAPI_I,...
                    ratioManualD);
end

