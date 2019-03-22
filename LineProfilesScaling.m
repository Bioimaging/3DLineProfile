% 3D line profile scaling
%
%   Nicolas Liaudet
%
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
%
%   v3.0 07-Feb-2019 NL
%   v3.1 14-Mar-2019 NL
%   Change names, removed features

clc
clear
close all

cd(fileparts(which('LineProfilesScaling.m')))
addpath(genpath('mfiles'))

DefaultOptions = load([pwd filesep 'mfiles' filesep 'DefaultOptions.mat']);
LastFolderPath = DefaultOptions.LastFolderPath;
LBL_DAPI       = DefaultOptions.LBL_DAPI;
LBL_HURP       = DefaultOptions.LBL_HURP;
LBL_TUBULIN    = DefaultOptions.LBL_TUBULIN;
NB_FACTOR      = DefaultOptions.NB_FACTOR;


answer = questdlg('Do you want to create a new empty symetry file ?','','Yes','No','Yes');
switch answer
    case ''        
    case 'No'        
    case 'Yes'
        makeSymetryXls
end


answer = questdlg('Do you want to load a LineProfile.mat ?','','Yes','No','Yes');
switch answer
    case ''
        return
    case 'Yes'
        [file,filepath] = uigetfile({'*LineProfile*.mat',...
            'Line Profile analysis file (*LineProfile*.mat)'},...
            'Select a file',LastFolderPath);
        if isnumeric(file)
            return
        end
        LastFolderPath = fileparts(filepath);
        DATA = load([LastFolderPath filesep file]);   
        DATA.filepath = filepath;
        save([filepath file],...
        'filepath','-append');

    case 'No'
        DATA = buildlineProfMeasure(LBL_DAPI,LBL_HURP,LBL_TUBULIN,LastFolderPath);       
        LastFolderPath = DATA.filepath;
end


LBL_DAPI       = DATA.LBL_DAPI;
LBL_HURP       = DATA.LBL_HURP;
LBL_TUBULIN    = DATA.LBL_TUBULIN;
LBL_FACTOR     = DATA.LBL_FACTOR;
LBL_CONDITION  = DATA.LBL_CONDITION;
LineProf       = DATA.LineProf;
NB_FACTOR      = length(LBL_FACTOR);

save([pwd filesep 'mfiles' filesep 'DefaultOptions.mat'],...
    'LastFolderPath','LBL_DAPI','LBL_HURP','LBL_TUBULIN','NB_FACTOR',...
    '-append');


if length(LBL_FACTOR) == 2
    [tbl,chi2,p,labels] = crosstab(LineProf.Factor1,LineProf.Factor2);
    for idx1 = 1:size(tbl,1)
        for idx2 = 1:size(tbl,2)
            disp(['Number of ' labels{idx1,1} ' - ' labels{idx2,2} ' cells: ',...
                num2str(tbl(idx1,idx2))])
        end
    end
    disp(['Chi2 test of independance: P-value = ' num2str(p) ])
else
    if length(LBL_CONDITION) == 2
        tbl = tabulate(LineProf.Factor1);
        f = cat(1,tbl{:,3});
        chi2 = (f - 50).^ 2 ./ 50;
        chi2 = sum(chi2(:));
        df = 1;        
        p = chi2cdf(chi2,df,'upper');
        disp(['Number of ' tbl{1,1} ' cells: ' num2str(tbl{1,2})])
        disp(['Number of ' tbl{2,1} ' cells: ' num2str(tbl{2,2})])
        disp(['Chi2 test of independance: P-value = ' num2str(p) ])
    end    
end
clearvars('-except','DATA')
    



%%
%1
LineProfile3D(DATA)
%2

%3
Intensity_vs_stripe_corr(DATA)
%4
LineProfile3D_Comparison(DATA)
%5

%6

%7
HURP_Intensity_ratio(DATA)
%8
Tubulin_Intensity_ratio(DATA)

%9
HURP_over_Tubulin_intensity_ratio(DATA)
%10
Depletion_efficiency(DATA)
%11
HURP_Iratio_vs_HalfSpindleratio_corr(DATA)