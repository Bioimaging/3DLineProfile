function makeSymetryXls
%makeSymetryXls Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
%
%   v1.0 14-Mar-2019 NL



path   =  uigetdir;
fnames = dir([path filesep '*.dv']);
fnames = {fnames.name}';

xlsdata(1,:) = {'file name' ,'Half-spindle ratio'};
xlsdata(2:length(fnames)+1,1) = fnames;
xlswrite([path filesep 'symmetry.xlsx'],xlsdata)

disp('Done!')


