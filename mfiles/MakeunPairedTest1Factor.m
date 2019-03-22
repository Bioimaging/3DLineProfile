function [ p,test_txt ] = MakeunPairedTest1Factor(x,y,alpha,varargin)
%MAKEUNPAIREDTEST1FACTOR make unpaired t-test 1 factor
%   
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   1.0 01-Feb-2018 NL


if nargin == 4
    Tail = varargin{1};
else
    Tail = 'both';
end
[h1,p1] = jbtest(x);
[h2,p2] = jbtest(y);

if ~h1&&~h2
    disp('normal distribution of the 2 samples')
     disp(['p1 = ' num2str(p1)])
     disp(['p2 = ' num2str(p2)])
     [h,p] = vartest2(x,y);
     test_txt = 'unPaired-sample t-test';
     if h == 0
         disp('Homoscedasticity: YOUPI!!!')
         disp(['p = ' num2str(p)])
         [h,p] = ttest2(x,y,'Vartype','equal','Tail',Tail,'alpha',alpha);
         if h == 0
            disp('Failure to reject the null hypothesis')
            disp(['p = ' num2str(p)])
         else
            disp('Rejection of the null hypothesis')
            disp(['p = ' num2str(p)])
         end
     else
         disp('non-Homoscedastic: c''est la vie!!!')
         disp(['p = ' num2str(p)])
         [h,p] = ttest2(x,y,'Vartype','unequal','Tail',Tail,'alpha',alpha);
         if h == 0
            disp('Failure to reject the null hypothesis')
            disp(['p = ' num2str(p)])
         else
            disp('Rejection of the null hypothesis')
            disp(['p = ' num2str(p)])
         end
     end
else
    disp('non-normal distribution of at least 1 sample')
    [p,h,stats] = ranksum(x,y,'alpha',alpha,'Tail',Tail);
    test_txt = 'Mann-Whitney U-test';
    if h == 1
        disp(['Rejection of the null hypothesis at the ' num2str(100*alpha) '% significance level.'])
        disp(['p = ' num2str(p)])
    else
        disp(['Failure to reject the null hypothesis at the ' num2str(100*alpha) '% significance level.'])
        disp(['p = ' num2str(p)])        
    end
end


end

