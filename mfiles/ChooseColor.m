function color = ChooseColor(Condition,rgb)
%   Nicolas Liaudet
%
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
%
%   v1.0 07-Feb-2019 NL

color = [];
colorfig = figure('Units','Pixel',...
    'MenuBar','none',...
    'Position',[10 10 650 325],...
    'Resize','off',...
    'WindowStyle','normal',...
    'Name',Condition,...
    'NumberTitle','off',...
    'DeleteFcn',@GetColor);
    
    cc = javax.swing.JColorChooser;
    cc.setColor(rgb(1),rgb(2),rgb(3))
    
    [jColorChooser,container] = javacomponent(cc,[1,1,650,325],colorfig);
    movegui(colorfig,'center')
    uiwait(colorfig)
    r = get(color,'Red');    
    g = get(color,'Green');
    b = get(color,'Blue');
    color = [r g b]/255;
    
    function GetColor(hObj,evt,hFigure)
        color = jColorChooser.getColor;    
    end

end