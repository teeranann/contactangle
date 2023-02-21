function myOutput=checkboxtest()
% width and height of the figure
width = 500;
height = 310;

% screen size
sz = get( 0, 'ScreenSize');
% center position
x = mean( sz( [1, 3]));
y = mean( sz( [2, 4]));
% Create a figure window:
fig = uifigure('Position',[x - width/2, y - height/2, width, height]);

% Create a button group and radio buttons:
bg = uibuttongroup('Parent',fig,...
    'Position',[50 70 175 120]);
label2 = uilabel(bg,...
    'Position',[10 90 175 20],...
    'Text','Crop image at the drop?');
rb4 = uiradiobutton(bg,'text','No Crop','Position',[10 60 175 20]);
rb5 = uiradiobutton(bg,'text','Crop','Position',[10 30 175 20]);

% Create a button group and radio buttons:
bg3 = uibuttongroup('Parent',fig,...
    'Position',[50 200 175 90]);
label1 = uilabel(bg3,...
    'Position',[10 60 175 20],...
    'Text','Rotated stage:');
rb6 = uicheckbox(bg3,'text','Trim 1st frame','Position',[10 30 175 20]);

% Create a button group and radio buttons:
bg2 = uibuttongroup('Parent',fig,...
    'Position',[275 70 175 220]);
label2 = uilabel(bg2,...
    'Position',[10 190 175 20],...
    'Text','Baseline selection methods:');
rb1 = uiradiobutton(bg2,'text','Manual','Position',[10 160 175 20]);
rb2 = uiradiobutton(bg2,'text','Reflection','Position',[10 130 175 20]);
rb3 = uiradiobutton(bg2,'text','Needle','Position',[10 100 175 20]);
% % Create a check box:
% cbx = uicheckbox(fig,'Position',[55 217 102 15],...
%     'ValueChangedFcn',@(cbx,event) cBoxChanged(cbx,rb3));
% end
%
% Create the function for the ValueChangedFcn callback:
% function cBoxChanged(rb4,rb3)
% val = rb4.Value;
% if val
%     rb3.Enable = 'off';
% else
%     rb3.Enable = 'on';
% end
% end
% btn = uibutton(fig,'state',...
%                'Text', 'OK',...
%                'Position',[width/2-100-10,20, 100, 30]);

btn2 = uibutton(fig,...
    'Push',...
    'Text', 'Cancel',...
    'Position',[width/2+10,20, 100, 30],...
    'ButtonPushedFcn',@(~, ~) pushBtCancel());

btnOk = uibutton(fig,...
    'push',...
    'Text','OK',...
    'Position',[width/2-100-10,20, 100, 30],...
    'ButtonPushedFcn',@(~, ~) pushBtOk(rb6,rb4,rb5,rb1,rb2,rb3));

fig.Visible='on';

uiwait(fig);
    function pushBtOk(rb6,rb4,rb5,rb1,rb2,rb3)
        myOutput(1:6)=[rb6.Value rb4.Value rb5.Value rb1.Value rb2.Value rb3.Value];
        %%put Values into myOutput
        close(fig)
    end
    function pushBtCancel()
        close(fig)
    end

end