function BCCT_WTA_SetupCal
D.fig = figure('Name','SetUp&Cal: Brain Covariance Connectivity Cor2SubCor&Winner Take All',...       
    'units','normalized',...      
    'menubar','none',...       
    'numbertitle','off',...      
    'color',[0.95 0.95 0.95],...
    'position',[0.20 0.15 0.5 0.5]);
movegui(D.fig,'center'); 
D.OutputText = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.05,0.8,0.1,0.08],...
    'style','text',...
    'string',{'Output','Directory'},...
    'Fontsize',10);
D.OutputEDIT = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.15,0.8,0.7,0.08],...
    'style','edit',...
    'string','NULL');
D.OutputSel = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.86,0.8,0.09,0.08],...
    'style','pushbutton',...
    'string','...');
%%
D.InputTEXT = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.05,0.7,0.1,0.08],...
    'style','text',...
    'string',{'Input','Directory'},...
    'Fontsize',10);
D.InputEDIT = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.15,0.7,0.7,0.08],...
    'style','edit',...
    'string','NULL');
D.InputSel = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.86,0.7,0.09,0.08],...
    'style','pushbutton',...
    'string','...');
%%
D.SeedROIText = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.05,0.6,0.1,0.08],...
    'style','text',...
    'string',{'Seed','ROI'},...
    'Fontsize',10);
D.SeedROIEDIT = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.15,0.6,0.7,0.08],...
    'style','edit',...
    'string','NULL');
D.SeedROISel = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.86,0.6,0.09,0.08],...
    'style','pushbutton',...
    'string','...');
%%
D.TargetROITEXT = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.05,0.5,0.1,0.08],...
    'style','text',...
    'string',{'Target','ROI'},...
    'Fontsize',10);
D.TargetROIEDIT = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.15,0.5,0.7,0.08],...
    'style','edit',...
    'string','NULL');
D.TargetROISel = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.86,0.5,0.09,0.08],...
    'style','pushbutton',...
    'string','...');
%%
D.COVTEXT = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.05,0.4,0.1,0.08],...
    'style','radiobutton',...
    'string','COV',...
    'Fontsize',10,...
    'value',1);
D.COVEDIT = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.15,0.4,0.7,0.08],...
    'style','edit',...
    'string','NULL');
D.COVSel = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.86,0.4,0.09,0.08],...
    'style','pushbutton',...
    'string','...');
%%
D.methods(1) = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.05,0.3,0.15,0.08],...
    'style','radiobutton',...
    'string','Pearson correlation',...
    'Fontsize',10,...
    'value',1);
D.methods(2) = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.25,0.3,0.15,0.08],...
    'style','radiobutton',...
    'string','partial correlation',...
    'Fontsize',10,...
    'value',0);
%%
D.Compute = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.05,0.1,0.2,0.2],...
    'style','pushbutton',...
    'string','Compute',...
    'Fontsize',10);
D.ShowRes = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.3,0.1,0.2,0.2],...
    'style','pushbutton',...
    'string','Show Results',...
    'Fontsize',10,...
    'enable','off');
D.ShowResRadar = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.55,0.1,0.2,0.2],...
    'style','pushbutton',...
    'string','Show Radar figure',...
    'Fontsize',10,...
    'enable','off');
D.Exit = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.8,0.1,0.15,0.2],...
    'style','pushbutton',...
    'string','Exit',...
    'Fontsize',10);

set(D.OutputSel,'callback',{@OutputSel,D});
set(D.InputSel,'callback',{@InputSel,D});
set(D.SeedROISel,'callback',{@SeedROISel,D});
set(D.TargetROISel,'callback',{@TargetROISel,D});
set(D.COVTEXT,'callback',{@COVselp,D});
set(D.COVSel,'callback',{@COVSel,D});
set(D.methods(1),'callback',{@methodsel1,D});
set(D.methods(2),'callback',{@methodsel2,D});
set(D.Compute,'callback',{@Compute,D});
set(D.ShowRes,'callback',{@ShowRes,D});
set(D.ShowResRadar,'callback',{@ShowResRadar,D});
set(D.Exit,'callback',{@Exit,D});
end
function Exit(varargin)
D = varargin{3};
close(D.fig);
BCCT_WTA_GUI;
end
function OutputSel(varargin)
D = varargin{3};
PG = uigetdir(pwd,'Output Directory Selection');
set(D.OutputEDIT,'string',PG)
end
function InputSel(varargin)
D = varargin{3};
PG = uigetdir(pwd,'Input Directory Selection');
set(D.InputEDIT,'string',PG)
end
function SeedROISel(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'SeedROI selection');
PGseed = fullfile(PathName,FileName);
set(D.SeedROIEDIT,'string',PGseed);
end
function TargetROISel(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'TargetROI selection');
PGtarget = fullfile(PathName,FileName);
set(D.TargetROIEDIT,'string',PGtarget);
end
function COVselp(varargin)
D = varargin{3};
val = get(D.COVTEXT,'value');
if val==1
    set(D.COVSel,'Enable','on');
    set(D.COVEDIT,'Enable','on');
else
    set(D.COVSel,'Enable','off');
    set(D.COVEDIT,'Enable','off');
end
end
function COVSel(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile('*.txt','COVtext selection');
PGtxt = fullfile(PathName,FileName);
set(D.COVEDIT,'string',PGtxt);
end
function methodsel1(varargin)
D = varargin{3};
set(D.methods(1),'val',1);
set(D.methods(2),'val',0);
end
function methodsel2(varargin)
D = varargin{3};
set(D.methods(1),'val',0);
set(D.methods(2),'val',1);
end
function Compute(varargin)
D = varargin{3};
cal_types = 'No';
if 0 % consider it in master version.
    cal_types = questdlg('Using ScatterOutlier?','ScatterOutlier','No');
end
Parameter.cal_types = cal_types;
Outputdir = get(D.OutputEDIT,'string');
if strcmp(Outputdir,'NULL')
    errordlg('Please Select Outputdir');
    error('Please Select Outputdir');
end
Parameter.Outputdir = Outputdir;
Inputdir = get(D.InputEDIT,'string');
if strcmp(Inputdir,'NULL')
    errordlg('Please Select Inputdir');
    error('Please Select Inputdir');
end
Parameter.Inputdir = Inputdir;
SeedROI = get(D.SeedROIEDIT,'string');
if strcmp(SeedROI,'NULL')
    errordlg('Please Select SeedROI');
    error('Please Select SeedROI');
end
Parameter.SeedROI = SeedROI;
TargetROI = get(D.TargetROIEDIT,'string');
if strcmp(TargetROI,'NULL')
    errordlg('Please Select TargetROI');
    error('Please Select TargetROI');
end
Parameter.TargetROI = TargetROI;
covs = get(D.COVTEXT,'val');
Parameter.covs = covs;
if covs==1
    COVtext = get(D.COVEDIT,'string');
    if strcmp(COVtext,'NULL');
        errordlg('Please select covariate text file');
        error('Please select covariate text file');
    end
    Parameter.COVtext = COVtext;
else
    Parameter.COVtext = '';
end
methods = get(D.methods(1),'val');
if methods==1
    methodused = 1;
else
    methodused = 2;
end
Parameter.methodused = methodused;
outputmat = fullfile(Outputdir,'SetUpparameter.mat');
save(outputmat,'Parameter');
BCCT_WTA_computemain(Parameter)
set(D.ShowRes,'enable','on')
set(D.ShowResRadar,'enable','on')
end
function ShowRes(varargin)
D = varargin{3};
% Outputdir = get(D.OutputEDIT,'string');
% ShowResSep_in(Outputdir);
% if ~isempty(dir(fullfile(Outputdir,'Scatter_LabedVal.mat')));
%     ShowResSep_Scatter_in(Outputdir);
% end
BCCT_ShowWTAGUI;
end
function ShowResRadar(varargin)
D = varargin{3};
% Outputdir = get(D.OutputEDIT,'string');
% ShowResRadarSep_in(Outputdir)
% 
% if ~isempty(dir(fullfile(Outputdir,'Scatter_LabedVal.mat')));
%     ShowResRadarSep_Scatter_in(Outputdir);
% end
BCCT_ShowWTAGUI
end
