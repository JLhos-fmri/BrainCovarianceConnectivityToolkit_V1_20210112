function BCCT_CaSCN_Mat_GUI
D.fig = figure('Name','Causal network of Structural Covaraince Connectivity(ROI-wise, volume-based)',...       
    'units','normalized',...      
    'menubar','none',...       
    'numbertitle','off',...      
    'color',[0.95 0.95 0.95],...
    'position',[0.20 0.15 0.5 0.5]);
movegui(D.fig,'center'); 

D.OutputText = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.05,0.9,0.1,0.08],...
    'style','text',...
    'string',{'Output','Directory'},...
    'Fontsize',10);
D.OutputEDIT = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.15,0.9,0.7,0.08],...
    'style','edit',...
    'string','NULL');
D.OutputSel = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.86,0.9,0.09,0.08],...
    'style','pushbutton',...
    'string','...');
%%
D.InputTEXT = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.05,0.8,0.1,0.08],...
    'style','text',...
    'string',{'Input','Directory'},...
    'Fontsize',10);
D.InputEDIT = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.15,0.8,0.7,0.08],...
    'style','edit',...
    'string','NULL');
D.InputSel = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.86,0.8,0.09,0.08],...
    'style','pushbutton',...
    'string','...');
%%
D.COVTEXT = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.05,0.7,0.1,0.08],...
    'style','radiobutton',...
    'string','COV',...
    'Fontsize',10,...
    'value',1);
D.COVEDIT = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.15,0.7,0.7,0.08],...
    'style','edit',...
    'string','NULL');
D.COVSel = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.86,0.7,0.09,0.08],...
    'style','pushbutton',...
    'string','...');

%%
D.mask = uicontrol('parent',D.fig,...
    'units','norm',...
    'pos',[0.05,0.61,0.1,0.08],...
    'style','text',...
    'Fontsize',10,...
    'string','mask');
D.maskedit = uicontrol('parent',D.fig,...
    'units','norm',...
    'pos',[0.15,0.61,0.7,0.08],...
    'style','edit',...
    'string','Defaults');
D.masksel = uicontrol('parent',D.fig,...
    'units','norm',...
    'pos',[0.86,0.61,0.09,0.08],...
    'style','pushbutton',...
    'string','...');

set(D.masksel,'callback',{@SelMASK,D});
%%
D.Seedtype1 = uibuttongroup('parent',D.fig,...
    'units','norm',...
    'pos',[0.05,0.52,0.9,0.08]);
D.Seedtype2 = uibuttongroup('parent',D.fig,...
    'units','norm',...
    'pos',[0.05,0.44,0.9,0.08]);
D.Seedtype3 = uibuttongroup('parent',D.fig,...
    'units','norm',...
    'pos',[0.05,0.36,0.9,0.08]);
D.Seedtype4 = uibuttongroup('parent',D.fig,...
    'units','norm',...
    'pos',[0.05,0.28,0.9,0.08]);

D.seedtypes(1) = uicontrol('parent',D.Seedtype1,...
    'units','norm',...
    'pos',[0.01,0.05,0.25,0.9],...
    'style','radiobutton',...
    'string','seperated ROIs',...
    'Fontsize',10,...
    'value',1);
D.seedtypes(2) = uicontrol('parent',D.Seedtype2,...
    'units','norm',...
    'pos',[0.01,0.05,0.25,0.9],...
    'style','radiobutton',...
    'string','combined ROIs',...
    'Fontsize',10,...
    'value',0);
D.seedtypes(3) = uicontrol('parent',D.Seedtype3,...
    'units','norm',...
    'pos',[0.01,0.05,0.25,0.9],...
    'style','radiobutton',...
    'string','Mat file(*.mat)',...
    'Fontsize',10,...
    'value',0);
D.seedtypes(4) = uicontrol('parent',D.Seedtype4,...
    'units','norm',...
    'pos',[0.01,0.05,0.25,0.9],...
    'style','radiobutton',...
    'string','Text File(*.txt)',...
    'Fontsize',10,...
    'value',0);

D.seedtype1_edit = uicontrol('parent',D.Seedtype1,...
    'units','norm',...
    'pos',[0.3,0.05,0.6,0.9],...
    'style','edit',...
    'string','NULL');
D.seedtype1_sel = uicontrol('parent',D.Seedtype1,...
    'units','norm',...
    'pos',[0.9,0.05,0.08,0.9],...
    'style','pushbutton',...
    'string','...',...
    'Tooltipstring','Select nifti images(one ROI one nifti).');


D.seedtype2_edit = uicontrol('parent',D.Seedtype2,...
    'units','norm',...
    'pos',[0.3,0.05,0.6,0.9],...
    'style','edit',...
    'string','NULL');
D.seedtype2_sel = uicontrol('parent',D.Seedtype2,...
    'units','norm',...
    'pos',[0.9,0.05,0.08,0.9],...
    'style','pushbutton',...
    'string','...',...
    'Tooltipstring','Select multiple ROIs in one nifti file.');

D.seedtype3_edit = uicontrol('parent',D.Seedtype3,...
    'units','norm',...
    'pos',[0.3,0.05,0.4,0.9],...
    'style','edit',...
    'string','NULL');
D.seedtype3_sel = uicontrol('parent',D.Seedtype3,...
    'units','norm',...
    'pos',[0.7,0.05,0.08,0.9],...
    'style','pushbutton',...
    'string','...');
D.seedtype3_matname = uicontrol('parent',D.Seedtype3,...
    'units','norm',...
    'pos',[0.8,0.55,0.18,0.4],...
    'style','text',...
    'string','variable name');
D.seedtype3_matnameedit = uicontrol('parent',D.Seedtype3,...
    'units','norm',...
    'pos',[0.8,0.05,0.18,0.4],...
    'style','edit');

D.seedtype4_edit = uicontrol('parent',D.Seedtype4,...
    'units','norm',...
    'pos',[0.3,0.05,0.6,0.9],...
    'style','edit',...
    'string','NULL');
D.seedtype4_sel = uicontrol('parent',D.Seedtype4,...
    'units','norm',...
    'pos',[0.9,0.05,0.08,0.9],...
    'style','pushbutton',...
    'string','...');

set(D.seedtypes(2),'val',0);
set(D.seedtype2_edit,'enable','off');
set(D.seedtype2_sel,'enable','off');
set(D.seedtypes(3),'val',0);
set(D.seedtype3_edit,'enable','off');
set(D.seedtype3_sel,'enable','off');
set(D.seedtype3_matname,'enable','off');
set(D.seedtype3_matnameedit,'enable','off');
set(D.seedtypes(4),'val',0);
set(D.seedtype4_edit,'enable','off');
set(D.seedtype4_sel,'enable','off');
%%

D.OrderEdit = uibuttongroup('parent',D.fig,...
    'units','norm',...
    'pos',[0.05,0.2,0.9,0.08]);


D.ORD1 = uicontrol('parent',D.OrderEdit,...
    'units','norm',...
    'pos',[0.05,0.05,0.1,0.9],...
    'style','text',...
    'Fontsize',10,...
    'string',{'Image', 'Orders'});
D.ORDedit1 = uicontrol('parent',D.OrderEdit,...
    'units','norm',...
    'pos',[0.15,0.05,0.3,0.9],...
    'style','edit',...
    'string','Defaults');
D.ORDsel1 = uicontrol('parent',D.OrderEdit,...
    'units','norm',...
    'pos',[0.46,0.05,0.09,0.9],...
    'style','pushbutton',...
    'string','...');
D.calord = uicontrol('parent',D.OrderEdit,...
    'units','norm',...
    'pos',[0.56,0.05,0.09,0.9],...
    'style','text',...
    'string',{'GCA', 'order'});
D.calordedit = uicontrol('parent',D.OrderEdit,...
    'units','norm',...
    'pos',[0.66,0.05,0.09,0.9],...
    'style','edit',...
    'string','1');
D.Calmethod = uibuttongroup('parent',D.OrderEdit,...
    'units','norm',...
    'pos',[0.76 0.05 0.19 0.9]);
D.Calmethod1 = uicontrol('parent',D.Calmethod,...
    'units','norm',...
    'pos',[0.01,0.01,0.48,0.98],...
    'style','rad',...
    'string','residual');
D.Calmethod2 = uicontrol('parent',D.Calmethod,...
    'units','norm',...
    'pos',[0.51,0.01,0.48,0.98],...
    'style','rad',...
    'string','Coefficient');
set(D.Calmethod1,'val',1);
set(D.Calmethod2,'val',0);

set(D.ORDsel1,'callback',{@SelImgord_I,D});
%%
D.GlobalScale = uicontrol('parent',D.fig,...
    'units','norm',...
    'style','rad',...
    'string','GlobalScale?',...
    'pos',[0.05,0.05,0.1,0.08],...
    'val',0);
D.perm  = uicontrol('parent',D.fig,...
    'units','norm',...
    'style','rad',...
    'string','permutation',...
    'pos',[0.2,0.05,0.1,0.08],...
    'val',0);
D.permnum = uicontrol('parent',D.fig,...
    'units','norm',...
    'style','edit',...
    'pos',[0.3,0.05,0.1,0.08],...
    'string','5000');

D.Comp = uicontrol('parent',D.fig,...
    'units','norm',...
    'pos',[0.6,0.05,0.15,0.08],...
    'style','pushbutton',...
    'string','Compute',...
    'Fontsize',10);
D.Exit = uicontrol('parent',D.fig,...
    'units','norm',...
    'pos',[0.8,0.05,0.15,0.08],...
    'style','pushbutton',...
    'string','Exit',...
    'Fontsize',10);
%

set(D.OutputSel,'callback',{@OutputSel,D});
set(D.InputSel,'callback',{@InputSel,D});
set(D.COVTEXT,'callback',{@COVselp,D});
set(D.COVSel,'callback',{@COVSel,D});
set(D.seedtypes(1),'callback',{@ChangeTosepROI,D});
set(D.seedtypes(2),'callback',{@ChangeToROI,D});
set(D.seedtypes(3),'callback',{@ChangeToMAT,D});
set(D.seedtypes(4),'callback',{@ChangeToTEXT,D});
set(D.seedtype1_sel,'callback',{@SelROIsep,D});
set(D.seedtype2_sel,'callback',{@SelROI,D});
set(D.seedtype3_sel,'callback',{@Selmat,D});
set(D.seedtype4_sel,'callback',{@Seltext,D});
% set(D.partial,'callback',{@PartialCom,D});
set(D.Comp,'callback',{@Compute,D});
set(D.Exit,'callback',{@Exit,D});
end
% function PartialCom(varargin);
% 
% end
function Compute(varargin)
D = varargin{3};
cal_types = 'No';
if 0
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
Parameter.GlobalScal = get(D.GlobalScale,'val');
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
Seedtype = get(D.seedtypes,'val');
if Seedtype{1}
    Parameter.SeedROItype = 1;
    pgs = get(D.seedtype1_edit,'string');
    Parameter.SeedROI_seproi = pgs;
elseif Seedtype{2}
    Parameter.SeedROItype = 2;
    roipath = get(D.seedtype2_edit,'string');
    Parameter.SeedROI_nifti = roipath;
elseif Seedtype{3}
    Parameter.SeedROItype = 3;
    matpath = get(D.seedtype3_edit,'string');
    Parameter.SeedROI_mat = matpath;
    varnam = get(D.seedtype3_matnameedit,'string');
    Parameter.SeedROI_varname = varnam;
else
    Parameter.SeedROItype = 4;
    txtpath = get(D.seedtype4_edit,'string');
    Parameter.SeedROI_txt = txtpath;
end

Parameter.imageorder = get(D.ORDedit1,'string');
Parameter.GCAorder = str2num(get(D.calordedit,'string'));
Parameter.Calmethod1 = get(D.Calmethod1,'val');

Parameter.perm = get(D.perm,'val');
Parameter.permnum = str2num(get(D.permnum,'string'));


masks = get(D.maskedit,'string');
Parameter.masks = masks;

outputmat = fullfile(Outputdir,'SetUpparameter.mat');
save(outputmat,'Parameter');
BCCT_CaSCN_matrixcomputemain(Parameter)
uiwait(msgbox('Finished'));
end
function SelROIsep(varargin)
D = varargin{3};
fold_name = uigetdir(pwd,'select the file directory which rois stored');
set(D.seedtype1_edit,'string',fold_name);
end
function SelROI(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'SeedROI selection');
PGseed = fullfile(PathName,FileName);
set(D.seedtype2_edit,'string',PGseed);
end
function Selmat(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.mat'},'mat file selection');
PGseed = fullfile(PathName,FileName);
set(D.seedtype3_edit,'string',PGseed);
end
function Seltext(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.txt'},'mat file selection');
PGseed = fullfile(PathName,FileName);
set(D.seedtype4_edit,'string',PGseed);
end
function ChangeTosepROI(varargin)
D = varargin{3};
set(D.seedtypes(1),'val',1);
set(D.seedtype1_sel,'enable','on');
set(D.seedtype1_edit,'enable','on');
set(D.seedtypes(2),'val',0);
set(D.seedtype2_edit,'enable','off');
set(D.seedtype2_sel,'enable','off');
set(D.seedtypes(3),'val',0);
set(D.seedtype3_edit,'enable','off');
set(D.seedtype3_sel,'enable','off');
set(D.seedtype3_matname,'enable','off');
set(D.seedtype3_matnameedit,'enable','off');
set(D.seedtypes(4),'val',0);
set(D.seedtype4_edit,'enable','off');
set(D.seedtype4_sel,'enable','off');
end
function ChangeToROI(varargin)
D = varargin{3};
set(D.seedtypes(1),'val',0);
set(D.seedtype1_sel,'enable','off');
set(D.seedtype1_edit,'enable','off');
set(D.seedtypes(2),'val',1);
set(D.seedtype2_edit,'enable','on');
set(D.seedtype2_sel,'enable','on');
set(D.seedtypes(3),'val',0);
set(D.seedtype3_edit,'enable','off');
set(D.seedtype3_sel,'enable','off');
set(D.seedtype3_matname,'enable','off');
set(D.seedtype3_matnameedit,'enable','off');
set(D.seedtypes(4),'val',0);
set(D.seedtype4_edit,'enable','off');
set(D.seedtype4_sel,'enable','off');
end
function ChangeToMAT(varargin)
D = varargin{3};
set(D.seedtypes(1),'val',0);
set(D.seedtype1_sel,'enable','off');
set(D.seedtype1_edit,'enable','off');
set(D.seedtypes(2),'val',0);
set(D.seedtype2_edit,'enable','off');
set(D.seedtype2_sel,'enable','off');
set(D.seedtypes(3),'val',1);
set(D.seedtype3_edit,'enable','on');
set(D.seedtype3_sel,'enable','on');
set(D.seedtype3_matname,'enable','on');
set(D.seedtype3_matnameedit,'enable','on');
set(D.seedtypes(4),'val',0);
set(D.seedtype4_edit,'enable','off');
set(D.seedtype4_sel,'enable','off');
end
function ChangeToTEXT(varargin)
D = varargin{3};
set(D.seedtypes(1),'val',0);
set(D.seedtype1_sel,'enable','off');
set(D.seedtype1_edit,'enable','off');
set(D.seedtypes(2),'val',0);
set(D.seedtype2_edit,'enable','off');
set(D.seedtype2_sel,'enable','off');
set(D.seedtypes(3),'val',0);
set(D.seedtype3_edit,'enable','off');
set(D.seedtype3_sel,'enable','off');
set(D.seedtype3_matname,'enable','off');
set(D.seedtype3_matnameedit,'enable','off');
set(D.seedtypes(4),'val',1);
set(D.seedtype4_edit,'enable','on');
set(D.seedtype4_sel,'enable','on');
end
function Exit(varargin)
D = varargin{3};
close(D.fig);
BCCT_CaSCN_GUI;
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
function SelImgord_I(varargin)
D = varargin{3};

[FileName,PathName,FilterIndex] = uigetfile({'*.txt'},'image orders');
PGseed = fullfile(PathName,FileName);
set(D.ORDedit1,'string',PGseed);
end

function SelMASK(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'maskfile selection');
PGseed = fullfile(PathName,FileName);
set(D.maskedit,'string',PGseed);
end