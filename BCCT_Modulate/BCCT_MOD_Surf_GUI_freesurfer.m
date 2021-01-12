function BCCT_MOD_Surf_GUI_freesurfer

D.fig = figure('Name','Modulate effect of Covaraince Connectivity(Seed to Whole Brain, Surfacce-based, (standard Freesurfer process))',...       
    'units','normalized',...      
    'menubar','none',...       
    'numbertitle','off',...      
    'color',[0.95 0.95 0.95],...
    'position',[0.20 0.15 0.5 0.5]);
movegui(D.fig,'center'); 

D.SpaceSel = uibuttongroup('parent',D.fig,...
    'units','norm',...
    'pos',[0.05,0.9,0.9,0.09]);
D.SpaceType1 = uicontrol('parent',D.SpaceSel,...
    'units','norm',...
    'pos',[0.05,0.1,0.15,0.8],...
    'style','rad',...
    'string','fsaverage');
D.SpaceType2 = uicontrol('parent',D.SpaceSel,...
    'units','norm',...
    'pos',[0.2,0.1,0.15,0.8],...
    'style','rad',...
    'string','fsaverage3');
D.SpaceType3 = uicontrol('parent',D.SpaceSel,...
    'units','norm',...
    'pos',[0.35,0.1,0.15,0.8],...
    'style','rad',...
    'string','fsaverage4');
D.SpaceType4 = uicontrol('parent',D.SpaceSel,...
    'units','norm',...
    'pos',[0.5,0.1,0.15,0.8],...
    'style','rad',...
    'string','fsaverage5');
D.SpaceType5 = uicontrol('parent',D.SpaceSel,...
    'units','norm',...
    'pos',[0.65,0.1,0.15,0.8],...
    'style','rad',...
    'string','fsaverage6');
D.SpaceType6 = uicontrol('parent',D.SpaceSel,...
    'units','norm',...
    'pos',[0.8,0.1,0.15,0.8],...
    'style','rad',...
    'string','fsaverage_sym');
%%
D.IO = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.05,0.55,0.9,0.35]);
D.outputdirtitle = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.05,0.81,0.1,0.18],...
    'style','text',...
    'string','Outputdir');
D.outputdirtext = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.2,0.81,0.6,0.18],...
    'style','edit',...
    'string','Null');
D.outputdirpb = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.85,0.81,0.1,0.18],...
    'style','pushbutton',...
    'string','...');
set(D.outputdirpb,'callback',{@OutputSel,D});

D.inputdirtitle = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.05,0.61,0.1,0.18],...
    'style','text',...
    'string','Inputdir');
D.inputdirtext = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.2,0.61,0.6,0.18],...
    'style','edit',...
    'string','Null');
D.inputdirpb = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.85,0.61,0.1,0.18],...
    'style','pushbutton',...
    'string','...');
set(D.inputdirpb,'callback',{@InputSel,D});

D.markertitle = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.05,0.41,0.1,0.18],...
    'style','text',...
    'string','markerName');
D.markertext = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.2,0.41,0.6,0.18],...
    'style','edit',...
    'string','Null');
D.markerpb = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.85,0.41,0.1,0.18],...
    'style','pushbutton',...
    'string','...');
set(D.markerpb,'callback',{@MarkerSel,D});

D.maskdirtitle = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.05,0.21,0.1,0.18],...
    'style','text',...
    'string','Maskdir');
D.maskdirtext = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.2,0.21,0.6,0.18],...
    'style','edit',...
    'string','Default');
D.maskdirpb = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.85,0.21,0.1,0.18],...
    'style','pushbutton',...
    'string','...');
set(D.maskdirpb,'callback',{@MaskSel,D});

D.covtitle1 = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.05,0.11,0.1,0.09],...
    'style','rad',...
    'string','with COV');

D.covtitle2 = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.05,0.01,0.1,0.09],...
    'style','rad',...
    'string','without COV');
D.covdirtext = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.2,0.01,0.6,0.18],...
    'style','edit',...
    'string','Null');
D.covdirpb = uicontrol('parent',D.IO,...
    'unit','norm',...
    'pos',[0.85,0.01,0.1,0.18],...
    'style','pushbutton',...
    'string','...');
set(D.covdirpb,'callback',{@CovSel,D});
%%
D.Seedtype1 = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.05,0.46,0.9,0.08]);
D.Seedtype2 = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.05,0.38,0.9,0.08]);
D.Seedtype3 = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.05,0.3,0.9,0.08]);
D.Seedtype4 = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.05,0.22,0.9,0.08]);


D.seedtypes(1) = uicontrol('parent',D.Seedtype1,...
    'units','norm',...
    'pos',[0.01,0.05,0.25,0.9],...
    'style','radiobutton',...
    'string','MNI',...
    'Fontsize',10,...
    'value',1);
D.seedtypes(2) = uicontrol('parent',D.Seedtype2,...
    'units','norm',...
    'pos',[0.01,0.05,0.25,0.9],...
    'style','radiobutton',...
    'string','Nifti Image ROIs',...
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

D.seedtype1_xtext = uicontrol('parent',D.Seedtype1,...
    'units','norm',...
    'pos',[0.3,0.3,0.05,0.4],...
    'style','text',...
    'string','x = ');
D.seedtype1_xedit = uicontrol('parent',D.Seedtype1,...
    'units','norm',...
    'pos',[0.35,0.3,0.1,0.4],...
    'style','edit');
D.seedtype1_ytext = uicontrol('parent',D.Seedtype1,...
    'units','norm',...
    'pos',[0.45,0.3,0.05,0.4],...
    'style','text',...
    'string','y = ');
D.seedtype1_yedit = uicontrol('parent',D.Seedtype1,...
    'units','norm',...
    'pos',[0.5,0.3,0.1,0.4],...
    'style','edit');
D.seedtype1_ztext = uicontrol('parent',D.Seedtype1,...
    'units','norm',...
    'pos',[0.6,0.3,0.05,0.4],...
    'style','text',...
    'string','z = ');
D.seedtype1_zedit = uicontrol('parent',D.Seedtype1,...
    'units','norm',...
    'pos',[0.65,0.3,0.1,0.4],...
    'style','edit');
% D.seedtype1_radtext = uicontrol('parent',D.Seedtype1,...
%     'units','norm',...
%     'pos',[0.75,0.3,0.1,0.4],...
%     'style','text',...
%     'string','rad(mm) = ');
% D.seedtype1_radedit = uicontrol('parent',D.Seedtype1,...
%     'units','norm',...
%     'pos',[0.85,0.3,0.1,0.4],...
%     'style','edit');

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
    'Tooltipstring','Select single ROI or multiple ROIs in one nifti file.');

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
D.Modu_Mod = uibuttongroup('parent',D.fig,...
    'units','norm',...
    'pos',[0.05,0.1,0.9,0.10]);

D.Modu_TEXT = uicontrol('parent',D.Modu_Mod,...
    'units','normalized',...
    'pos',[0.01,0.05,0.14,0.9],...
    'style','text',...
    'string',{'Modulate','Variable'},...
    'Fontsize',10);
D.Modu_EDIT = uicontrol('parent',D.Modu_Mod,...
    'units','normalized',...
    'pos',[0.15,0.05,0.85,0.9],...
    'style','edit',...
    'string','NULL');
D.Modu_Sel = uicontrol('parent',D.Modu_Mod,...
    'units','normalized',...
    'pos',[0.9,0.05,0.09,0.9],...
    'style','pushbutton',...
    'string','...');

set(D.Modu_Sel,'callback',{@ModufactSel,D});

%%
D.GlobalScale = uicontrol('parent',D.fig,...
    'unit','norm',...
    'pos',[0.05,0.05,0.25,0.05],...
    'style','rad',...
    'string','GlobalScale?',...
    'val',0);
D.Comp = uicontrol('parent',D.fig,...
    'unit','norm',...
    'pos',[0.4,0.05,0.25,0.05],...
    'style','pushbutton',...
    'string','Compute');
D.Exit = uicontrol('parent',D.fig,...
    'unit','norm',...
    'pos',[0.7,0.05,0.25,0.05],...
    'style','pushbutton',...
    'string','Exit');

set(D.Comp,'callback',{@Compute,D});
set(D.Exit,'callback',{@Exit,D});
%%
set(D.seedtypes(1),'callback',{@ChangeToMNI,D});
set(D.seedtypes(2),'callback',{@ChangeToROI,D});
set(D.seedtypes(3),'callback',{@ChangeToMAT,D});
set(D.seedtypes(4),'callback',{@ChangeToTEXT,D});
set(D.seedtype2_sel,'callback',{@SelROI,D});
set(D.seedtype3_sel,'callback',{@Selmat,D});
set(D.seedtype4_sel,'callback',{@Seltext,D});

end
function Compute(varargin)
D = varargin{3};
Outdir = get(D.outputdirtext,'string');
Parameter.Outdir = Outdir;
Indir = get(D.inputdirtext,'string');
Parameter.Indir = Indir;
FSlabs(1) = get(D.SpaceType1,'val');
FSlabs(2) = get(D.SpaceType2,'val');
FSlabs(3) = get(D.SpaceType3,'val');
FSlabs(4) = get(D.SpaceType4,'val');
FSlabs(5) = get(D.SpaceType5,'val');
FSlabs(6) = get(D.SpaceType6,'val');
Parameter.FSlabs = FSlabs;
GlobalScale = get(D.GlobalScale,'val');
Parameter.GlobalScale = GlobalScale;
covlab(1) = get(D.covtitle1,'val');
covlab(2) = get(D.covtitle2,'val');
Parameter.covlab = covlab;
masksdir = get(D.maskdirtext,'string');
Parameter.masksdir = masksdir;
markerdir = get(D.markertext,'string');
Parameter.markerdir = markerdir;
COVdir = get(D.covdirtext,'string');
Parameter.COVdir = COVdir;

Seedtype = get(D.seedtypes,'val');
if Seedtype{1}
    Parameter.SeedROItype = 1;
    tempv1 = get(D.seedtype1_xedit,'string');
    mni(1) = str2num(tempv1);
    tempv2 = get(D.seedtype1_yedit,'string');
    mni(2) = str2num(tempv2);
    tempv3 = get(D.seedtype1_zedit,'string');
    mni(3) = str2num(tempv3);
%     tempr = get(D.seedtype1_radedit,'string');
%     radius = str2num(tempr);
    Parameter.SeedROI_mni = mni;
%     Parameter.SeedROI_radius = radius;
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

Parameter.ModuFactor = get(D.Modu_EDIT,'string');
save([Outdir,filesep,'SetUpParameter.mat'],'Parameter');
BCCT_MOD_Surf_compute_freesurfer(Parameter);
end
function Exit(varargin)
D = varargin{3};
close(D.fig);
BCCT_MOD_Surf_GUI;
end

function OutputSel(varargin)
D = varargin{3};
PG = uigetdir(pwd,'Output Directory Selection');
set(D.outputdirtext,'string',PG);
end
function InputSel(varargin)
D = varargin{3};
uiwait(msgbox('Make sure that the input directory was created by freesurfer!'))
PG = uigetdir(pwd,'Input Directory Selection');
set(D.inputdirtext,'string',PG);
folds = dir(PG);
subnum = 1;
for i = 1:length(folds)-2
    if folds(i+2).isdir % make sure it's directory
        foldtemp = folds(i+2).name;
        if length(foldtemp)>=9 % fsaverage length is 9
            if ~strcmp(foldtemp(1:9),'fsaverage') % it's not fsaverage or fsaverage*
                if ~isempty(dir([PG,filesep,foldtemp,filesep,'surf'])) % contain surf fold
                    sublist{subnum,1} = foldtemp;
                    subnum = subnum+1;
                end
            end
        else
            if ~isempty(dir([PG,filesep,foldtemp,filesep,'surf'])) % contain surf fold
                sublist{subnum,1} = foldtemp;
                subnum = subnum+1;
            end
        end
    end
end
for i = 1:length(sublist)
    disp(sublist{i,1});
end
uiwait(msgbox('Please check the subjects in command window'));
end
function MaskSel(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.mgh';'*.txt'},'maskfile selection');
PGseed = fullfile(PathName,FileName);
set(D.maskdirtext,'string',PGseed);
end
function MarkerSel(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.mgh'},'marker example file selection (lh.*.mgh)');
PGseed = fullfile(PathName,FileName);
[pt,na,ex] = fileparts(PGseed);
set(D.markertext,'string',[na,ex]);
end
function CovSel(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.txt'},'Covariate of No interest');
PGseed = fullfile(PathName,FileName);
set(D.covdirtext,'string',PGseed);
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
function ChangeToMNI(varargin)
D = varargin{3};
set(D.seedtypes(1),'val',1);
set(D.seedtype1_xtext,'enable','on');
set(D.seedtype1_xedit,'enable','on');
set(D.seedtype1_ytext,'enable','on');
set(D.seedtype1_yedit,'enable','on');
set(D.seedtype1_ztext,'enable','on');
set(D.seedtype1_zedit,'enable','on');
% set(D.seedtype1_radtext,'enable','on');
% set(D.seedtype1_radedit,'enable','on');
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
set(D.seedtype1_xtext,'enable','off');
set(D.seedtype1_xedit,'enable','off');
set(D.seedtype1_ytext,'enable','off');
set(D.seedtype1_yedit,'enable','off');
set(D.seedtype1_ztext,'enable','off');
set(D.seedtype1_zedit,'enable','off');
% set(D.seedtype1_radtext,'enable','off');
% set(D.seedtype1_radedit,'enable','off');
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
set(D.seedtype1_xtext,'enable','off');
set(D.seedtype1_xedit,'enable','off');
set(D.seedtype1_ytext,'enable','off');
set(D.seedtype1_yedit,'enable','off');
set(D.seedtype1_ztext,'enable','off');
set(D.seedtype1_zedit,'enable','off');
% set(D.seedtype1_radtext,'enable','off');
% set(D.seedtype1_radedit,'enable','off');
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
set(D.seedtype1_xtext,'enable','off');
set(D.seedtype1_xedit,'enable','off');
set(D.seedtype1_ytext,'enable','off');
set(D.seedtype1_yedit,'enable','off');
set(D.seedtype1_ztext,'enable','off');
set(D.seedtype1_zedit,'enable','off');
% set(D.seedtype1_radtext,'enable','off');
% set(D.seedtype1_radedit,'enable','off');
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


function ModufactSel(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile('*.txt','Modulate Factor selection');
PGtxt = fullfile(PathName,FileName);
set(D.Modu_EDIT,'string',PGtxt);
end