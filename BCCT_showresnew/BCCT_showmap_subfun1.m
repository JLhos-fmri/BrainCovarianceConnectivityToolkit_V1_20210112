function BCCT_showmap_subfun1(Parameter)
[Vbg Dbg] = Dynamic_read_dir_NIFTI(Parameter.BG);
D_BG = reshape(Dbg,Vbg.dim);
D_BG(isnan(D_BG)) = 0;
D_BG(D_BG<30) = 30;
D_BG(D_BG>85) = 85;
D_BG = (D_BG-30)/55*0.99+0.005;
if Parameter.Mod==1 % TRFZother
    [VDatShow Dat] = Dynamic_read_dir_NIFTI(Parameter.Mod1.TRFZmap);
    DatShow = reshape(Dat,VDatShow.dim);
    Intv = Parameter.Mod1.Int;
    Pval = Parameter.Mod1.Pval;
    
    DatShow_orig = dynamicBC_Reslice_forshow(Parameter.Mod1.TRFZmap,1,Parameter.BG);
    DatShow_orig(isnan(DatShow_orig)) = 0;
    DatShow_orig = reshape(DatShow_orig,Vbg.dim);    
    DatShow_F2 = DatShow_orig.*(DatShow_orig>Intv)+DatShow_orig.*(DatShow_orig<-Intv);
    
    
    ColormapOut = AFNICOLORMAP(64);
    ColormapOutPos = ColormapOut(1:32,:);
    ColormapOutNeg = ColormapOut(33:64,:);
    
    minv = min(abs(DatShow_F2(DatShow_F2>0)));
    maxv = max(abs(DatShow_F2(:)));
    
    indpos = find(DatShow_F2>0);
    indneg = find(DatShow_F2<0);
    Dall = D_BG;
    %
    if Parameter.show==0 %% both
        Dall(indpos) = 1.5+(DatShow_F2(indpos)-minv)/(maxv-minv)*0.49+0.01;
        Dall(indneg) = 1.0+(DatShow_F2(indneg)+minv)*-1/(maxv-minv)*0.49+0.01;
    elseif Parameter.show==1 %% pos only
        Dall(indpos) = 1.5+(DatShow_F2(indpos)-minv)/(maxv-minv)*0.49+0.01;
    elseif Parameter.show==2
        Dall(indneg) = 1.0+(DatShow_F2(indneg)+minv)*-1/(maxv-minv)*0.49+0.01;
    end
        
    
    
    A = rot90(squeeze(Dall(:,:,round(Vbg.dim(3)/2))),1);
%     A = squeeze(Dall(:,:,round(Vbg.dim(3)/2)));
    B = rot90(squeeze(Dall(:,round(Vbg.dim(2)/2),:)),1);
    C = rot90(squeeze(Dall(round(Vbg.dim(1)/2),:,:)),1);
    ABC(1:size(C,1),1:size(C,2)) = C;
    ABC(1:size(B,1),size(C,2)+1:size(C,2)+size(B,2)) = B;
    ABC(size(C,1)+1:size(C,1)+size(A,1),1:size(A,2)) = A;
    
    
    
elseif Parameter.Mod==2 % pmap
    
    
end


gray2 = gray(64);
COLORS = [gray2;ColormapOut];

Hsize = get(0,'ScreenSize');
Bsize = min(Hsize(3),Hsize(4));
rat1 = Bsize/size(ABC,1);
rat2 = Bsize/size(ABC,2);
minrat = min(rat1,rat2);
Figsize = [minrat*size(ABC,1),minrat*size(ABC,2)]*0.8;
con.Hfig = figure('pos',[(Hsize(3)-Figsize(1))*0.6, (Hsize(4)-Figsize(2))*0.5, Figsize(1), Figsize(2)]);
con.Hfigpos = [(Hsize(3)-Figsize(1))*0.6, (Hsize(4)-Figsize(2))*0.5, Figsize(1), Figsize(2)];
con.axes = axes('parent',con.Hfig,'unit','norm','pos',[0,0,1,1]);
axis(con.axes,'off');
con.COLORS = COLORS;
imshow(ABC,[0,2],'parent',con.axes);axis off;colormap(con.axes,COLORS);
con.fig = figure('unit','norm',...
    'pos',[0.1,0.4,0.2,0.2]);
mni_coord = cor2mni(round(Vbg.dim/2), Vbg.mat);
con.xtext = uicontrol('parent',con.fig,...
    'unit','norm',...
    'pos',[0,2/3,1/3,1/3],...
    'style','text',...
    'string',{'X',['mni: ',num2str(mni_coord(1))],['cor: ',num2str(round(Vbg.dim(1)/2))]});
con.xadd = uicontrol('parent',con.fig,...
    'unit','norm',...
    'pos',[1/3,2/3,1/6,1/3],...
    'style','pushbutton',...
    'string','Add');
con.xminus = uicontrol('parent',con.fig,...
    'unit','norm',...
    'pos',[1/2,2/3,1/6,1/3],...
    'style','pushbutton',...
    'string','Minus');

con.ytext = uicontrol('parent',con.fig,...
    'unit','norm',...
    'pos',[0,1/3,1/3,1/3],...
    'style','text',...
    'string',{'Y',['mni: ',num2str(mni_coord(2))],['cor: ',num2str(round(Vbg.dim(2)/2))]});
con.yadd = uicontrol('parent',con.fig,...
    'unit','norm',...
    'pos',[1/3,1/3,1/6,1/3],...
    'style','pushbutton',...
    'string','Add');
con.yminus = uicontrol('parent',con.fig,...
    'unit','norm',...
    'pos',[1/2,1/3,1/6,1/3],...
    'style','pushbutton',...
    'string','Minus');

con.ztext = uicontrol('parent',con.fig,...
    'unit','norm',...
    'pos',[0,0/3,1/3,1/3],...
    'style','text',...
    'string',{'Z',['mni: ',num2str(mni_coord(3))],['cor: ',num2str(round(Vbg.dim(3)/2))]});
con.zadd = uicontrol('parent',con.fig,...
    'unit','norm',...
    'pos',[1/3,0/3,1/6,1/3],...
    'style','pushbutton',...
    'string','Add');
con.zminus = uicontrol('parent',con.fig,...
    'unit','norm',...
    'pos',[1/2,0/3,1/6,1/3],...
    'style','pushbutton',...
    'string','Minus');

con.printC = uicontrol('parent',con.fig,...
    'unit','norm',...
    'pos',[2/3,2/3,1/3,1/3],...
    'style','pushbutton',...
    'string','PrintCurrent');
con.printA = uicontrol('parent',con.fig,...
    'unit','norm',...
    'pos',[2/3,1/3,1/3,1/3],...
    'style','pushbutton',...
    'string','PrintAll');
con.exit = uicontrol('parent',con.fig,...
    'unit','norm',...
    'pos',[2/3,0/3,1/3,1/3],...
    'style','pushbutton',...
    'string','exit');
set(con.xadd,'callback',{@Mod1Xadd,con,Dall,Vbg})
set(con.yadd,'callback',{@Mod1Yadd,con,Dall,Vbg})
set(con.zadd,'callback',{@Mod1Zadd,con,Dall,Vbg})

set(con.xminus,'callback',{@Mod1Xminus,con,Dall,Vbg})
set(con.yminus,'callback',{@Mod1Yminus,con,Dall,Vbg})
set(con.zminus,'callback',{@Mod1Zminus,con,Dall,Vbg})

set(con.printA,'callback',{@Mod1PrintAll,con,Dall});
set(con.printC,'callback',{@Mod1PrintCurrent,con});

set(con.exit,'callback',{@Mod1Exit,con});
end
function Mod1Exit(varargin)
con = varargin{3};
close(con.Hfig)
close(con.fig)
% BCCT_ShowMap_GUI;
end

function Mod1Xadd(varargin)
con = varargin{3};
Dall = varargin{4};
Vbg = varargin{5};
dims = Vbg.dim;
strings = get(con.xtext,'string');
mnixt = strings{2};
corxt = strings{3};
ind1 = find(mnixt==':');
ind2 = find(corxt==':');
mnicx = str2num(mnixt(ind1+1:end));

strings = get(con.ytext,'string');
mniyt = strings{2};
coryt = strings{3};
ind1 = find(mniyt==':');
ind2 = find(coryt==':');
mnicy = str2num(mniyt(ind1+1:end));

strings = get(con.ztext,'string');
mnizt = strings{2};
corzt = strings{3};
ind1 = find(mnizt==':');
ind2 = find(corzt==':');
mnicz = str2num(mnizt(ind1+1:end));

mnic = [mnicx,mnicy,mnicz];
corc = mni2cor(mnic,Vbg.mat);
mnicxnew = mnicx+1;
corcnew = mni2cor([mnicxnew,mnicy,mnicz],Vbg.mat);

if corcnew(1)<1||corcnew(1)>dims(1)-1
    uiwait(msgbox('To the data edges'));
    mnicU = mnic;
    corcU = corc;
else
    mnicU = [mnicxnew,mnicy,mnicz];
    corcU = corcnew;
end
set(con.xtext,'string',{'X',['mni: ',num2str(mnicU(1))],['cor: ',num2str(corcU(1))]});
% close(con.Hfig);
% con.Hfig = figure('pos',con.Hfigpos);

A = rot90(squeeze(Dall(:,:,corcU(3))),1);
B = rot90(squeeze(Dall(:,corcU(2),:)),1);
C = rot90(squeeze(Dall(corcU(1),:,:)),1);
ABC(1:size(C,1),1:size(C,2)) = C;
ABC(1:size(B,1),size(C,2)+1:size(C,2)+size(B,2)) = B;
ABC(size(C,1)+1:size(C,1)+size(A,1),1:size(A,2)) = A;
imshow(ABC,[0,2],'parent',con.axes);axis off;colormap(con.axes,con.COLORS);

end
function Mod1Yadd(varargin)
con = varargin{3};
Dall = varargin{4};
Vbg = varargin{5};
dims = Vbg.dim;
strings = get(con.xtext,'string');
mnixt = strings{2};
corxt = strings{3};
ind1 = find(mnixt==':');
ind2 = find(corxt==':');
mnicx = str2num(mnixt(ind1+1:end));

strings = get(con.ytext,'string');
mniyt = strings{2};
coryt = strings{3};
ind1 = find(mniyt==':');
ind2 = find(coryt==':');
mnicy = str2num(mniyt(ind1+1:end));

strings = get(con.ztext,'string');
mnizt = strings{2};
corzt = strings{3};
ind1 = find(mnizt==':');
ind2 = find(corzt==':');
mnicz = str2num(mnizt(ind1+1:end));

mnic = [mnicx,mnicy,mnicz];
corc = mni2cor(mnic,Vbg.mat);
mnicynew = mnicy+1;
corcnew = mni2cor([mnicx,mnicynew,mnicz],Vbg.mat);

if corcnew(2)<1||corcnew(2)>dims(2)-1
    uiwait(msgbox('To the data edges'));
    mnicU = mnic;
    corcU = corc;
else
    mnicU = [mnicx,mnicynew,mnicz];
    corcU = corcnew;
end
set(con.ytext,'string',{'Y',['mni: ',num2str(mnicU(2))],['cor: ',num2str(corcU(2))]});
% close(con.Hfig);
% con.Hfig = figure('pos',con.Hfigpos);

A = rot90(squeeze(Dall(:,:,corcU(3))),1);
B = rot90(squeeze(Dall(:,corcU(2),:)),1);
C = rot90(squeeze(Dall(corcU(1),:,:)),1);
ABC(1:size(C,1),1:size(C,2)) = C;
ABC(1:size(B,1),size(C,2)+1:size(C,2)+size(B,2)) = B;
ABC(size(C,1)+1:size(C,1)+size(A,1),1:size(A,2)) = A;
imshow(ABC,[0,2],'parent',con.axes);axis off;colormap(con.axes,con.COLORS);
end
function Mod1Zadd(varargin)
con = varargin{3};
Dall = varargin{4};
Vbg = varargin{5};
dims = Vbg.dim;
strings = get(con.xtext,'string');
mnixt = strings{2};
corxt = strings{3};
ind1 = find(mnixt==':');
ind2 = find(corxt==':');
mnicx = str2num(mnixt(ind1+1:end));

strings = get(con.ytext,'string');
mniyt = strings{2};
coryt = strings{3};
ind1 = find(mniyt==':');
ind2 = find(coryt==':');
mnicy = str2num(mniyt(ind1+1:end));

strings = get(con.ztext,'string');
mnizt = strings{2};
corzt = strings{3};
ind1 = find(mnizt==':');
ind2 = find(corzt==':');
mnicz = str2num(mnizt(ind1+1:end));

mnic = [mnicx,mnicy,mnicz];
corc = mni2cor(mnic,Vbg.mat);
mnicznew = mnicz+1;
corcnew = mni2cor([mnicx,mnicy,mnicznew],Vbg.mat);

if corcnew(3)<1||corcnew(3)>dims(3)-1
    uiwait(msgbox('To the data edges'));
    mnicU = mnic;
    corcU = corc;
else
    mnicU = [mnicx,mnicy,mnicznew];
    corcU = corcnew;
end
set(con.ztext,'string',{'Z',['mni: ',num2str(mnicU(3))],['cor: ',num2str(corcU(3))]});
% close(con.Hfig);
% con.Hfig = figure('pos',con.Hfigpos);

A = rot90(squeeze(Dall(:,:,corcU(3))),1);
B = rot90(squeeze(Dall(:,corcU(2),:)),1);
C = rot90(squeeze(Dall(corcU(1),:,:)),1);
ABC(1:size(C,1),1:size(C,2)) = C;
ABC(1:size(B,1),size(C,2)+1:size(C,2)+size(B,2)) = B;
ABC(size(C,1)+1:size(C,1)+size(A,1),1:size(A,2)) = A;
imshow(ABC,[0,2],'parent',con.axes);axis off;colormap(con.axes,con.COLORS);
end

function Mod1Xminus(varargin)

con = varargin{3};
Dall = varargin{4};
Vbg = varargin{5};
dims = Vbg.dim;
strings = get(con.xtext,'string');
mnixt = strings{2};
corxt = strings{3};
ind1 = find(mnixt==':');
ind2 = find(corxt==':');
mnicx = str2num(mnixt(ind1+1:end));

strings = get(con.ytext,'string');
mniyt = strings{2};
coryt = strings{3};
ind1 = find(mniyt==':');
ind2 = find(coryt==':');
mnicy = str2num(mniyt(ind1+1:end));

strings = get(con.ztext,'string');
mnizt = strings{2};
corzt = strings{3};
ind1 = find(mnizt==':');
ind2 = find(corzt==':');
mnicz = str2num(mnizt(ind1+1:end));

mnic = [mnicx,mnicy,mnicz];
corc = mni2cor(mnic,Vbg.mat);
mnicxnew = mnicx-1;
corcnew = mni2cor([mnicxnew,mnicy,mnicz],Vbg.mat);

if corcnew(1)<1||corcnew(1)>dims(1)-1
    uiwait(msgbox('To the data edges'));
    mnicU = mnic;
    corcU = corc;
else
    mnicU = [mnicxnew,mnicy,mnicz];
    corcU = corcnew;
end
set(con.xtext,'string',{'X',['mni: ',num2str(mnicU(1))],['cor: ',num2str(corcU(1))]});
% close(con.Hfig);
% con.Hfig = figure('pos',con.Hfigpos);

A = rot90(squeeze(Dall(:,:,corcU(3))),1);
B = rot90(squeeze(Dall(:,corcU(2),:)),1);
C = rot90(squeeze(Dall(corcU(1),:,:)),1);
ABC(1:size(C,1),1:size(C,2)) = C;
ABC(1:size(B,1),size(C,2)+1:size(C,2)+size(B,2)) = B;
ABC(size(C,1)+1:size(C,1)+size(A,1),1:size(A,2)) = A;
imshow(ABC,[0,2],'parent',con.axes);axis off;colormap(con.axes,con.COLORS);
end
function Mod1Yminus(varargin)
con = varargin{3};
Dall = varargin{4};
Vbg = varargin{5};
dims = Vbg.dim;
strings = get(con.xtext,'string');
mnixt = strings{2};
corxt = strings{3};
ind1 = find(mnixt==':');
ind2 = find(corxt==':');
mnicx = str2num(mnixt(ind1+1:end));

strings = get(con.ytext,'string');
mniyt = strings{2};
coryt = strings{3};
ind1 = find(mniyt==':');
ind2 = find(coryt==':');
mnicy = str2num(mniyt(ind1+1:end));

strings = get(con.ztext,'string');
mnizt = strings{2};
corzt = strings{3};
ind1 = find(mnizt==':');
ind2 = find(corzt==':');
mnicz = str2num(mnizt(ind1+1:end));

mnic = [mnicx,mnicy,mnicz];
corc = mni2cor(mnic,Vbg.mat);
mnicynew = mnicy-1;
corcnew = mni2cor([mnicx,mnicynew,mnicz],Vbg.mat);

if corcnew(2)<1||corcnew(2)>dims(2)-1
    uiwait(msgbox('To the data edges'));
    mnicU = mnic;
    corcU = corc;
else
    mnicU = [mnicx,mnicynew,mnicz];
    corcU = corcnew;
end
set(con.ytext,'string',{'Y',['mni: ',num2str(mnicU(2))],['cor: ',num2str(corcU(2))]});
% close(con.Hfig);
% con.Hfig = figure('pos',con.Hfigpos);

A = rot90(squeeze(Dall(:,:,corcU(3))),1);
B = rot90(squeeze(Dall(:,corcU(2),:)),1);
C = rot90(squeeze(Dall(corcU(1),:,:)),1);
ABC(1:size(C,1),1:size(C,2)) = C;
ABC(1:size(B,1),size(C,2)+1:size(C,2)+size(B,2)) = B;
ABC(size(C,1)+1:size(C,1)+size(A,1),1:size(A,2)) = A;
imshow(ABC,[0,2],'parent',con.axes);axis off;colormap(con.axes,con.COLORS);
end
function Mod1Zminus(varargin)
con = varargin{3};
Dall = varargin{4};
Vbg = varargin{5};
dims = Vbg.dim;
strings = get(con.xtext,'string');
mnixt = strings{2};
corxt = strings{3};
ind1 = find(mnixt==':');
ind2 = find(corxt==':');
mnicx = str2num(mnixt(ind1+1:end));

strings = get(con.ytext,'string');
mniyt = strings{2};
coryt = strings{3};
ind1 = find(mniyt==':');
ind2 = find(coryt==':');
mnicy = str2num(mniyt(ind1+1:end));

strings = get(con.ztext,'string');
mnizt = strings{2};
corzt = strings{3};
ind1 = find(mnizt==':');
ind2 = find(corzt==':');
mnicz = str2num(mnizt(ind1+1:end));

mnic = [mnicx,mnicy,mnicz];
corc = mni2cor(mnic,Vbg.mat);
mnicznew = mnicz-1;
corcnew = mni2cor([mnicx,mnicy,mnicznew],Vbg.mat);

if corcnew(3)<1||corcnew(3)>dims(3)-1
    uiwait(msgbox('To the data edges'));
    mnicU = mnic;
    corcU = corc;
else
    mnicU = [mnicx,mnicy,mnicznew];
    corcU = corcnew;
end
set(con.ztext,'string',{'Z',['mni: ',num2str(mnicU(3))],['cor: ',num2str(corcU(3))]});
% close(con.Hfig);
% con.Hfig = figure('pos',con.Hfigpos);

A = rot90(squeeze(Dall(:,:,corcU(3))),1);
B = rot90(squeeze(Dall(:,corcU(2),:)),1);
C = rot90(squeeze(Dall(corcU(1),:,:)),1);
ABC(1:size(C,1),1:size(C,2)) = C;
ABC(1:size(B,1),size(C,2)+1:size(C,2)+size(B,2)) = B;
ABC(size(C,1)+1:size(C,1)+size(A,1),1:size(A,2)) = A;
imshow(ABC,[0,2],'parent',con.axes);axis off;colormap(con.axes,con.COLORS);
end
function Mod1PrintAll(varargin)
con = varargin{3};
Dall = varargin{4};
pg = uigetdir('OutputdirForPictures');
for i = 1:size(Dall,1)
    outnametemp = fullfile(pg,['x-',num2str(i),'.tif']);
    A = rot90(squeeze(Dall(i,:,:)),1);
    h = figure('pos',[10,10,size(A,1),size(A,2)]);
    hax = axes('parent',h,'pos',[0,0,1,1]);axis(hax,'off');
    imshow(A,[0,2],'parent',hax);colormap(con.COLORS);axis off;
    print(h,'-dtiff','-r300',outnametemp);
    close(h)
end
for i = 1:size(Dall,2)
    outnametemp = fullfile(pg,['y-',num2str(i),'.tif']);    
    A = rot90(squeeze(Dall(:,i,:)),1);
    h = figure('pos',[10,10,size(A,1),size(A,2)]);
    hax = axes('parent',h,'pos',[0,0,1,1]);axis(hax,'off');
    imshow(A,[0,2],hax);colormap(con.COLORS);axis off;
    print(h,'-dtiff','-r300',outnametemp);
    close(h)
end
for i = 1:size(Dall,3)
    outnametemp = fullfile(pg,['z-',num2str(i),'.tif']);
    A = rot90(squeeze(Dall(:,:,i)),1);
    h = figure('pos',[10,10,size(A,1),size(A,2)]);
    hax = axes('parent',h,'pos',[0,0,1,1]);axis(hax,'off');
    imshow(A,[0,2],hax);colormap(con.COLORS);axis off;
    print(h,'-dtiff','-r300',outnametemp);
    close(h)    
end
end
function Mod1PrintCurrent(varargin)
con = varargin{3};
pg = uigetdir('OutputdirForPictures');
cloc = clock;
print(con.Hfig,'-dtiff','-r300',fullfile(pg,['Current',...
    num2str(cloc(1)),num2str(cloc(2)),num2str(cloc(3)),...
    num2str(cloc(4)),num2str(cloc(5)),'.tif']));
end
