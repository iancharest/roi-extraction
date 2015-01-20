function roi_extraction(varargin)
% ROI_EXTRACTION toolbox v0.1
%
% SYNTAX
%   roi_extraction()
%
%
% DESCRIPTION
%
%   This toolbox allows a semi-automated roi definition procedure.
%   The graphics user interface requires that you input:
%
%   1 - the folder in which your subjects are located (e.g. /home/myproject/mridata)
%   2 - the SPM.mat of the first subject for the functional localiser from which you want to expand ROIs (e.g. /home/myproject/mridata/subject1/functionallocaliserstats/SPM.mat)
%   3 - the native space anatomical of the first subject of the study (e.g. /home/myproject/mridata/subject1/structurals/sMPRAGE.nii)
%   4 - whether you want to use contiguous expansion of the ROIs or discontiguous 
%       (discontiguous is not recommended for now and an option to add an anatomically defined inclusion mask will be ported in future releases)
%   5 - the Number of voxels that you wish to save in your binary mask (this also handles matlab regular expression e.g. linspace(20,400,20) or vectors e.g. [20 40 60 80 100])
%   6 - the activation contrast that you entered in SPM at the first level to define this ROI
%   7 - A name that you want to use for your ROI masks (e.g. FFA). The name will be used in twofolds:
%            a: a folder will be created to store binary .img masks with this name (e.g. /home/myproject/mridata/subject1/masks/FFA/)
%            b: the name will be used as a suffix for the number of voxels chosen (e.g. /home/myproject/mridata/subject1/masks/FFA/leftFFA_20.img
%                                                                                       /home/myproject/mridata/subject1/masks/FFA/leftFFA_40.img ...)
%   
%   As hinted above, the toolbox will proceed with left and right hemisphere definitions independently, one after the other.
%   The toolbox will recursively select the subjects contained in your mridata folder so it is important that your folder architecture
%   remain the same across all the subjects.
%    
% REQUIREMENTS:
%
%   - SPM 8 OR ABOVE NEEDS TO BE DEFINED IN MATLAB PATH
%
% Medical Research Council - Cognition and Brain Sciences Unit
% Ian Charest - August 2014
% 
% Last Modified by Ian Charest 03-Aug-2014 10:45:20

warnstate = warning;
warning off;


%% Check for SPM
try
    spmdir = spm('dir');
    %spm('defaults', 'fmri');
catch
    disp('Please add spm path.');
    warning(warnstate(1).state);
    return
end

%% define GUI and handle the case of a GUI callback 
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @roi_extraction_OpeningFcn, ...
                   'gui_OutputFcn',  @roi_extraction_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

gui_mainfcn(gui_State, varargin{:});


function roi_extraction_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to roi_extraction (see VARARGIN)

% Choose default command line output for roi_extraction
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes roi_extraction wait for user response (see UIRESUME)
% uiwait(handles.figure1);



function pushbuttonfolder_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.studydir = uigetdir('Please select your study folder');
try
d = dir(handles.studydir);
isub = [d(:).isdir]; % returns logical vector
handles.subjects = {d(isub).name}';
handles.subjects(ismember(handles.subjects,{'.','..'})) = [];
handles.subjects(cellfun(@isempty,strfind(handles.subjects,'CBU'))) = [];
end
guidata(hObject, handles);


% --- Executes on button press in pushbuttonspm
function pushbuttonspm_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
if isfield(handles,'studydir')
    FilterSpec = fullfile(handles.studydir,'*SPM.mat');
else
    FilterSpec = fullfile('*SPM.mat');
end
[jnk,PathName] = uigetfile(FilterSpec,'Please select the first subject''s SPM.mat');

try
    [name, subject, studydir] = getpathinfo(PathName);

    if ~(strcmp(studydir,handles.studydir))
        warndlg({'Warining, that SPM.mat seems to come from a different study';'Things could go berserk from here!'});
    end
    if ~(ismember(subject,handles.subjects))
        warndlg({'Warining, either you have not yet defined the study folder or';'that subject (SPM.mat) is not in study dir'});
    end

    handles.statsdir = name;
end
guidata(hObject, handles);


% --- Executes on button press in pushbuttonanatomy.
function pushbuttonanatomy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonanatomy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'studydir')
    FilterSpec = fullfile(handles.studydir,'*.nii');
else
    FilterSpec = pwd;
end
[jnk,PathName] = uigetfile(FilterSpec,'Please select the first subject''s native space anatomical image');

try
    [name, subject, studydir] = getpathinfo(PathName);

    if ~(strcmp(studydir,handles.studydir))
        warndlg({'Warining, that anatomy seems to come from a different study';'Things could go berserk from here!'});
    end
    if ~(ismember(subject,handles.subjects))
        warndlg({'Warining, either you have not yet defined the study folder or';'that subject anatomy is not in study dir'});
    end

    handles.anatomydir = name;
end
guidata(hObject, handles);

% --- Executes on button press in pushbuttonanatomymask.
function pushbuttonanatomymask_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonanatomymask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'studydir')
    FilterSpec = fullfile(handles.studydir,'*.nii');
else
    FilterSpec = pwd;
end
[jnk,PathName] = uigetfile(FilterSpec,'Please select the first subject''s anatomical mask to constrain the search');

try
    [name, maskfolder, studydir] = getpathinfo(PathName);

    
    handles.anatomymaskdir = fullfile(maskfolder,name);
end
guidata(hObject, handles);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
tref         = str2double(get(hObject,'String')) ;
handles.tmap = tref;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
try
    nvoxels = eval(get(hObject,'String'));
catch
    nvoxels = str2double(get(hObject,'String'));
end
handles.nvoxels = nvoxels;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MyCheckboxCallback(hObject, eventdata, handles)
other = setdiff(get(handles.uipanel1,'Children'),hObject);
for ii = 1:length(other)
    set(other(ii),'Value',get(other(ii),'Min'));
end


function checkbox1_Callback(hObject, eventdata, handles)
MyCheckboxCallback(hObject, eventdata, handles)
handles.extendmethod = 'contiguous';
guidata(hObject, handles);


function checkbox2_Callback(hObject, eventdata, handles)
MyCheckboxCallback(hObject, eventdata, handles)
handles.extendmethod = 'discontiguous';
guidata(hObject, handles);

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.roiname = get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonstart.
function pushbuttonstart_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if allfieldscomplete(handles)==1
    extractROI(handles)
else
    warndlg('You need to fill in each fields')
end



%% --------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = roi_extraction_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = hObject;
varargout{2} = eventdata;
varargout{3} = handles;


function [name, subject, studydir] = getpathinfo(filepath)

delimiters = strfind(filepath,filesep);

name = deblank(filepath(delimiters(end-1)+1:delimiters(end)-1));

subject = filepath(delimiters(end-2)+1:delimiters(end-1)-1);

studydir = filepath(1:delimiters(end-2)-1);


function answer=allfieldscomplete(handles)

% the toolbox need to have

% studydir
% subjects
% statsdir
% contrast
% anatomydir
% nvoxels

if isfield(handles,'studydir') && isfield(handles,'subjects') ...
       && isfield(handles,'statsdir') && isfield(handles,'tmap') && isfield(handles,'extendmethod') ...
       && isfield(handles,'nvoxels') && isfield(handles,'anatomydir')  && isfield(handles,'roiname')
   answer=1;
else
   answer=0;   
end

function extractROI(options)


% roiExtraction version 0.0
%
%
% http://www.mrc-cbu.cam.ac.uk/people/ian.charest
%
% by Ian Charest 


mridataPath = fullfile(options.studydir);
subjects = options.subjects;
extendmethod = options.extendmethod;
nSubjects = numel(subjects);

if isfield(options,'anatomymaskdir')
    anatomymask = options.anatomymaskdir;
else
    anatomymask = 0;
end

% which ROIs do you want to extract?
roiIdentifier = options.roiname;
	
	
% define the spm t-maps related to the ROIs you want to extract
tIdentifier = options.tmap;

if anatomymask~=0
    hemiIdentifiers = {...
        'lh'...
        'rh'...
        'bilateral'...
        };
else
    hemiIdentifiers = {...
        'lh'...
        'rh'...
        };
end
%% define the number of voxels (I start with 20 voxels and grow ROIs until I reach 1000 voxels. Modify this according to your taste and your ROIs. 

nHemi = length(hemiIdentifiers);
nVoxels=options.nvoxels;

% flags for reslicing the anatomical
flags.mean=0;
flags.which=1;
flags.wrap=[1 1 0];

options.lessIsMore=0;
options.showNeg=0;
options.threshold=0;
options.transparency=0.3;
options.title='';    


%% loop over subjects, hemispheres, and ROIs
for subI = 1:nSubjects
    
    % define thisSubject  
    thisSubject.name  = subjects{subI};
    thisSubject.root  = fullfile(mridataPath,thisSubject.name);
    thisSubject.masks = fullfile(thisSubject.root,'masks',roiIdentifier);    
    thisSubject.stats = fullfile(thisSubject.root,options.statsdir);
    thisSubject.roi   = roiIdentifier;
    
    % this will create the directory if not existent
    gotoDir(thisSubject.masks)
    
    % get headers of a stats volume (this is used to reslice anat to dims
    % of stats map
    V = spm_vol(fullfile(thisSubject.stats,'spmT_0001.img'));

    % find the fname
    thisSubject.tmap=V.fname;
    
    % define the subject's structural
    thisSubject.structural = get_files(fullfile(thisSubject.root,options.anatomydir),'s*.nii');
    
    thisSubjectWitness = fullfile(thisSubject.masks,['right' roiIdentifier '_allpages.ps']);
    if (~exist(thisSubjectWitness,'file')) 
        fprintf('***\t now working on %s for %s \t***\n',thisSubject.name,thisSubject.roi)
        % reslice the anat to fit the T-map
        P=char(thisSubject.tmap,thisSubject.structural);
        spm_reslice(P,flags)     
 
        % now we have a structural image that fits the t-maps (64x64x32)
        thisSubject.reslicedStructural = get_files(fullfile(thisSubject.root,'structurals'),'rs*.nii');
        thisSubject.anatMap            = spm_read_vols(spm_vol(thisSubject.reslicedStructural));
        thisSubject.extendmethod       = extendmethod;
        % perform the roi extraction
        brainVol   = cell(1,length(hemiIdentifiers));
        for hemiI = 1:nHemi
            if anatomymask~=0
                thisSubject.anatomicalMask = fullfile(thisSubject.root,anatomymask,[hemiIdentifiers{hemiI} '.' roiIdentifier '.nii']);
            end
            thisSubject.hemisphere = hemiIdentifiers{hemiI};
            thisSubject.tMap       = tIdentifier;
            roiExtraction(thisSubject,nVoxels);
        end % hemiI
    else
        fprintf('***\t work already done for %s for %s \t***\n',thisSubject.name,thisSubject.roi)
    end
end% subI
return


% local FUNCTION: brainVol = roiExtraction(userOptions,nVoxels)
% This is the core function of the toolbox.
% % USAGE
%           brainVol = roiExtraction(userOptions,nVoxels)
%
% FUNCTION
%           this function will write down a given number of binary
%           maps (defined by nVoxels) to the location of your choice
%           as defined by userOptions.masks 
%           It works with some automatised functionality as well as visual
%           inspection of the t-maps that define your ROIs.
%           It allows navigation using the spm visualisation tools on the
%           subject's native space anatomical map.
%           
%
% ARGUMENTS
% userOptions.root       the root folder location for that given subject
%                        (absoulte path passed as a string)
%           
% userOptions.masks      the absolute folder location where you want the
%                        masks for this subject to be written
%
% userOptions.stats      the absolution path and filename of the t-map
%                        which you want to display on the subject's anatomical map
%                        for the roi extraction
% 
% userOptions.structural the absolution path and filename of the native
%                        space anatomical map 
% 
% userOptions.roi        for each ROI that you want to define, this program
%                        first tells you what to look for (i.e. which
%                        functional ROI)
%
% userOptions.hemisphere the program also specifies which hemisphere it is
%                        currently working on
%
% nVoxels                a vector of scalars defining the number of voxels
%                        you want to include in the resulting binary mask.
%                           
% Written by Ian Charest - MRC-CBSU - Updated - May 2012
%
% Thanks to the developers of SPM, a neuroimaging software that is needed
% for this tool to work. Some of the code was also inspired from the SPM
% code. 
% Please visit 
% www.fil.ion.ucl.ac.uk/spm/software/spm8/
% 
% Also thanks to Xu Cui for developing xjview which inspired some of this
% code
%
% Please visit 
% www.alivelearn.net/xjview8/

function outputmasks = roiExtraction(userOptions,nVoxels)

spmdir = spm('dir');
try
    spm('defaults', 'fmri');
catch
    [];
end

thisSubjectMaskLocation = userOptions.masks;
gotoDir(thisSubjectMaskLocation)

anatMap        = userOptions.anatMap;
thisLocStats   = userOptions.stats;
thisStruct     = userOptions.structural;
thisHemisphere = userOptions.hemisphere;
thisRoi        = userOptions.roi;
structimg      = thisStruct;
extendmethod   = userOptions.extendmethod;

if isfield(userOptions,'anatomicalMask')
    anatomymask = userOptions.anatomicalMask;
else
    anatomymask = 0;
end

if strcmp(extendmethod,'contiguous');
    contiguous=1;
else
    contiguous=0;
end

cd(thisLocStats)
message        = {['navigate to ' thisHemisphere,' ',thisRoi];'you can use spm right-click controls in glass brain';'press pause when you are on the spot to start growing rois'};
uiwait(warndlg(message,'Instructions'));

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Stats: Results');
xSPM = getxSPM(userOptions);
[SPM,xSPM]=spm_getSPM(xSPM);
M = SPM.xVol.M;
DIM =  SPM.xVol.DIM;
options.figI = Finter;
%-Space units
%----------------------------------------------------------------------
try
    try
        units = SPM.xVol.units;
    catch
        units = xSPM.units;
    end
catch
    try
        if strcmp(spm('CheckModality'),'EEG')
            datatype = {...
                'Volumetric (2D/3D)',...
                'Scalp-Time',...
                'Scalp-Frequency',...
                'Time-Frequency',...
                'Frequency-Frequency'};
            selected = spm_input('Data Type: ','+1','m',datatype);
            datatype = datatype{selected};
        else
            datatype = 'Volumetric (2D/3D)';
        end
    catch
        datatype     = 'Volumetric (2D/3D)';
    end
    
    switch datatype
        case 'Volumetric (2D/3D)'
            units    = {'mm' 'mm' 'mm'};
        case 'Scalp-Time'
            units    = {'mm' 'mm' 'ms'};
        case 'Scalp-Frequency'
            units    = {'mm' 'mm' 'Hz'};
        case 'Time-Frequency'
            units    = {'Hz' 'ms' ''};
        case 'Frequency-Frequency'
            units    = {'Hz' 'Hz' ''};
        otherwise
            error('Unknown data type.');
    end
end
if DIM(3) == 1, units{3} = ''; end
xSPM.units      = units;
SPM.xVol.units  = units;


spm_clf(Finter);
spm('FigName',['SPM{',xSPM.STAT,'}: Results'],Finter,CmdLine);
hReg  = spm_results_ui('SetupGUI',M,DIM,xSPM,Finter);


hMIPax = axes('Parent',Fgraph,'Position',[0.05 0.60 0.55 0.36],'Visible','off');
hMIPax = spm_mip_ui(xSPM.Z,xSPM.XYZmm,M,DIM,hMIPax,units);
spm_XYZreg('XReg',hReg,hMIPax,'spm_mip_ui')
hold on

spm_sections_noClear(xSPM,hReg,structimg)

% this will wait for space bar input
pause

% get the xyz vox location
xyz_vox = mm_to_vox(hReg,xSPM);
if userOptions.tMap<10
    spmT_map=fullfile(thisLocStats,['spmT_000' num2str(userOptions.tMap)]);
elseif userOptions.tMap<100
    spmT_map=fullfile(thisLocStats,['spmT_00'  num2str(userOptions.tMap)]);
end

roiName = [thisHemisphere,thisRoi];

%% flow chart
filename = fullfile([thisSubjectMaskLocation,filesep,roiName,'_allpages']);
V = spm_vol([spmT_map,'.img']);
Y = spm_read_vols(V);

if strcmp(thisHemisphere,'bilateral')
   nVoxels = [nVoxels nVoxels(end)+nVoxels(2) - nVoxels(1):nVoxels(2)-nVoxels(1):nVoxels(end)*2];
end

nRoisToCreate = numel(nVoxels);

spm_progress_bar('Init',nRoisToCreate,'Expanding Rois','Rois Complete');

for roiI = 1:nRoisToCreate
    
    % call this for contiguous voxels:
    if anatomymask~=0
        
        maskV = spm_vol(anatomymask);
        maskData = spm_read_vols(maskV);
        
        if contiguous==1
            newRoi     = resizeRoi(single(xyz_vox)',Y, nVoxels(roiI),maskData);
        else
            % call this for discontiguous voxels:
            newRoi     = discontiguousRoi(Y,nVoxels(roiI),maskData);
        end
        
    else
        if contiguous==1
            newRoi     = resizeRoi(single(xyz_vox)',Y, nVoxels(roiI));
        else
            % call this for discontiguous voxels:
            newRoi     = discontiguousRoi(Y,nVoxels(roiI));
        end
    end
    
    newRoimask = roi2mask(newRoi,size(Y));
    
    % define options
    options.title   = [roiName,'_',num2str(nVoxels(roiI)),'_voxels'];
    options.nVoxels = nVoxels(roiI);    
    Fgraph = spm_figure('GetWin','Graphics');
    options.figI = Fgraph;
    showRoiOnAnatomy(newRoimask.*Y,map2vol(anatMap),options)
    spm_print(filename)
    % update metaData structure
    newRoiMapMetadataStruct = V;
    newRoiMapMetadataStruct.fname = fullfile([thisSubjectMaskLocation,filesep,roiName,'_',num2str(nVoxels(roiI)),'.img']);
    newRoiMapMetadataStruct.descrip = ['binary ROI mask for ' roiName];
    newRoiMapMetadataStruct.dim = size(newRoimask);
    
    spm_write_vol(newRoiMapMetadataStruct, newRoimask);
    
    spm_progress_bar('Set',roiI); % update the progress bar
end
return


% local FUNCTION: getxSPM(userOptions)
% Thanks to Alex Woolgar for providing this very useful function.
function xSPM = getxSPM(userOptions)

% set up xSPM for visual localisation of various visual regions in SPM (to 
% save on GUI interaction)
% Alex Woolgar 22 Feb 2012
%
% Ian Charest Modified it to use the last digit of the passed t-map
% identifier

xSPM.swd       = userOptions.stats; % dir with SPM.mat
xSPM.title     = [userOptions.hemisphere ' ' userOptions.roi];
xSPM.Ic        = userOptions.tMap;
xSPM.n         = 1;
xSPM.Im        = []; %no masking contrast
xSPM.pm        = [];
xSPM.Ex        = [];
xSPM.u         = 1.69; %approx p < 0.01 (uncor)
xSPM.k         = 0;
xSPM.thresDesc = 'none';
return

% local FUNCTION: gotoDir(varargin)
% Thanks to Nikolaus Kriegeskorte and Cai Wingfield for providing this very useful function.
function gotoDir(varargin)

% gotoDir(path, dir)
%     Goes to path; then goes to dir, making it if necessary.
%
% gotoDir(path)
%     Goes to path, making all required directories on the way.

switch nargin
	case 2
		path = varargin{1};
		dir = varargin{2};

		try
			cd(path)
		catch
			error('gotoDir:NonexistantPath', ['The path "' path '" does not refer!']);
		end%try
		
		try
			cd(dir);
		catch
			fprintf(['The directory "' dir '" doesn''t exist at "' path '"; making it.\n']);
			mkdir(dir);
			cd(dir);
		end%try

	case 1

		path = varargin{1};
		sIndices = strfind(path, filesep);
		
		for i = 1:numel(sIndices)
		
			if i == 1 && sIndices(i) == 1
				continue;
			end%if
			
			try
				cd(path(1:sIndices(i)-1));
			catch
				fprintf(['The directory "' path(1:sIndices(i)-1) '" doesn''t exist... making it.\n']);
				mkdir(path(1:sIndices(i)-1));
                cd(path(1:sIndices(i)-1));
			end%try
		
		end%for:i
        
		% cleanup final directory!
        try
            cd(path);
        catch
			fprintf(['The directory "' path '" doesn''t exist... making it.\n']);
            mkdir(path);
        end%try
		
		cd(path);

otherwise
	error('gotoDir:BadNargin', 'Only 1 or 2 arguments allowed.');
end%switch:nargin
return

% local FUNCTION: files = get_files(direc, filt)
% Thanks to Russell Thompson's get_files function from which this is a
% modified version. modified by Ian Charest
function files = get_files(direc, filt)
% =========================================================================
% return a list of files
% filt = filter string
% direc = cell array of directory names
% revised 07-2011 Ian Charest
if nargin~=2, error('get_files:missing inputs, Please input folder(s) and file filter.'); end%if
files = [];
if ischar(direc) % if direc is already a character array
    currDir = direc;
    tmp = dir(fullfile(currDir,filt)); % find all files matching f*.nii
    tmp = [repmat([currDir filesep],size(tmp,1),1) char(tmp.name)]; % build the full path name for these files
    files = char(files,tmp);
else % if direc is a cell array
    if size(direc,1)>size(direc,2)
        nRuns=size(direc,1);
    else
        nRuns=size(direc,2);
    end
    for runI=1:nRuns % loop through each EPI session
        currDir = char(direc{runI});
        tmp = dir(fullfile(currDir,filt)); % find all files matching f*.nii
        tmp = [repmat([currDir filesep],size(tmp,1),1) char(tmp.name)]; % build the full path name for these files
        files = char(files,tmp);
    end
end
files = files(~all(files'==' ')',:);
return

% local FUNCTION:  spm_sections_noClear(xSPM,hReg,img)
% Thanks to John Ashburner and all the genius people at the FIL.
function spm_sections_noClear(xSPM,hReg,img)
%
% modified version of spm_sections that inhibs the clearing of the Fgraph.
% Rendering of regional effects [SPM{.}] on orthogonal sections
% FORMAT spm_sections(xSPM,hReg,img)
%
% xSPM  - structure containing details of excursion set (see spm_getSPM)
% hReg  - handle of MIP register
% img   - filename of background image
%__________________________________________________________________________
%
% spm_sections is called by spm_results_ui and uses variable img to
% create three orthogonal sections through a background image.
% Regional foci from the selected xSPM are rendered on this image.
%__________________________________________________________________________
% Copyright (C) 1994-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_sections.m 4199 2011-02-10 20:07:17Z guillaume $
%
% modified by 
% Ian Charest - MRC-CBSU - July 2011 - 
% so that the glass brain and anatomical 
% is not cleared with the call to spm_results_ui('Clear',Fgraph)


if ~nargin, [SPM,xSPM] = spm_getSPM; end
if nargin < 2, hReg = []; end
if nargin < 3 || isempty(img)
    [img, sts] = spm_select(1,'image','Select image for rendering on');
    if ~sts, return; end
end

Fgraph = spm_figure('GetWin','Graphics');
%spm_results_ui('Clear',Fgraph);
spm_orthviews('Reset');

global st prevsect
st.Space = spm_matrix([0 0 0  0 0 -pi/2]) * st.Space;
prevsect = img;

h = spm_orthviews('Image', img, [0.05 0.05 0.9 0.45]);
spm_orthviews('AddContext', h); 
spm_orthviews('MaxBB');
if ~isempty(hReg), spm_orthviews('Register', hReg); end
spm_orthviews('AddBlobs', h, xSPM.XYZ, xSPM.Z, xSPM.M);
spm_orthviews('Redraw');
return


% local FUNCTION: newRoi=resizeRoi(roi, map, nVox, mask)
% Thanks to Nikolaus Kriegeskorte for providing this very useful function.
function newRoi=resizeRoi(roi, map, nVox, mask)
% USAGE
%           newRoi=resizeRoi(roi, map, nVox[, mask])
%
% FUNCTION
%           redefine a roi to make it conform to the size
%           given by nVox. the new roi is defined by a region
%           growing process, which (1) is seeded at the voxel
%           that has the maximal statistical parameter within
%           the passed statistical map and (2) is prioritized
%           by the map's values.
%
% ARGUMENTS
% roi       (r)egion (o)f (i)nterest
%           a matrix of voxel positions
%           each row contains ONE-BASED coordinates (x, y, z) of a voxel.
%
% map       a 3D statistical-parameter map
%           the map must match the volume, relative to which
%           the roi-voxel coords are specified in roi.
%
% nVox      number of voxels the resized roi is to have
%
% [mask]    optional binary mask. if present, the region growing is
%           restricted to the nonzero entries of it.

% PARAMETERS
%maxNlayers=2;   %defines how many complete layers are maximally added
if ~exist('mask','var') || (exist('mask','var') && isempty(mask))
    mask=ones(size(map));
end
map(~mask)=min(map(:));
% DEFINE THE VOLUME
vol=zeros(size(map));
% FIND THE SEED (A MAXIMAL MAP VALUE IN ROI&mask)
mapINDs=sub2ind(size(vol),roi(:,1),roi(:,2),roi(:,3)); %single indices to MAP specifying voxels in the roi
roimap=map(mapINDs);                                      %column vector of statistical-map subset for the roi
[roimax,roimax_roimapIND]=max(roimap);                    %the maximal statistical map value in the roi and its index within roimap
seed_mapIND=mapINDs(roimax_roimapIND);                    %seed index within map
newRoi=[];
if nVox==0; return; end;
if mask(seed_mapIND)
    vol(seed_mapIND)=1;
    [x,y,z]=ind2sub(size(vol),seed_mapIND);
    newRoi=[newRoi;[x,y,z]];
end
if nVox==1 || isempty(newRoi)
    return;
end
% GROW THE REGION
for i=2:nVox    
    % DEFINE THE FRINGE
    cFringe=vol;
    [ivolx,ivoly,ivolz]=ind2sub(size(vol),find(vol));
    superset=[ivolx-1,ivoly,ivolz;
        ivolx+1,ivoly,ivolz;
        ivolx,ivoly-1,ivolz;
        ivolx,ivoly+1,ivolz;
        ivolx,ivoly,ivolz-1;
        ivolx,ivoly,ivolz+1];
    % exclude out-of-volume voxels
    outgrowths = superset(:,1)<1 | superset(:,2)<1 | superset(:,3)<1 | ...
        superset(:,1)>size(vol,1) | superset(:,2)>size(vol,2) | superset(:,3)>size(vol,3);
    superset(find(outgrowths),:)=[];
    % draw the layer (excluding multiply defined voxels)
    cFringe(sub2ind(size(vol),superset(:,1),superset(:,2),superset(:,3)))=1;
    cFringe=cFringe&mask;
    cFringe=cFringe-vol;    
    if size(find(cFringe),1)==0
        break; % exit the loop (possible cause of empty fringe: the whole volume is full)
    end    
    % FIND A MAXIMAL-MAP-VALUE FRINGE VOXEL...
    mapINDs=find(cFringe);                                    %single indices to MAP specifying voxels in the fringe
    fringemap=map(mapINDs);                                   %column vector of statistical-map subset for the fringe
    [fringemax,fringemax_fringemapIND]=max(fringemap);        %the maximal statistical map value in the roi and its index within roimap
    fringemax_mapIND=mapINDs(fringemax_fringemapIND);         %seed index within map    
    % ...INCLUDE IT
    vol(fringemax_mapIND)=1;
    [x,y,z]=ind2sub(size(vol),fringemax_mapIND);
    newRoi=[newRoi;[x,y,z]];    
end
return

% local FUNCTION: newRoi=discontiguousRoi(prioritizationMap,nVox,mask)
% Thanks to Nikolaus Kriegeskorte for providing this very useful function.
function roi=discontiguousRoi(prioritizationMap,nVox,mask)

mask(mask~=0)=1; % set all nonzero entries to one (mask with continuous values are correctly handled (only 0 counts as vacant)

if nVox>sum(mask(:))
    disp('Function discontiguousRoi, WARNING : the mask passed does not allow an ROI of the requested size. All mask voxels are returned as the ROI.')
    roi=mask2roi(mask);
else % mask has at least nVox voxels
    prioritizationMap(~mask)=min(prioritizationMap(:)); % sets voxels not covered by mask to zero.
    
    % retrieve the threshold value for the nVox maximum entries inside mask.
    threshold=autothresholdMap(prioritizationMap,nVox); 
    roiMask=prioritizationMap>=threshold & mask;
    roi=mask2roi(roiMask);
end
return


% local FUNCTION: mask=roi2mask(roi,volSize_vox)
% Thanks to Nikolaus Kriegeskorte for providing this very useful function.

% roi2mask is a function with two arguments:
%  ARGUMENTS:
%    roi: this is a 3xn matrix where each row contains the coordinates for
%         a point inside the roi
%    volSize_vox: this is a 1x3 vector containing the dimensions of the
%                 scanned volume.  E.g., [32 64 64]
%  RETURNS:
%    mask: a volume of size volSize_vox which is all 0s, except for the
%          points indicated by roi, which are 1s.

function mask=roi2mask(roi,volSize_vox)
roi_INDs=sub2ind(volSize_vox,roi(:,1),roi(:,2),roi(:,3)); %single indices to MAP specifying voxels in the roi
mask=false(volSize_vox);
mask(roi_INDs)=true;
return

% local FUNCTION: mask=mask2roi(mask)
% Thanks to Nikolaus Kriegeskorte for providing this very useful function.

function roi=mask2roi(mask)

[x,y,z]=ind2sub(size(mask),find(mask));
roi=[x,y,z];
return

% local FUNCTION: xyz_vox = mm_to_vox(hReg,xSPM)
% Thanks to Johan Carlin for providing this very useful function.

% Returns the voxel coordinates for the current position in the
% SPM results viewer.
% J Carlin 10/6/2011

%function xyz_vox = mm_to_vox(hReg,xSPM)
function xyz_vox = mm_to_vox(hReg,xSPM)

%hReg = findobj('Tag','hReg'); % get results figure handle
xyz_mm = spm_XYZreg('GetCoords',hReg); % get mm coordinates for current locationassignin('base','xyz_mm',xyz_mm); % puts xyz_mm in base workspace
xyz_vox = xSPM.iM*[xyz_mm;1]; %evalin('base','xSPM.iM*[xyz_mm;1]'); % avoid passing xSPM explicitly
xyz_vox = xyz_vox(1:3);
return


% local FUNCTION: vol=map2vol(map)
% Thanks to Nikolaus Kriegeskorte for providing this very useful function.
function vol=map2vol(map)

if strcmp(class(map),'single') || strcmp(class(map),'double')
    % scale into RGB range [0,1] for display
    map=map-min(map(:)); % move min to zero
    map=map/max(map(:)); % scale max to one
    vol=double(permute(repmat(map,[1 1 1 3]),[1 2 4 3]));
else
    % make vol same class as map (e.g. logical)
    vol=permute(repmat(map,[1 1 1 3]),[1 2 4 3]);
end
return


% local FUNCTION: showRoiOnAnatomy(roiMap,anatVol,options)
% Inspired from Nikolaus Kriegeskorte's function showMapOnAnatomy, 
% thanks for providing this very useful function.
% written by Ian Charest - MRC-CBSU - May 2012
function showRoiOnAnatomy(roiMap,anatVol,options)


options=setIfUnset(options,'lessIsMore',0);
options=setIfUnset(options,'showNeg',0);
options=setIfUnset(options,'threshold',0);
options=setIfUnset(options,'transparency',0.3);
options=setIfUnset(options,'title','');
options=setIfUnset(options,'figI',51);
extremeVal=max(roiMap(:));

lessIsMore=options.lessIsMore;
showNeg=options.showNeg;
threshold=options.threshold;
transparency=options.transparency;
nVoxels=options.nVoxels;
figI=options.figI;

title={options.title;[' (',num2str(nVoxels),' voxels marked)']};
vol=addStatMapToVol(anatVol, roiMap, threshold, extremeVal, lessIsMore, showNeg, transparency);

%showVol(vol,[title,[' (',num2str(nMarkedVoxels),' marked)']],figI,'right');

showVol(vol,title,figI,'right');
return

% local FUNCTION: options=setIfUnset(options,field,value)
function options=setIfUnset(options,field,value)
% if options.(field) is empty or doesn't exist, this function sets options.(field) to value.

if ~isfield(options, field) || isempty(options.(field))
       options.(field)=value;
end
return

% local FUNCTION: vol=addStatMapToVol(vol, map, critVal, extremeVal,lessIsMore, showNeg, transparency)
% Thanks to Nikolaus Kriegeskorte for providing this function
function vol=addStatMapToVol(vol, map, critVal, extremeVal, lessIsMore, showNeg, transparency)

% USAGE:        vol=addStatMapToVol(vol, map, [critVal, extremeVal, lessIsMore=0, showNeg=1, transparency=0])
%
% FUNCTION:     to superimpose the statistical map "map" to the true-color
%               volume "vol".
%
% ARGUMENTS:
% vol           the volume as a stack of true-color slices: X by Y by 3 by Z.
%               the third dimension encodes the color component (red, green,
%               blue).
%               anatomically, the X axis points to the left, the Y axis to
%               the back of the brain, and the Z axis up.
%               if vol contains a scalar zero, the map itself is used as a
%               grayscale background to its colored peaks.
%
% map           a statistical map as a 3D array of double-precision floats
%
% critVal       threshold determining what part of the map is superimposed:
%               only voxels whose ABSOLUTE map value exceeds critVal are
%               marked
%
% extremeVal    absolute value of the map entries above which the color
%               scale does not differentiate anymore
%
% [lessIsMore=0]  if this optional argument is nonzero, then instead of
%               marking voxels whose absolute map value EXCEEDS the
%               threshold, the function marks the voxels whose absolute map
%               value IS SMALLER than the threshold. this should be nonzero
%               for p maps.
%
% [showNeg=1]   if this optional argument is nonzero or missing (defaults to 1),
%               negative map values whose absolute value exceeds the
%               threshold are highlighted. if showNeg is zero, only
%               positive values exceeding the threshold are highlighted.
%
% RETURN VALUE: true-color volume (X by Y by 3 by Z) with the thresholded
%               map superimposed
%
% GENERAL IDEA: Given a volume of colour information and a statistical map
%               (for example, r values), this function will overwrite the
%               colour information in the original volume at the voxels
%               where the statistics in the map are significant.  In these
%               voxels, the colour will reflect the significance level, but
%               can be capped using extremeVal.  So, for statistics in the
%               map which are less than critVal, no colour will be changed.
%               For values between critVal and extremeVal, the colour will
%               be given by the statistics; and for values greater than
%               extremeVal, the colour will be at the maximum but
%               unchanging.  Finally, if lessIsMore is used, the relations
%               are reversed (since, for example using p values, to be
%               significant is to be *less* than a given threshold).




%% RESCALE THE STATISTICAL MAP
if ~exist('critVal','var')
    critVal=0;
end

if ~exist('extremeVal','var')
    extremeVal=max(abs(map(:)));
end

if ~exist('lessIsMore','var')
    lessIsMore=0;
end

if ~exist('showNeg','var')
    showNeg=1;
end

if ~exist('transparency','var')
    transparency=0;
end

if numel(vol)==1 && vol==0
    vol=map2vol(map);
end

if lessIsMore
    % flip positive and negative part of the map value axis
    % such that the maximum and minimum fall on zero
    % and the extreme values (positive and negative)
    % correspond to small values in the original map
    mx=max(abs(map(:)));
    mn=-mx;
    map(map>=0)=mx-map(map>=0);
    map(map<0)=mn-map(map<0);
    
    critVal=mx-critVal;
    extremeVal=mx-extremeVal;
end


NaN_flags=isnan(map);
if sum(NaN_flags(:))
    disp('NaNs were found in the map and set to 0.')
    map(NaN_flags)=0;
end


imagemap=double(map);

%shift to-be-shown portions (abs(map) between critVal and extremeVal) to span
%[1,31] (negative map values) and [33,64] (positive map values) 
imagemap(map>critVal)=(imagemap(map>critVal)-critVal)/(extremeVal-critVal)*(64-33)+33;
imagemap(map<-critVal)=(imagemap(map<-critVal)+critVal)/(extremeVal-critVal)*(31-1)+31;
%figure(iFig-2); clf; imagesc(imagemap(:,:,5));

imagemap(abs(map)<=critVal)=32;
imagemap(map>extremeVal)=64;
imagemap(map<-extremeVal)=1;
%figure(iFig-1); clf; imagesc(imagemap(:,:,5));

if ~showNeg
    imagemap(map<0)=32;
end


%% DEFINE THE COLORMAP
bgc=[0 0 0]; %background color
critposc=[1 0.3 0]; %color at critical positive value
critnegc=[0 0.3 1]; %color at critical negative value
maxposc=[1 1 0]; %color at maximum positive value
maxnegc=[0 1 0.3]; %color at maximum negative value

for i=1:31
    distToMin=(i-1)/30;
    distToMid=(31-i)/30;
    cm(i,:)=distToMin*critnegc+distToMid*maxnegc;
end

cm(32,:)=bgc;

for i=33:64
    distToMax=(64-i)/31;
    distToMid=(i-33)/31;
    cm(i,:)=distToMax*critposc+distToMid*maxposc;
end


%% MARK THE VOXELS
toBeMarked_INDs=find(imagemap~=32);

volXYZ3=permute(vol,[1 2 4 3]);

  redmap=volXYZ3(:,:,:,1);
greenmap=volXYZ3(:,:,:,2);
 bluemap=volXYZ3(:,:,:,3);

  redmap(toBeMarked_INDs)=  redmap(toBeMarked_INDs)*transparency + (1-transparency)*cm(round(imagemap(toBeMarked_INDs)),1);
greenmap(toBeMarked_INDs)=greenmap(toBeMarked_INDs)*transparency + (1-transparency)*cm(round(imagemap(toBeMarked_INDs)),2);
 bluemap(toBeMarked_INDs)= bluemap(toBeMarked_INDs)*transparency + (1-transparency)*cm(round(imagemap(toBeMarked_INDs)),3);

volXYZ3(:,:,:,1)=  redmap;
volXYZ3(:,:,:,2)=greenmap;
volXYZ3(:,:,:,3)= bluemap;

vol=permute(volXYZ3,[1 2 4 3]);
return


% local FUNCTION: showVol(volORmap, title, figI, right, skipNslices,nHorPanels, nVerPanels, labels, labelPos)
% Thanks to Nikolaus Kriegeskorte for providing this function.
function showVol(volORmap, title, figI, right, skipNslices, nHorPanels, nVerPanels, labels, labelPos)

% USAGE:        showVol(vol, [title, figI, right, skipNslices, nHorPanels, nVerPanels, labels, labelPos])
%
% FUNCTION:     to display a volume as a set of slices
%
% ARGUMENTS:
% vol           the volume as a stack of true-color slices:
%               X by Y by 3 by Z. the third dimension encodes 
%               the color component (red, green, blue).
%               anatomically, the X axis points to the left,
%               the Y axis to the back of the brain, and the
%               Z axis up.
%
% [title]       string to be used as a title (optional)
%
% [figI]        figure number (optional), generates autonumbered next
%               figure if this argumant is missing or empty.
%
% [right]       a string specifying the orientation of the
%               slices. if right is 'right', then right is 
%               right, i.e. the view onto the slices is from
%               above. if right is 'left' or simply 'wrong',
%               then right is left, i.e. the view onto the 
%               slices is from below.
%
% [skipNslices] optional number of slices to be skipped between slices to
%               be shown. defaults to 0 (all slices shown).
%
% labels        structured array of string labels (optional)
%
% labelPos      nLabels by 3 matrix of label positions (optional)

% DEBUG
% vol(1,1,:,1)=[1 0 0];
% vol(10,1,:,1)=[1 0 0];

if ndims(volORmap)==3 % if its a map...
    vol=map2vol(volORmap);
else % must be a vol
    vol=volORmap;
end

if ~exist('skipNslices','var')
    sliceIs=size(vol,4):-1:1;
else
    sliceIs=1:skipNslices+1:size(vol,4);
end

vol=vol(:,:,:,sliceIs);

if nargin<7 || nHorPanels==0
    nHorPanels=ceil(sqrt(size(vol,4)));
    nVerPanels=ceil(sqrt(size(vol,4)));
end

if ~exist('right','var')
    right='right';
end

if ~exist('figI','var') || (exist('figI','var') && isempty(figI))
    figI=0;
end

if ~exist('title','var')
    title='';
end

if figI(1)
    figure(figI(1));
    if numel(figI)>1
        subplot(figI(2),figI(3),figI(4)); cla;
        if numel(figI)>4
            % vol=vol(:,:,:,figI(5));
            vol=vol(5:75,40:75,:,figI(5));
        end
    else
        clf;
    end
else
    figI=figure;
end

set(figI(1),'Color','w');
%set(figI,'WindowStyle','docked');

margin=0;
width=(1-(nHorPanels-1)*margin)/nHorPanels;
height=(1-(nVerPanels-1)*margin)/nVerPanels;

[sx,sy,sc,sz]=size(vol);

for sliceI=1:size(vol,4)
    %subplot(nVerPanels,nHorPanels,sliceI);
    left=mod(sliceI-1,nHorPanels)*(width+margin);
    bottom=1-ceil(sliceI/nHorPanels)*(height+margin);
    
    if numel(figI)<2
        subplot('Position',[left bottom width height]);
    end
    
    if strcmp(right,'right')
        im=flipdim(permute(vol(:,:,:,sliceI),[2 1 3]),2);
    else
        im=permute(vol(:,:,:,sliceI),[2 1 3]);
    end
    image(im);
    
    axis equal;
    axis([0 size(im,2)-1 0 size(im,1)-1])

    set(gca,'Visible','off');
end

if exist('labels','var') && strcmp(right,'right')
    labelPos(:,1)=sx-labelPos(:,1); % flip x
end

h=axes('Parent',gcf); hold on;
set(h,'Visible','off');
axis([0 1 0 1]);

text(0.5,-0.04,['Right is ',right,'.'],'HorizontalAlignment','Center','Color',[0.7 0.7 0.7],'FontSize',7);
% if strcmp(right,'right')
%     text(0.5,-0.04,'Right is right.','HorizontalAlignment','Center','Color',[0.7 0.7 0.7],'FontSize',7);
% else
%     text(0.5,-0.04,'Right is left.','HorizontalAlignment','Center','Color',[0.7 0.7 0.7],'FontSize',7);
% end


% text(0.5,-0.08,title,'HorizontalAlignment','Center','FontSize',7,'FontWeight','bold','Color',[.4 .4 .4]);
text(0.5,0.05,title,'HorizontalAlignment','Center','FontSize',14,'FontWeight','bold','Color',[.8 .8 .8]);
%text(0.5,1.05,title,'HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','Color','k');

if exist('labels','var')
    for labelI=1:length(labels);
        sliceI=labelPos(labelI,3);
        left=mod(sliceI-1,nHorPanels)*(width+margin);
        bottom=1-ceil(sliceI/nHorPanels)*(height+margin);
        subplot('Position',[left bottom width height]);
        
        % text(labelPos(labelI,1),labelPos(labelI,2),labels(labelI),'BackgroundColor','w','HorizontalAlignment','center');
        text(labelPos(labelI,1),labelPos(labelI,2),labels(labelI),'Color','w','FontWeight','bold','HorizontalAlignment','center');
        text(labelPos(labelI,1),labelPos(labelI,2),labels(labelI),'Color','k','FontWeight','normal','HorizontalAlignment','center');
    end
end
return

function [threshold,binMap]=autothresholdMap(map,nVox,twosided)

% returns the threshold that highlights nVox voxels in the statistical map
% map. if twosided is passed and 1, then the returned threshold highlights
% nVox voxels when applied to the absolute values of map.

sz=size(map);
map=map(:); % convert to a column vector (explicit conversion only needed because map can be a row vector)

if ~exist('twosided','var'), twosided=0; end;

if isnan(sum(map(:)))
    disp('autothresholdMap: NaNs were found and set to zero.');
    map(isnan(map))=0;
end

if twosided
    map=abs(map);
end

% to speed this up...
% eligibilityFactor=0.1; 
%eligibilityFactor=0.05; % should always work
%eligibleVoxelI=find(map>eligibilityFactor*std(map(:)));
% if size(eligibleVoxelI,1)<=nVox
%     disp('nEligibleVoxels:'); size(eligibleVoxelI,1)
%     disp('nVox:'); nVox
%     disp('eligibilityFactor:'); eligibilityFactor
%     error('ERROR: number of eligible voxels needs to be larger than nVox. please decrease the eligibilityFactor.');
% end

%eligibleVoxelI=find(map); % consider all voxels with a stat value greater than zero
eligibleVoxelI=(1:numel(map))'; % consider all voxels
if nVox>numel(eligibleVoxelI)
    nVox=numel(eligibleVoxelI);
end

eligibleVoxelI_statVal=[eligibleVoxelI, map(eligibleVoxelI)];
clear eligibleVoxelI;

eligibleVoxelI_statVal_sorted=flipud(sortrows(eligibleVoxelI_statVal,2));
eligibleVoxelI_statVal_sorted=[eligibleVoxelI_statVal_sorted;[nVox+1,-10e10]]; % add one below the lowest in case the threshold is supposed to mark all voxels
clear eligibleVoxelI_statVal;
%disp('sorted')

% select nVox many voxels from the top of the list
%topVoxelI=eligibleVoxelI_statVal_sorted(1:nVox,1);

if nVox>0
    threshold=(eligibleVoxelI_statVal_sorted(nVox,2)+eligibleVoxelI_statVal_sorted(nVox+1,2))/2;
else % nVox==0
    threshold=max(map(:))+1;
end
    
clear eligibleVoxelI_statVal_sorted;

map=reshape(map,sz);
binMap=map>threshold;
return;
