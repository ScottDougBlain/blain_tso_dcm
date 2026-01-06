%% dcm_config_template.m
% DCM Pipeline Configuration Template
% ====================================
% Copy this file to config/[study_name]_config.m and modify for your study
%
% This configuration file defines all parameters needed to run a complete
% DCM analysis from VOI extraction through group-level PEB inference.
%
% Usage:
%   1. Copy this template: cp dcm_config_template.m config/mystudy_config.m
%   2. Edit the copy with your study parameters
%   3. Run: results = dcm_run_pipeline('config/mystudy_config.m')

%% Initialize configuration structure
config = struct();

%% ========================================================================
%  STUDY INFORMATION
%  ========================================================================
config.study_name = 'StudyName';
config.study_description = 'Brief description of DCM analysis';
config.analyst = 'Your Name';
config.date_created = datestr(now, 'yyyy-mm-dd');

%% ========================================================================
%  PATH CONFIGURATION
%  ========================================================================
% Base paths - modify these for your system/study
config.paths.project_root = '/path/to/your/project';
config.paths.spm_path = '/path/to/SPM12';

% Output directories (will be created if they don't exist)
config.paths.output_base = '[project_root]/DCM_Analysis';
config.paths.log_dir = '[output_base]/Logs';

% Template paths - use [Variable] for substitution
% Available variables: project_root, output_base, Subject, Run, OutputName
config.paths.templates.images = '[project_root]/Subjects/[Subject]/func/task/run_0[Run]/';
config.paths.templates.firstlevel = '[project_root]/FirstLevel/[Subject]/';
config.paths.templates.dcm_output = '[output_base]/DCM/[Subject]/';
config.paths.templates.voi_output = '[output_base]/VOI/[Subject]/';

% Image file filter (regex for spm_select)
config.paths.image_filter = 's6w3ut.*nii';

%% ========================================================================
%  SUBJECT DEFINITIONS
%  ========================================================================
% Option 1: Load from CSV file (recommended)
% config.subjects_file = '[project_root]/subjects.csv';
%
% CSV format: SubjectID,Run1,Run2,...,RunN,Group,Exclude,Notes
% Example:
%   SubjectID,Run1,Run2,Run3,Run4,Run5,Run6,Group,Exclude,Notes
%   1001,1,1,1,1,1,1,HC,0,
%   1002,1,1,1,1,1,1,HC,0,
%   2001,1,1,1,1,1,1,SZ,0,
%   2002,1,1,1,1,1,1,SZ,1,Excessive motion

% Option 2: Define subjects inline (for smaller studies)
config.subjects = {
    % SubjectID,  Runs,           Group, Exclude
    '1001',       [1 2 3 4 5 6],  'HC',  0;
    '1002',       [1 2 3 4 5 6],  'HC',  0;
    '2001',       [1 2 3 4 5 6],  'SZ',  0;
    '2002',       [1 2 3 4 5 6],  'SZ',  0;
};

%% ========================================================================
%  ACQUISITION PARAMETERS
%  ========================================================================
config.acquisition.TR = 2.0;              % Repetition time (seconds)
config.acquisition.TE = 0.03;             % Echo time (seconds)
config.acquisition.n_slices = 32;         % Number of slices
config.acquisition.ref_slice = 16;        % Reference slice for timing

%% ========================================================================
%  ROI DEFINITIONS
%  ========================================================================
% Define regions of interest for VOI extraction
% Format: {Name, Center [x y z], Radius, SubjSpec_Center, SubjSpec_Radius, SubjSpec_Mode}
%
% For standard spherical ROIs, leave last 3 columns as [],[],[]
% Minimum 2 ROIs required, maximum 8 recommended

config.rois = {
    % Name          Center [x y z]      Radius  SubjSpec params
    'ROI1',         [0 0 0],            5,      [], [], [];
    'ROI2',         [10 10 10],         5,      [], [], [];
    'ROI3',         [20 20 20],         5,      [], [], [];
};

config.n_rois = size(config.rois, 1);

%% ========================================================================
%  EXPERIMENTAL CONDITIONS
%  ========================================================================
% List conditions to include in DCM (must match names in first-level GLM)
config.conditions = {
    'Condition1';
    'Condition2';
};

config.n_conditions = size(config.conditions, 1);

%% ========================================================================
%  VOI EXTRACTION OPTIONS
%  ========================================================================
config.voi.contrast_num = 1;              % Contrast index for thresholding
config.voi.threshold_p = 1;               % p-value threshold (1 = unthresholded)
config.voi.threshold_extent = 0;          % Minimum cluster size
config.voi.correction = 'none';           % 'none' or 'FWE'
config.voi.eoi_contrast = 0;              % F-contrast for EOI adjustment (0 = skip)

%% ========================================================================
%  DCM MODEL SPECIFICATION
%  ========================================================================
% Define the model using A, B, C matrices
% Values: 0 = exclude, 1 = optional (for future BMR), 2 = required

nROI = config.n_rois;
nCond = config.n_conditions;

% --- A Matrix: Intrinsic (fixed) connectivity ---
% A(target, source) = connection from source to target
% Diagonal = self-connections
config.dcm.A = zeros(nROI, nROI);
config.dcm.A = config.dcm.A + 2*eye(nROI);  % Self-connections required
% Add your connections here, e.g.:
% config.dcm.A(2, 1) = 2;   % ROI1 -> ROI2

% --- B Matrices: Modulatory effects ---
% B{condition}(target, source) = modulation of connection by condition
config.dcm.B = cell(nCond, 1);
for c = 1:nCond
    config.dcm.B{c} = zeros(nROI, nROI);
end
% Add modulations here, e.g.:
% config.dcm.B{2}(2, 1) = 1;   % Condition2 modulates ROI1->ROI2

% --- C Matrix: Driving inputs ---
% C(region, condition) = input to region from condition
config.dcm.C = zeros(nROI, nCond);
config.dcm.C(1, :) = 2;   % First ROI receives all inputs (typical for visual)

%% ========================================================================
%  DCM OPTIONS
%  ========================================================================
config.dcm.options.nonlinear = 0;         % Nonlinear DCM (0 = linear)
config.dcm.options.two_state = 0;         % Two-state neuronal model
config.dcm.options.stochastic = 0;        % Stochastic DCM
config.dcm.options.centre = 1;            % Centre inputs

%% ========================================================================
%  PEB OPTIONS
%  ========================================================================
config.peb.fields = {'A', 'B'};           % Parameters to analyze at group level
config.peb.Q = 'all';                     % Random effects specification

%% ========================================================================
%  PIPELINE CONTROL
%  ========================================================================
% Set to 0 to skip a stage
config.steps.extract_voi = 1;             % Stage 1: VOI extraction
config.steps.specify_dcm = 1;             % Stage 2: DCM specification
config.steps.estimate_dcm = 1;            % Stage 3: DCM estimation
config.steps.run_peb = 1;                 % Stage 4: PEB analysis
config.steps.export_results = 1;          % Stage 5: Export to CSV

% Error handling: 'continue' to proceed despite errors, 'stop' to halt
config.error_handling = 'continue';

%% ========================================================================
%  OUTPUT OPTIONS
%  ========================================================================
config.output.save_workspace = true;      % Save full MATLAB workspace
config.output.save_figures = true;        % Save diagnostic figures
config.output.export_csv = true;          % Export parameters to CSV
config.output.verbose = 2;                % 0=quiet, 1=normal, 2=verbose
