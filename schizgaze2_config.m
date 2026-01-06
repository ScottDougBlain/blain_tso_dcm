%% schizgaze2_config.m
% DCM Pipeline Configuration for SchizGaze2 Study
% ================================================
% Gaze perception DCM analysis comparing HC and SZ groups
%
% Authors: Scott Blain, Ivy Tso
% Created: November 2025

%% Initialize configuration structure
config = struct();

%% ========================================================================
%  STUDY INFORMATION
%  ========================================================================
config.study_name = 'SchizGaze2';
config.study_description = 'DCM analysis of gaze perception in schizophrenia';
config.analyst = 'Scott Blain';
config.date_created = '2025-11-15';

%% ========================================================================
%  PATH CONFIGURATION
%  ========================================================================
% Base paths - MODIFY THESE FOR YOUR SYSTEM
config.paths.project_root = '/nfs/net/fatty/Data/SchizGaze2_16';
config.paths.spm_path = '/nii/bionic/MethodsCore12/SPM/SPM12/spm12_with_R7219';

% Output directories
config.paths.output_base = '[project_root]/ScottB_SchizGaze/DCM_Pipeline_Output';
config.paths.log_dir = '[output_base]/Logs';

% Template paths
config.paths.templates.images = '[project_root]/Subjects/[Subject]/func/gaze/run_0[Run]/';
config.paths.templates.firstlevel = '[project_root]/ScottB_SchizGaze/DCMConcat2/[Subject]/';
config.paths.templates.dcm_output = '[output_base]/DCM/[Subject]/';
config.paths.templates.voi_output = '[output_base]/VOI/[Subject]/';

% Image file filter
config.paths.image_filter = 's6w3ut.*nii';

%% ========================================================================
%  SUBJECT DEFINITIONS
%  ========================================================================
% Load subjects from CSV file
config.subjects_file = '[project_root]/ScottB_SchizGaze/DCM_Pipeline_Output/schizgaze2_subjects.csv';

% Alternative: Copy the CSV to the config directory and use:
% config.subjects_file = fullfile(fileparts(mfilename('fullpath')), 'schizgaze2_subjects.csv');

%% ========================================================================
%  ACQUISITION PARAMETERS
%  ========================================================================
config.acquisition.TR = 0.8;              % Repetition time (seconds)
config.acquisition.TE = 0.03;             % Echo time (seconds)
config.acquisition.n_slices = 16;         % Number of slices
config.acquisition.ref_slice = 8;         % Reference slice for timing

%% ========================================================================
%  ROI DEFINITIONS
%  ========================================================================
% ROIs for gaze perception network
% Format: {Name, Center [x y z], Radius, SubjSpec_Center, SubjSpec_Radius, SubjSpec_Mode}

config.rois = {
    % Primary visual cortex
    'V1',           [-6 -88 0],         5,      [], [], [];
    
    % Fusiform face area
    'FFA',          [42 -52 -18],       5,      [], [], [];
    
    % Posterior superior temporal sulcus (gaze processing)
    'pSTS',         [54 -42 6],         5,      [], [], [];
    
    % Cerebellum (optional - included in some analyses)
    'Cerebellum',   [-26 -66 -26],      5,      [], [], [];
};

config.n_rois = size(config.rois, 1);

%% ========================================================================
%  EXPERIMENTAL CONDITIONS
%  ========================================================================
% Conditions from the gaze task (must match first-level GLM)
config.conditions = {
    'Face';     % Face presentation
    'Gaze';     % Gaze shift
};

config.n_conditions = size(config.conditions, 1);

%% ========================================================================
%  VOI EXTRACTION OPTIONS
%  ========================================================================
config.voi.contrast_num = 1;              % Effects of interest contrast
config.voi.threshold_p = 1;               % Unthresholded extraction
config.voi.threshold_extent = 0;          % No cluster threshold
config.voi.correction = 'none';           % No correction
config.voi.eoi_contrast = 0;              % F-contrast for adjustment (0 = skip)

%% ========================================================================
%  DCM MODEL SPECIFICATION
%  ========================================================================
nROI = config.n_rois;
nCond = config.n_conditions;

% --- A Matrix: Intrinsic (fixed) connectivity ---
% Full model with all possible connections (to be pruned by PEB/BMR)
% Diagonal = self-connections (required)
% A(target, source) = connection from source to target

config.dcm.A = zeros(nROI, nROI);

% Self-connections (always on)
config.dcm.A = config.dcm.A + 2*eye(nROI);

% Feedforward pathway: V1 -> FFA -> pSTS
config.dcm.A(2, 1) = 2;   % V1 -> FFA
config.dcm.A(3, 2) = 2;   % FFA -> pSTS

% Feedback connections
config.dcm.A(1, 2) = 1;   % FFA -> V1 (optional)
config.dcm.A(2, 3) = 1;   % pSTS -> FFA (optional)

% Cerebellar connections (if using 4 ROIs)
if nROI >= 4
    config.dcm.A(4, 3) = 1;   % pSTS -> Cerebellum
    config.dcm.A(3, 4) = 1;   % Cerebellum -> pSTS
    config.dcm.A(4, 2) = 1;   % FFA -> Cerebellum
end

% --- B Matrices: Modulatory effects ---
config.dcm.B = cell(nCond, 1);

% Condition 1: Face - no modulation (baseline)
config.dcm.B{1} = zeros(nROI, nROI);

% Condition 2: Gaze - modulates feedforward connections
config.dcm.B{2} = zeros(nROI, nROI);
config.dcm.B{2}(2, 1) = 1;   % Gaze modulates V1 -> FFA
config.dcm.B{2}(3, 2) = 1;   % Gaze modulates FFA -> pSTS
if nROI >= 4
    config.dcm.B{2}(4, 3) = 1;   % Gaze modulates pSTS -> Cerebellum
end

% --- C Matrix: Driving inputs ---
% Visual stimuli drive V1
config.dcm.C = zeros(nROI, nCond);
config.dcm.C(1, :) = 2;   % All conditions drive V1

%% ========================================================================
%  DCM OPTIONS
%  ========================================================================
config.dcm.options.nonlinear = 0;         % Linear DCM
config.dcm.options.two_state = 0;         % Single-state model
config.dcm.options.stochastic = 0;        % Deterministic
config.dcm.options.centre = 1;            % Centre inputs
config.dcm.options.induced = 0;           % Not spectral DCM

%% ========================================================================
%  PEB OPTIONS
%  ========================================================================
config.peb.fields = {'A', 'B'};           % Analyze A and B matrices
config.peb.Q = 'all';                     % Random effects on all parameters

%% ========================================================================
%  PIPELINE CONTROL
%  ========================================================================
% Set to 0 to skip a stage
config.steps.extract_voi = 1;             % Stage 1: VOI extraction
config.steps.specify_dcm = 1;             % Stage 2: DCM specification
config.steps.estimate_dcm = 1;            % Stage 3: DCM estimation
config.steps.run_peb = 1;                 % Stage 4: PEB analysis
config.steps.export_results = 1;          % Stage 5: Export to CSV

% Error handling
config.error_handling = 'continue';       % 'continue' or 'stop'

%% ========================================================================
%  OUTPUT OPTIONS
%  ========================================================================
config.output.save_workspace = true;
config.output.save_figures = true;
config.output.export_csv = true;
config.output.verbose = 2;

%% ========================================================================
%  USAGE EXAMPLES
%  ========================================================================
% To run this configuration:
%
% 1. Full pipeline:
%    results = dcm_run_pipeline('config/schizgaze2_config.m');
%
% 2. Exclude specific subjects:
%    results = dcm_run_pipeline('config/schizgaze2_config.m', ...
%                               'exclude', {'1023', '2020'});
%
% 3. Run only healthy controls:
%    results = dcm_run_pipeline('config/schizgaze2_config.m', ...
%                               'groups', {'HC'});
%
% 4. Exclude from PEB only:
%    results = dcm_run_pipeline('config/schizgaze2_config.m', ...
%                               'exclude_peb', {'1050'});
