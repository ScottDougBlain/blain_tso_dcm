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
config.paths.project_root = '/fs/ess/PAS2478/';
config.paths.spm_path = '/fs/ess/PAS2478/Lab_Resources/matlab_toolboxes/spm12';

% Output directories (will be created if they don't exist)
config.paths.output_base = '[project_root]/Projects/DCM_Analysis_Test';
config.paths.log_dir = '[output_base]/Logs';

% Template paths - use [Variable] for substitution
% Available variables: project_root, output_base, data_base, Subject, Run
%
% IMPORTANT: firstlevel is where VOI extraction looks for SPM.mat files.
% This path is automatically updated by the pipeline if:
%   - concat.enabled=true: Updated to concat_output after concatenation
%   - glm.run_dcm=true: Updated to glm.output.dcm after GLM creation
config.paths.templates.firstlevel = '[project_root]/FirstLevel/[Subject]/';  % Where SPM.mat files are
config.paths.templates.dcm_output = '[output_base]/DCM/[Subject]/';
config.paths.templates.voi_output = '[output_base]/VOI/[Subject]/';

%% ========================================================================
%  DATA LOCATIONS
%  ========================================================================
% Define where your functional images, timing files, and motion parameters
% live. These are defined ONCE here and automatically propagated to the
% timing concatenation utility (config.timing) and the GLM pipeline
% (config.glm) via the auto-populate block at the bottom of this file.
%
% You only need to edit this section — do NOT duplicate these values in
% the timing or GLM sections below.

%% --- Functional Images ---
% Template path to functional image directories (one per run per subject)
% Uses [Variable] substitution: [project_root], [Subject], [Run]
config.data.images_dir = '[project_root]/Datasets/SchizGaze2/MRI/Subjects/[Subject]/func/gaze/run_[Run]/';
config.data.images_filter = '^s6w3.*\.nii$';  % Regex filter for spm_select

%% --- Timing CSV ---
% Master timing file containing trial onsets, conditions, etc.
config.data.timing_dir = '[project_root]/Datasets/SchizGaze2/MRI/MRI_Timing';
config.data.timing_file = 'SchizGaze2_masterfile_fMRI_event.csv';

% Column names in the timing CSV
config.data.timing_columns.subject = 'Subject';       % Subject ID column
config.data.timing_columns.run = 'RunList';            % Run number column
config.data.timing_columns.onset = 'TrialOnset';       % Trial onset time column
config.data.timing_columns.duration = 'Duration';      % Trial duration column
config.data.timing_columns.concat_onset = 'Onset_Concat';  % Concatenated onset column (if pre-computed)

%% --- Motion Parameters ---
% Template path to motion parameter files (same directory structure as images)
config.data.motion_dir = '[project_root]/Datasets/SchizGaze2/MRI/Subjects/[Subject]/func/gaze/run_[Run]/';
config.data.motion_patterns = {'rp_utrun_*.txt', 'rp_*.txt', 'rp*.txt'};  % Filename patterns to try (in order)
config.data.motion_derivatives = true;  % Include temporal derivatives (6 params + 6 derivatives = 12)

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
    '1009',       [1 2 3 4 5 6],  'HC',  0;
    '1014',       [1 2 3 4 5 6],  'HC',  0;
    '2002',       [1 2 3 4 5 6],  'SZ',  0;
    '2003',       [1 2 3 4 5 6],  'SZ',  0;
};

%% ========================================================================
%  ACQUISITION PARAMETERS
%  ========================================================================
config.acquisition.TR = 0.8;              % Repetition time (seconds)
config.acquisition.TE = 0.03;             % Echo time (seconds)
config.acquisition.n_slices = 60;         % Number of slices
config.acquisition.ref_slice = 10;        % Reference slice for timing
config.acquisition.volumes_per_run = 275; % Volumes per run (for concatenation)

%% ========================================================================
%  RUN CONCATENATION (for DCM)
%  ========================================================================
% DCM requires a single-session model. If your first-level has multiple
% runs/sessions, enable concatenation to merge them before VOI extraction.
%
% This stage:
%   1. Copies SPM.mat to concat output directory
%   2. Runs spm_fmri_concatenate to merge sessions
%   3. Re-estimates the concatenated model
%   4. Adds an F-contrast for Effects of Interest (EOI)
%
config.concat.enabled = false;            % Set to true to enable concatenation
config.concat.overwrite = false;          % Overwrite existing concatenated models

% Output path for concatenated models (use [Variable] substitution)
config.paths.templates.concat_output = '[output_base]/Concat/[Subject]/';

%% ========================================================================
%  TIMING FILE CONCATENATION (for DCM)
%  ========================================================================
% Optional preprocessing step: reads the timing CSV (defined in DATA
% LOCATIONS above) and creates a new version with a ConcatOnset column
% where onsets accumulate across runs rather than resetting to 0.
%
% Formula: ConcatOnset = TrialOnset + (volumes_per_run * TR) * (Run - 1)
%
% Example with volumes_per_run=228 and TR=2 (456s per run):
%   Run 1, Trial at 2s    -> ConcatOnset = 2
%   Run 2, Trial at 2s    -> ConcatOnset = 458  (2 + 456)
%   Run 3, Trial at 2s    -> ConcatOnset = 914  (2 + 912)
%
config.timing.enabled = false;            % Set to true to create concatenated timing
config.timing.overwrite = false;          % Overwrite existing output file

% Output: Path for concatenated timing CSV (will add ConcatOnset column)
config.timing.output_file = '[output_base]/Timing/masterfile_DCM.csv';

% Input file and column names are auto-populated from DATA LOCATIONS above.
% Override here ONLY if the timing concat utility uses a different file or
% different column names than the GLM (uncommon):
%   config.timing.master_file = '/path/to/different/file.csv';
%   config.timing.onset_column = 'DifferentOnsetColumn';
%   config.timing.run_column = 'DifferentRunColumn';

%% ========================================================================
%  FIRST-LEVEL GLM CONFIGURATION
%  ========================================================================
% CONTEXT: DCM requires a SINGLE-SESSION first-level GLM. If your fMRI data
% has multiple runs, you must merge them before VOI extraction. There are two
% ways to do this:
%
%   (A) config.concat (above) - Merge EXISTING multi-session GLMs
%   (B) config.glm.run_dcm   - Build a NEW single-session GLM from scratch
%
% Enable ONE of these, not both. Whichever runs last redirects the firstlevel
% path that all downstream pipeline stages use.
%
% Additionally, a standard GLM can be run for ROI identification:
%
%   (C) config.glm.run_standard - Standard multi-session GLM (independent,
%       does NOT feed into DCM, used only for activation maps/ROI finding)
%
% These are controlled by independent flags and can be run separately.
%
% NOTE: Image paths, timing file, and motion parameters are defined in the
% DATA LOCATIONS section above and auto-populated into config.glm.* at the
% bottom of this file. You do NOT need to set them here.

%% --- Enable/Disable GLM Stages ---
config.glm.run_standard = 0;          % Set to 1 to run standard GLM
config.glm.run_dcm = 0;               % Set to 1 to run DCM GLM
config.glm.overwrite = false;         % Overwrite existing GLM output

%% --- Standard GLM Conditions ---
% Condition mapping: {Name, Column, Values, Filter}
%   - Name: SPM condition name
%   - Column: CSV column to filter on
%   - Values: Value(s) that define this condition (vector for OR logic)
%   - Filter: Optional additional filter expression (empty for none)
config.glm.standard.conditions = {
    'Gaze',     'TaskType',  1,   '';
    'Gender',   'TaskType',  2,   '';
};

%% --- DCM GLM Conditions ---
config.glm.dcm.conditions = {
    'AllTrials', 'TaskType',  [1 2],  '';    % All face trials (Gaze + Gender)
    'Gaze',      'TaskType',  1,      '';    % Gaze trials only
};

%% --- Concatenated Onset Handling (DCM GLM only) ---
% When building the DCM GLM, onsets must be concatenated across runs.
% Two options: use a pre-existing column in the CSV, or compute on the fly.
config.glm.timing.compute_concat_onsets = false;  % true=compute from onset+run, false=use existing column

%% --- Parametric Modulators (optional) ---
config.glm.pmod.enabled = false;              % Toggle parametric modulators
config.glm.pmod.definitions = {
    % Condition   PmodName      Column        Polynomial order
    'Gaze',      'GazeAngle',  'GazeSignal',  1;
};

%% --- Model Settings ---
config.glm.model.hrf = 'canonical';           % HRF type: 'canonical', 'canonical+time', 'canonical+time+disp'
config.glm.model.high_pass_filter = 128;      % HPF cutoff in seconds

%% --- Contrast Settings ---
config.glm.contrasts.auto_generate = true;    % Auto-generate contrasts
config.glm.contrasts.include_pairwise = true; % Include pairwise comparisons (standard GLM)
config.glm.contrasts.eoi_contrast = true;     % Include EOI F-contrast (DCM GLM)

%% --- Output Paths ---
config.glm.output.standard = '[output_base]/GLM/Standard/[Subject]/';
config.glm.output.dcm = '[output_base]/GLM/DCM/[Subject]/';

%% ========================================================================
%  ROI DEFINITIONS
%  ========================================================================
% Define regions of interest for VOI extraction.
%
% TWO EXTRACTION MODES ARE SUPPORTED:
%
% 1. FIXED SPHERE MODE (3 columns):
%    Format: {Name, Center, Radius}
%    - Extracts eigenvariate from all voxels in a fixed sphere
%    - Center: MNI coordinates [x y z]
%    - Radius: Sphere radius in mm
%
% 2. SUBJECT-SPECIFIC PEAK MODE (6 columns):
%    Format: {Name, Search_Center, Search_Radius, Extract_Center, Extract_Radius, Mode}
%
%    Workflow:
%      1. Create SEARCH sphere (larger, e.g., 10mm) at group coordinates
%      2. Find subject-specific activation peak within search sphere
%      3. Move EXTRACTION sphere (smaller, e.g., 5mm) to that peak
%      4. Extract eigenvariate from moved extraction sphere
%
%    Parameters:
%      - Name:           ROI name (used in output filenames)
%      - Search_Center:  MNI coordinates [x y z] for search sphere
%      - Search_Radius:  Search sphere radius in mm (constrains peak search)
%      - Extract_Center: MNI coordinates [x y z] for extraction sphere
%                        (typically same as Search_Center)
%      - Extract_Radius: Extraction sphere radius in mm (smaller than search)
%      - Mode:           'local' = nearest peak to center
%                        'global' = highest peak in search sphere
%
% NOTES:
%   - Minimum 2 ROIs required, maximum 8 recommended
%   - For subject-specific mode, threshold settings (config.voi.threshold_*)
%     control what counts as a "peak" (e.g., p < 0.001 uncorrected)
%   - If no suprathreshold voxels exist in search sphere, extraction fails
%
% EXAMPLE - Subject-specific peak finding (recommended for group studies):
config.rois = {
    % Name         Search_Center      Search_R  Extract_Center     Extract_R  Mode
    'Cerebellum', [-8 -76 -29],       10,       [-8 -76 -29],      5,         'local';
    'Insula',     [35 20 -3],         10,       [35 20 -3],        5,         'local';
    'IPL',        [44 -45 45],        10,       [44 -45 45],       5,         'local';
    'Fusiform',   [40 -64 -10],       10,       [40 -64 -10],      5,         'local';
};

% EXAMPLE - Fixed sphere extraction (simpler, use when activation is consistent):
% config.rois = {
%     % Name         Center           Radius
%     'Cerebellum', [-8 -76 -29],     5;
%     'Insula',     [35 20 -3],       5;
%     'IPL',        [44 -45 45],      5;
%     'Fusiform',   [40 -64 -10],     5;
% };

config.n_rois = size(config.rois, 1);

%% ========================================================================
%  EXPERIMENTAL CONDITIONS
%  ========================================================================
% List conditions to include in DCM (must match names in first-level GLM)
config.conditions = {
    'Face';     % Condition 1: Face (baseline) - drives regions, not modulated
    'Gaze';     % Condition 2: Gaze - modulates connections
};

config.n_conditions = size(config.conditions, 1);

%% ========================================================================
%  VOI EXTRACTION OPTIONS
%  ========================================================================
% These settings control how VOI time-series are extracted from the GLM.
%
% CONTRAST THRESHOLDING:
%   For subject-specific peak finding, threshold settings control which
%   voxels are considered when searching for peaks. Common approaches:
%
%   - Unthresholded (p=1): Include all voxels in search sphere
%     Pros: Always finds a peak; Cons: May pick noise
%
%   - Liberal threshold (p=0.05, uncorrected): Include moderately active voxels
%     Good balance for most analyses
%
%   - Standard threshold (p=0.001, uncorrected): Conservative peak detection
%     Recommended for subject-specific extraction
%
%   - Corrected threshold (p=0.05, FWE): Very conservative
%     May fail if no voxels survive correction in search sphere
%
% EXTENT THRESHOLD:
%   Minimum cluster size (in voxels) for peak detection. Set to 0 for
%   subject-specific mode (peaks are already spatially constrained by
%   search sphere). Use higher values (e.g., 10) for whole-brain searches.

config.voi.contrast_num = 1;              % Contrast index for thresholding/peak finding
config.voi.threshold_p = 0.001;           % p-value threshold (1 = unthresholded)
config.voi.threshold_extent = 0;          % Minimum cluster size (voxels)
config.voi.correction = 'none';           % 'none' or 'FWE'

% EOI ADJUSTMENT:
%   The Effects of Interest (EOI) F-contrast removes variance from nuisance
%   regressors (motion, session constants) while preserving task-related
%   variance. This is HIGHLY RECOMMENDED for DCM, especially with
%   concatenated runs.
%
%   - eoi_contrast = 0: No adjustment (not recommended)
%   - eoi_contrast = N: Use existing F-contrast at index N
%   - eoi_contrast = -1: Auto-find existing EOI contrast
%   - add_eoi_contrast = true: Create EOI contrast if it doesn't exist

config.voi.eoi_contrast = -1;             % F-contrast for EOI adjustment (-1 = auto-find)
config.voi.add_eoi_contrast = true;       % Auto-add EOI F-contrast if missing

%% ========================================================================
%  DCM MODEL SPECIFICATION
%  ========================================================================
% DCM models effective connectivity using three matrices:
%
%   A = Intrinsic (fixed) connectivity - always active, baseline connections
%   B = Modulatory connectivity - how experimental conditions change connections
%   C = Driving inputs - which regions receive external task input
%
% MATRIX ENCODING:
%   Use binary values: 0 or 1
%     0 = exclude from model (connection not estimated)
%     1 = include in model (connection strength will be estimated)
%
%   During estimation, SPM estimates the actual connection strengths
%   (can be positive or negative). PEB/BMR can later prune connections
%   that don't contribute at the group level.
%
% INTERPRETATION:
%   - Positive A values = excitatory baseline connectivity
%   - Negative A values = inhibitory baseline connectivity
%   - Positive B values = condition increases connectivity strength
%   - Negative B values = condition decreases connectivity strength
%   - C values = how strongly task input drives each region
%
% MATRIX ORIENTATION:
%   A(i,j) = connection FROM region j TO region i
%   B(i,j,k) = how condition k modulates connection FROM j TO i
%   C(i,k) = input FROM condition k TO region i
%
%   Example: A(2,1) = 1 means connection from region 1 to region 2
%
% ROI ORDER: Defined by config.rois (here: 1=Cerebellum, 2=Insula, 3=IPL, 4=Fusiform)

nROI = config.n_rois;     % Number of ROIs (should be 4 in this example)
nCond = config.n_conditions;  % Number of conditions (should be 2)

% -------------------------------------------------------------------------
% A MATRIX: Intrinsic (fixed/endogenous) connectivity
% -------------------------------------------------------------------------
% These connections are ALWAYS active, regardless of task condition.
% Think of them as the "baseline wiring" of the network.
%
% A(target, source) = connection from source to target
%
% Common approaches:
%   - Full model: All connections = 1 (let PEB determine which matter)
%   - Hypothesis-driven: Only include theoretically plausible connections
%   - Hierarchical: Enforce feedforward/feedback structure
%
%        Source regions (columns)
%        Cereb  Ins   IPL   Fus
% Target   |     |     |     |
% Cereb   [1     1     1     1]   <- connections TO Cerebellum
% Ins     [1     1     1     1]   <- connections TO Insula
% IPL     [1     1     1     1]   <- connections TO IPL
% Fus     [1     1     1     1]   <- connections TO Fusiform

config.dcm.A = ones(nROI, nROI);  % Full connectivity
% This is a FULL model - all regions connected bidirectionally.
% Self-connections (diagonal) represent intrinsic inhibition within each region.
% PEB analysis will determine which connections are significant.

% -------------------------------------------------------------------------
% B MATRICES: Modulatory (bilinear) effects
% -------------------------------------------------------------------------
% These specify HOW experimental conditions CHANGE connectivity.
% B values are ADDED to A values when the condition is active.
%
% Example: If A(2,1)=0.5 and B(2,1,2)=0.3 for condition 2,
%          then during condition 2, effective connectivity = 0.5 + 0.3 = 0.8
%
% One B matrix per condition, same dimensions as A (nROI x nROI).
%
% Common approaches:
%   - Baseline condition: B = zeros (no modulation)
%   - Task condition: B specifies which connections are modulated
%   - Self-connections: Usually NOT modulated (diagonal = 0)

config.dcm.B = cell(nCond, 1);  % Cell array with one matrix per condition

% Condition 1: Face (baseline) - no modulation
config.dcm.B{1} = zeros(nROI, nROI);

% Condition 2: GAZE - modulates all between-region connections
% Diagonal = 0 because we don't modulate self-connections
config.dcm.B{2} = [
%   Cereb Ins  IPL  Fus
    0     1    1    1;    % TO Cerebellum: modulated by Ins, IPL, Fus
    1     0    1    1;    % TO Insula: modulated by Cereb, IPL, Fus
    1     1    0    1;    % TO IPL: modulated by Cereb, Ins, Fus
    1     1    1    0;    % TO Fusiform: modulated by Cereb, Ins, IPL
];
% Interpretation: Gaze condition can change all between-region connectivity.
% PEB will determine which modulations are significant.

% -------------------------------------------------------------------------
% C MATRIX: Driving inputs
% -------------------------------------------------------------------------
% Specifies which regions receive DIRECT input from experimental conditions.
% This is the "entry point" for task-related information into the network.
%
% C(region, condition) = input to region from condition
%
% Common approaches:
%   - Single entry point: One region receives all input, propagates to others
%   - Multiple entry points: Several regions receive direct input
%   - Sensory input: Often goes to primary sensory regions first
%
% NOTE: At least one region must receive driving input for the model to work.

config.dcm.C = zeros(nROI, nCond);  % Initialize as no inputs

% Face condition drives ALL 4 regions
config.dcm.C(:, 1) = 1;
% Gaze doesn't drive regions directly - it only modulates connections (via B)
%
% Alternative: Single entry point (e.g., only Cerebellum receives input)
% config.dcm.C(1, 1) = 1;  % Cerebellum receives Face input only

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
config.steps.specify_dcm = 0;             % Stage 2: DCM specification
config.steps.estimate_dcm = 0;            % Stage 3: DCM estimation
config.steps.run_peb = 0;                 % Stage 4: PEB analysis
config.steps.export_results = 0;          % Stage 5: Export to CSV

% Error handling: 'continue' to proceed despite errors, 'stop' to halt
config.error_handling = 'continue';

%% ========================================================================
%  PATH REMAPPING (for migrated data)
%  ========================================================================
% If SPM.mat files contain old paths (from data migration), specify
% the old and new base paths here. The pipeline will update SPM.mat
% files before VOI extraction.
config.path_remap.enabled = true;
config.path_remap.old_base = '/net/data4/SchizGaze2_16';
config.path_remap.new_base = '/fs/ess/PAS2478/Datasets/SchizGaze2/MRI';

% Additional string replacements for filename mismatches
% Each row: {old_pattern, new_pattern}
config.path_remap.replacements = {
    's6w3urun_', 's6w3utrun_';   % Fix filename pattern mismatch
};

%% ========================================================================
%  OUTPUT OPTIONS
%  ========================================================================
config.output.save_workspace = true;      % Save full MATLAB workspace
config.output.save_figures = true;        % Save diagnostic figures
config.output.export_csv = true;          % Export parameters to CSV
config.output.verbose = 2;                % 0=quiet, 1=normal, 2=verbose

%% ========================================================================
%  AUTO-POPULATE (do not edit below this line)
%  ========================================================================
% Propagates DATA LOCATIONS into the config.timing and config.glm fields
% that the pipeline scripts expect. This avoids duplicating paths and
% column names across sections.
%
% If you need the timing concatenation utility or the GLM to use DIFFERENT
% file paths or column names than what's in DATA LOCATIONS, override the
% specific field AFTER this block.

% --- Timing concatenation utility (dcm_create_concat_timing.m) ---
if ~isfield(config.timing, 'master_file')
    config.timing.master_file = fullfile(config.data.timing_dir, config.data.timing_file);
end
if ~isfield(config.timing, 'onset_column')
    config.timing.onset_column = config.data.timing_columns.onset;
end
if ~isfield(config.timing, 'run_column')
    config.timing.run_column = config.data.timing_columns.run;
end

% --- GLM pipeline (glm_load_timing.m, glm_build_batch.m, etc.) ---
config.glm.images.template = config.data.images_dir;
config.glm.images.filter = config.data.images_filter;

config.glm.timing.dir = config.data.timing_dir;
config.glm.timing.file = config.data.timing_file;
config.glm.timing.subject_column = config.data.timing_columns.subject;
config.glm.timing.run_column = config.data.timing_columns.run;
config.glm.timing.onset_column = config.data.timing_columns.onset;
config.glm.timing.duration_column = config.data.timing_columns.duration;
config.glm.timing.concat_onset_column = config.data.timing_columns.concat_onset;

config.glm.motion.dir_template = config.data.motion_dir;
config.glm.motion.patterns = config.data.motion_patterns;
config.glm.motion.include_derivatives = config.data.motion_derivatives;
