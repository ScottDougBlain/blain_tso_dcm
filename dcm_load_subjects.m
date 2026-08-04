function subjects = dcm_load_subjects(csv_file, varargin)
%DCM_LOAD_SUBJECTS Load subject information from CSV file
%
% Reads subject information from a CSV file and returns a structured
% array with subject IDs, runs, group labels, and exclusion flags.
%
% USAGE:
%   subjects = dcm_load_subjects('subjects.csv')
%   subjects = dcm_load_subjects('subjects.csv', 'exclude', {'1023', '2020'})
%   subjects = dcm_load_subjects('subjects.csv', 'groups', {'HC'})
%
% INPUTS:
%   csv_file - Path to CSV file with subject information
%
% OPTIONAL KEY-VALUE PAIRS:
%   'exclude'      - Cell array of subject IDs to exclude
%   'include_only' - Cell array of subject IDs to include (overrides exclude)
%   'groups'       - Cell array of group labels to include
%   'exclude_flagged' - Logical, exclude subjects with Exclude=1 in CSV (default: true)
%
% CSV FORMAT:
%   Required columns: SubjectID, Group
%   Optional columns: Run1, Run2, ..., RunN, Exclude, Notes
%
%   Example:
%   SubjectID,Run1,Run2,Run3,Run4,Run5,Run6,Group,Exclude,Notes
%   1001,1,1,1,1,1,1,HC,0,
%   1002,1,1,1,1,1,1,HC,0,
%   1003,1,1,0,1,1,1,HC,0,Missing run 3
%   2001,1,1,1,1,1,1,SZ,0,
%   2020,1,1,1,1,1,1,SZ,1,Excessive motion
%
% OUTPUTS:
%   subjects - Structure array with fields:
%       .id       - Subject ID (string)
%       .runs     - Vector of run numbers to include
%       .group    - Group label (string)
%       .excluded - Logical, whether subject is excluded
%       .notes    - Any notes from CSV
%
% EXAMPLES:
%   % Load all non-excluded subjects
%   subj = dcm_load_subjects('subjects.csv');
%
%   % Exclude additional subjects
%   subj = dcm_load_subjects('subjects.csv', 'exclude', {'1023', '2020'});
%
%   % Load only healthy controls
%   subj = dcm_load_subjects('subjects.csv', 'groups', {'HC'});
%
% SEE ALSO: dcm_run_pipeline, dcm_config_template

% Parse optional inputs
p = inputParser;
addRequired(p, 'csv_file', @ischar);
addParameter(p, 'exclude', {}, @iscell);
addParameter(p, 'include_only', {}, @iscell);
addParameter(p, 'groups', {}, @iscell);
addParameter(p, 'exclude_flagged', true, @islogical);
parse(p, csv_file, varargin{:});

opts = p.Results;

% Check file exists
if ~exist(csv_file, 'file')
    error('Subject CSV file not found: %s', csv_file);
end

% Read CSV file
try
    T = readtable(csv_file, 'TextType', 'string', 'VariableNamingRule', 'preserve');
catch err
    error('Failed to read CSV file: %s\n%s', csv_file, err.message);
end

% Validate required columns
required_cols = {'SubjectID', 'Group'};
for i = 1:length(required_cols)
    if ~ismember(required_cols{i}, T.Properties.VariableNames)
        error('CSV missing required column: %s', required_cols{i});
    end
end

% Find run columns (Run1, Run2, etc.)
run_cols = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'Run'));
n_runs = length(run_cols);

% Check for optional columns
has_exclude = ismember('Exclude', T.Properties.VariableNames);
has_notes = ismember('Notes', T.Properties.VariableNames);

% Initialize output structure
n_subj = height(T);
subjects = struct('id', {}, 'runs', {}, 'group', {}, 'excluded', {}, 'notes', {});

% Process each subject
for i = 1:n_subj
    subj = struct();

    % Get subject ID (ensure it's a string)
    % Handle both numeric and string/categorical SubjectID columns
    sid = T.SubjectID(i);
    if isnumeric(sid)
        subj.id = num2str(sid);
    elseif iscell(sid)
        subj.id = char(sid{1});
    else
        subj.id = char(sid);
    end
    
    % Get runs to include
    if n_runs > 0
        runs = [];
        for r = 1:n_runs
            col_name = run_cols{r};
            if T.(col_name)(i) == 1
                runs = [runs, r]; %#ok<AGROW>
            end
        end
        subj.runs = runs;
    else
        % If no run columns, assume all 6 runs
        subj.runs = 1:6;
    end
    
    % Get group
    subj.group = char(T.Group(i));
    
    % Get notes
    if has_notes && ~ismissing(T.Notes(i))
        subj.notes = char(T.Notes(i));
    else
        subj.notes = '';
    end
    
    % Determine exclusion status
    subj.excluded = false;
    
    % Check CSV exclude flag
    if opts.exclude_flagged && has_exclude
        if T.Exclude(i) == 1
            subj.excluded = true;
        end
    end
    
    % Check manual exclude list
    if ~isempty(opts.exclude) && ismember(subj.id, opts.exclude)
        subj.excluded = true;
    end
    
    % Check include_only list (overrides other exclusions)
    if ~isempty(opts.include_only)
        if ~ismember(subj.id, opts.include_only)
            subj.excluded = true;
        else
            subj.excluded = false;  % Include even if flagged
        end
    end
    
    % Check group filter
    if ~isempty(opts.groups)
        if ~ismember(subj.group, opts.groups)
            subj.excluded = true;
        end
    end
    
    subjects(i) = subj;
end

% Summary
n_included = sum(~[subjects.excluded]);
n_excluded = sum([subjects.excluded]);
fprintf('Loaded %d subjects: %d included, %d excluded\n', n_subj, n_included, n_excluded);

% List excluded subjects if any
if n_excluded > 0
    excluded_ids = {subjects([subjects.excluded]).id};
    fprintf('Excluded: %s\n', strjoin(excluded_ids, ', '));
end

end
