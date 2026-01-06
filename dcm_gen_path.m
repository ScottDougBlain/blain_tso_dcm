function path_out = dcm_gen_path(template, varargin)
%DCM_GEN_PATH Generate path from template with variable substitution
%
% Takes a template string containing [Variable] placeholders and replaces
% them with actual values from the configuration or provided key-value pairs.
%
% USAGE:
%   path = dcm_gen_path(template, config)
%   path = dcm_gen_path(template, 'Key1', Value1, 'Key2', Value2, ...)
%   path = dcm_gen_path(template, config, 'Key1', Value1, ...)
%
% INPUTS:
%   template  - Path template string with [Variable] placeholders
%               Example: '[project_root]/Subjects/[Subject]/func/'
%
%   config    - Configuration structure (optional)
%               Recognized fields from config.paths:
%               - project_root, output_base, spm_path, etc.
%
%   Key/Value - Additional variable substitutions (override config)
%               Common keys: Subject, Run, OutputName, iRun, iSubject
%
% OUTPUTS:
%   path_out  - Resolved path string
%
% SUPPORTED VARIABLES:
%   [project_root]  - Base project directory
%   [output_base]   - Output directory base
%   [spm_path]      - SPM installation path
%   [Subject]       - Current subject ID
%   [Run]           - Current run number (string)
%   [iRun]          - Current run index (number)
%   [iSubject]      - Current subject index (number)
%   [OutputName]    - Output directory name
%   [date]          - Current date (YYYYMMDD)
%
% EXAMPLES:
%   % Using config structure
%   config.paths.project_root = '/data/study1';
%   path = dcm_gen_path('[project_root]/Subjects/[Subject]/', config, 'Subject', '1001');
%   % Returns: '/data/study1/Subjects/1001/'
%
%   % Using key-value pairs only
%   path = dcm_gen_path('/data/[Subject]/run[Run]/', 'Subject', '1001', 'Run', '2');
%   % Returns: '/data/1001/run2/'
%
% SEE ALSO: dcm_run_pipeline, dcm_config_template

% Initialize variables map
vars = struct();

% Parse inputs
config = [];
extra_vars = {};

if nargin >= 2
    if isstruct(varargin{1})
        config = varargin{1};
        extra_vars = varargin(2:end);
    else
        extra_vars = varargin;
    end
end

% Extract variables from config.paths
if ~isempty(config) && isfield(config, 'paths')
    path_fields = fieldnames(config.paths);
    for i = 1:length(path_fields)
        field = path_fields{i};
        val = config.paths.(field);
        if ischar(val)
            vars.(field) = val;
        end
    end
end

% Add common derived variables
vars.date = datestr(now, 'yyyymmdd');

% Parse extra key-value pairs
for i = 1:2:length(extra_vars)
    if i+1 <= length(extra_vars)
        key = extra_vars{i};
        val = extra_vars{i+1};
        if isnumeric(val)
            val = num2str(val);
        end
        vars.(key) = val;
    end
end

% Perform substitutions (iterate to handle nested variables)
path_out = template;
max_iterations = 10;  % Prevent infinite loops
for iter = 1:max_iterations
    path_prev = path_out;
    
    % Find all [Variable] patterns
    pattern = '\[([^\]]+)\]';
    matches = regexp(path_out, pattern, 'tokens');
    
    for m = 1:length(matches)
        var_name = matches{m}{1};
        
        % Check if variable exists
        if isfield(vars, var_name)
            % Replace [Variable] with its value
            path_out = strrep(path_out, ['[' var_name ']'], vars.(var_name));
        end
    end
    
    % Stop if no more substitutions were made
    if strcmp(path_out, path_prev)
        break;
    end
end

% Warn about unresolved variables
unresolved = regexp(path_out, '\[([^\]]+)\]', 'tokens');
if ~isempty(unresolved)
    unresolved_names = cellfun(@(x) x{1}, unresolved, 'UniformOutput', false);
    warning('DCM:UnresolvedPath', 'Unresolved variables in path: %s', strjoin(unresolved_names, ', '));
end

% Clean up path (remove double slashes, except in protocol like http://)
path_out = regexprep(path_out, '([^:])//+', '$1/');

end
