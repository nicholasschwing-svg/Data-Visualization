function varargout = WorkspaceManager(action, varargin)
% WorkspaceManager Centralized workspace persistence helpers.
% Actions:
%   default         -> struct with default workspace values
%   normalize, ws   -> fill missing defaults and normalize derived fields
%   save, ws, path  -> persist workspace JSON
%   load, path      -> load workspace JSON
%   setlast, path   -> save last-workspace pointer
%   getlast         -> load last-workspace pointer

    switch lower(action)
        case 'default'
            varargout{1} = defaultWorkspace();
        case 'normalize'
            ws = varargin{1};
            varargout{1} = normalizeWorkspace(ws);
        case 'save'
            ws = varargin{1};
            outPath = varargin{2};
            txt = jsonencode(ws, 'PrettyPrint', true);
            fid = fopen(outPath, 'w');
            if fid < 0
                error('WorkspaceManager:WriteFailed', 'Unable to write workspace file: %s', outPath);
            end
            c = onCleanup(@() fclose(fid)); %#ok<NASGU>
            fwrite(fid, txt, 'char');
            varargout{1} = outPath;
        case 'load'
            inPath = varargin{1};
            raw = fileread(inPath);
            ws = jsondecode(raw);
            varargout{1} = normalizeWorkspace(ws);
        case 'setlast'
            wsPath = varargin{1};
            ptr = struct('lastWorkspacePath', wsPath, 'updatedAt', posixtime(datetime('now')));
            ptrPath = lastPointerPath();
            if ~isfolder(fileparts(ptrPath)), mkdir(fileparts(ptrPath)); end
            fid = fopen(ptrPath, 'w');
            c = onCleanup(@() fclose(fid)); %#ok<NASGU>
            fwrite(fid, jsonencode(ptr, 'PrettyPrint', true), 'char');
            varargout{1} = ptrPath;
        case 'getlast'
            ptrPath = lastPointerPath();
            if ~isfile(ptrPath)
                varargout{1} = '';
                return;
            end
            ptr = jsondecode(fileread(ptrPath));
            if isfield(ptr, 'lastWorkspacePath')
                varargout{1} = char(ptr.lastWorkspacePath);
            else
                varargout{1} = '';
            end
        otherwise
            error('WorkspaceManager:BadAction', 'Unknown action: %s', action);
    end
end

function ws = defaultWorkspace()
    cacheDir = fullfile(prefdir, 'timeline_app_cache');
    if ~isfolder(cacheDir), mkdir(cacheDir); end
    defaultDbPath = deriveIndexDbPath('', 'default-workspace');
    ws = struct( ...
        'name', 'default-workspace', ...
        'campaignRoot', '', ...
        'sources', struct('id', {}, 'type', {}, 'label', {}, 'rootPath', {}, 'enabled', {}, 'config', {}), ...
        'excludePatterns', {{'@tmp', '@eaDir', '.DS_Store', 'thumbs.db'}}, ...
        'includePatterns', {{}}, ...
        'indexDbPath', defaultDbPath, ...
        'lastOpenedAt', posixtime(datetime('now')));
end

function ws = normalizeWorkspace(ws)
    d = defaultWorkspace();
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(ws, f{i})
            ws.(f{i}) = d.(f{i});
        end
    end
    if isempty(ws.sources)
        ws.sources = d.sources;
    end
    ws.indexDbPath = normalizeIndexDbPath(ws, d.indexDbPath);
end

function p = lastPointerPath()
    p = fullfile(prefdir, 'timeline_app_cache', 'last_workspace.json');
end

function outPath = normalizeIndexDbPath(ws, fallbackPath)
    outPath = char(string(ws.indexDbPath));
    legacyPath = fullfile(prefdir, 'timeline_app_cache', 'timeline_index.sqlite');
    shouldDerive = isempty(strtrim(outPath)) || strcmpi(outPath, legacyPath);
    if shouldDerive
        wsName = '';
        if isfield(ws, 'name')
            wsName = ws.name;
        end
        rootPath = '';
        if isfield(ws, 'campaignRoot')
            rootPath = ws.campaignRoot;
        end
        outPath = deriveIndexDbPath(rootPath, wsName);
    elseif nargin >= 2 && isempty(outPath)
        outPath = fallbackPath;
    end
end

function dbPath = deriveIndexDbPath(campaignRoot, wsName)
    cacheDir = fullfile(prefdir, 'timeline_app_cache');
    if ~isfolder(cacheDir), mkdir(cacheDir); end
    fingerprint = lower(strtrim(char(string(campaignRoot))));
    if isempty(fingerprint)
        fingerprint = lower(strtrim(char(string(wsName))));
    end
    if isempty(fingerprint)
        fingerprint = 'default';
    end
    fingerprint = regexprep(fingerprint, '[^a-z0-9]+', '_');
    fingerprint = regexprep(fingerprint, '^_+|_+$', '');
    if isempty(fingerprint)
        fingerprint = 'workspace';
    end
    if strlength(fingerprint) > 80
        fingerprint = extractBefore(fingerprint, 81);
    end
    dbPath = fullfile(cacheDir, sprintf('timeline_index_%s.sqlite', fingerprint));
end
