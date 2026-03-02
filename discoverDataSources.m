function sources = discoverDataSources(campaignRoot, excludePatterns, maxDepth, progressFcn, cancelFcn)
% discoverDataSources Lightweight auto-discovery of known sensor roots.
% Uses conservative folder-name heuristics to avoid flooding the UI with
% leaf folders (e.g., "hsic" output folders).
    if nargin < 2 || isempty(excludePatterns)
        excludePatterns = {'@tmp','@eaDir','.DS_Store','thumbs.db'};
    end
    if nargin < 3 || isempty(maxDepth)
        maxDepth = 3;
    end
    if nargin < 4 || isempty(progressFcn)
        progressFcn = @(~)[];
    end
    if nargin < 5 || isempty(cancelFcn)
        cancelFcn = @()false;
    end

    sources = struct('id', {}, 'type', {}, 'label', {}, 'rootPath', {}, 'enabled', {}, 'config', {});
    if isempty(campaignRoot) || ~isfolder(campaignRoot)
        return;
    end

    q = {struct('path', campaignRoot, 'depth', 0)};
    seenPath = containers.Map('KeyType','char','ValueType','logical');
    seenTypeAtRoot = containers.Map('KeyType','char','ValueType','logical');
    dirsVisited = 0;

    while ~isempty(q)
        if cancelFcn(), return; end
        cur = q{1}; q(1) = [];
        dirsVisited = dirsVisited + 1;
        if mod(dirsVisited, 100) == 0
            progressFcn(struct('dirsVisited', dirsVisited, 'currentPath', cur.path));
        end

        listing = dir(cur.path);
        for i = 1:numel(listing)
            name = listing(i).name;
            if ~listing(i).isdir || strcmp(name,'.') || strcmp(name,'..')
                continue;
            end
            if shouldExclude(name, excludePatterns)
                continue;
            end
            child = fullfile(cur.path, name);

            [stype, label] = classifySource(name);
            if ~strcmp(stype, 'UNKNOWN') && ~isRedundantLeaf(name)
                typeRootKey = [stype '|' nearestTopLevel(campaignRoot, child)];
                if ~isKey(seenPath, child) && ~isKey(seenTypeAtRoot, typeRootKey)
                    seenPath(child) = true;
                    seenTypeAtRoot(typeRootKey) = true;
                    sid = lower(regexprep([stype '_' label], '[^a-zA-Z0-9_]', '_'));
                    sources(end+1) = struct('id', sid, 'type', stype, 'label', label, ...
                        'rootPath', child, 'enabled', true, 'config', struct()); %#ok<AGROW>
                end
            end

            if cur.depth < maxDepth
                q{end+1} = struct('path', child, 'depth', cur.depth+1); %#ok<AGROW>
            end
        end
    end

    progressFcn(struct('dirsVisited', dirsVisited, 'currentPath', campaignRoot));
end

function tf = shouldExclude(name, patterns)
    low = lower(name);
    tf = false;
    for p = 1:numel(patterns)
        pat = lower(string(patterns{p}));
        if contains(low, pat)
            tf = true;
            return;
        end
    end
end

function tf = isRedundantLeaf(name)
    lname = lower(string(name));
    % Common leaf folder names that should never appear as top-level source rows.
    tf = any(strcmp(lname, ["hsic","hdr","raw","cal","cubes","frames","output","outputs"]));
end

function key = nearestTopLevel(rootPath, childPath)
    rel = strrep(childPath, [rootPath filesep], '');
    parts = strsplit(rel, filesep);
    if isempty(parts)
        key = childPath;
    else
        key = parts{1};
    end
end

function [stype, label] = classifySource(name)
    lname = lower(string(name));

    if contains(lname, "fridge")
        stype = 'FRIDGE'; label = char(name); return;
    elseif strcmp(lname, "hsi")
        stype = 'HSI'; label = char(name); return;
    elseif contains(lname, "lwir")
        stype = 'LWIR'; label = char(name); return;
    elseif contains(lname, "mwir")
        stype = 'MWIR'; label = char(name); return;
    elseif startsWith(lname, "fast")
        stype = 'FAST'; label = char(name); return;
    elseif contains(lname, "cerberus") || contains(lname, "mx20") || contains(lname, "mx-20")
        stype = 'CERBERUS'; label = char(name); return;
    end

    stype = 'UNKNOWN';
    label = char(name);
end
