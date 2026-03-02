function sources = discoverDataSources(campaignRoot, excludePatterns, maxDepth)
% discoverDataSources Lightweight auto-discovery of known sensor roots.
    if nargin < 2 || isempty(excludePatterns)
        excludePatterns = {'@tmp','@eaDir','.DS_Store','thumbs.db'};
    end
    if nargin < 3 || isempty(maxDepth)
        maxDepth = 3;
    end
    sources = struct('id', {}, 'type', {}, 'label', {}, 'rootPath', {}, 'enabled', {}, 'config', {});
    if isempty(campaignRoot) || ~isfolder(campaignRoot)
        return;
    end

    q = {struct('path', campaignRoot, 'depth', 0)};
    seen = containers.Map('KeyType','char','ValueType','logical');

    while ~isempty(q)
        cur = q{1}; q(1) = [];
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
            if ~strcmp(stype, 'UNKNOWN') && ~isKey(seen, child)
                seen(child) = true;
                sid = lower(regexprep([stype '_' label '_' char(java.util.UUID.randomUUID)], '[^a-zA-Z0-9_]', '_'));
                sources(end+1) = struct('id', sid, 'type', stype, 'label', label, ...
                    'rootPath', child, 'enabled', true, 'config', struct()); %#ok<AGROW>
            end
            if cur.depth < maxDepth
                q{end+1} = struct('path', child, 'depth', cur.depth+1); %#ok<AGROW>
            end
        end
    end
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

function [stype, label] = classifySource(name)
    lname = lower(name);
    if contains(lname, 'fridge')
        stype = 'FRIDGE'; label = name; return;
    elseif contains(lname, 'hsi')
        stype = 'HSI'; label = name; return;
    elseif contains(lname, 'lwir')
        stype = 'LWIR'; label = name; return;
    elseif contains(lname, 'mwir')
        stype = 'MWIR'; label = name; return;
    elseif contains(lname, 'fast')
        stype = 'FAST'; label = name; return;
    elseif contains(lname, 'cerberus') || contains(lname, 'mx20')
        stype = 'CERBERUS'; label = name; return;
    end
    stype = 'UNKNOWN';
    label = name;
end
