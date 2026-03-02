function stats = updateWorkspaceIndex(workspace, progressFcn, cancelFcn)
% updateWorkspaceIndex Incremental per-directory indexing with cancellation.
    if nargin < 2 || isempty(progressFcn), progressFcn = @(~)[]; end
    if nargin < 3 || isempty(cancelFcn), cancelFcn = @()false; end
    IndexStore('init', workspace.indexDbPath);
    IndexStore('replacesources', workspace.indexDbPath, workspace.sources);

    stats = struct('dirsScanned',0,'filesScanned',0,'warnings',0);
    for si = 1:numel(workspace.sources)
        src = workspace.sources(si);
        if ~src.enabled || ~isfolder(src.rootPath)
            continue;
        end
        stack = {src.rootPath};
        while ~isempty(stack)
            if cancelFcn(), return; end
            dirPath = stack{end}; stack(end)=[];
            [isExc,~] = excludedPath(dirPath, workspace.excludePatterns);
            if isExc, continue; end

            info = dir(dirPath);
            if isempty(info), continue; end
            dirMtime = int64(info(1).datenum * 86400);
            prev = IndexStore('dirmtime', workspace.indexDbPath, dirPath);
            if ~isempty(prev) && int64(prev) == dirMtime
                continue;
            end

            listing = dir(dirPath);
            stats.dirsScanned = stats.dirsScanned + 1;
            progressFcn(struct('source',src.label,'currentPath',dirPath,'dirsScanned',stats.dirsScanned,'filesScanned',stats.filesScanned));
            for i = 1:numel(listing)
                nm = listing(i).name;
                if strcmp(nm,'.') || strcmp(nm,'..'), continue; end
                p = fullfile(dirPath, nm);
                if listing(i).isdir
                    stack{end+1} = p; %#ok<AGROW>
                    continue;
                end
                [isExcFile,~] = excludedPath(p, workspace.excludePatterns);
                if isExcFile, continue; end
                if ~matchesSourceFile(src.type, nm), continue; end
                [tsMs, kind, metadata] = parseTimestampForFile(src.type, p);
                item = struct();
                item.item_id = makeId(src.id, p);
                item.source_id = src.id;
                item.file_path = p;
                item.timestamp_utc = tsMs;
                item.duration_ms = [];
                item.kind = kind;
                item.metadata_json = jsonencode(metadata);
                item.file_mtime = int64(listing(i).datenum * 86400);
                item.file_size = int64(listing(i).bytes);
                IndexStore('upsertitem', workspace.indexDbPath, item);
                stats.filesScanned = stats.filesScanned + 1;
                if isempty(tsMs), stats.warnings = stats.warnings + 1; end
                if mod(stats.filesScanned, 200) == 0
                    progressFcn(struct('source',src.label,'currentPath',p,'dirsScanned',stats.dirsScanned,'filesScanned',stats.filesScanned));
                    if cancelFcn(), return; end
                end
            end
            IndexStore('touchdir', workspace.indexDbPath, src.id, dirPath, dirMtime);
        end
    end
end

function [tf, pat] = excludedPath(pathStr, pats)
    low = lower(pathStr); tf = false; pat = '';
    for i = 1:numel(pats)
        p = lower(string(pats{i}));
        if contains(low, p)
            tf = true; pat = p; return;
        end
    end
end

function tf = matchesSourceFile(stype, name)
    [~,~,ext] = fileparts(lower(name));
    switch upper(stype)
        case 'FRIDGE', tf = strcmp(ext,'.hdr') && contains(lower(name),'aaro');
        case {'HSI','CERBERUS','LWIR','MWIR','FAST'}, tf = strcmp(ext,'.hdr') || strcmp(ext,'.hsic') || strcmp(ext,'.raw');
        otherwise, tf = true;
    end
end

function [tsMs, kind, metadata] = parseTimestampForFile(stype, filePath)
    [~,fname,ext] = fileparts(filePath);
    metadata = struct('parse', 'ok');
    tsMs = [];
    kind = 'frame';
    switch upper(stype)
        case {'CERBERUS','LWIR','MWIR','HSI'}
            dt = parseCerbFilenameTime([fname ext], '(?<year>\d{4})-(?<month>\d{2})-(?<day>\d{2})_(?<hour>\d{2})-(?<min>\d{2})-(?<sec>\d{2})');
            if ~isnat(dt), tsMs = int64(posixtime(dt)*1000); end
        case 'FAST'
            [dt, modality, pointId] = parseFastFilenameTime([fname ext]);
            metadata.modality = modality; metadata.pointId = pointId;
            if ~isnat(dt), tsMs = int64(posixtime(dt)*1000); end
        case 'FRIDGE'
            [tStart, tEnd, waveLabel, captureKey] = parseFridgeHeader(filePath);
            metadata.wavelength = waveLabel; metadata.captureKey = captureKey;
            if ~isnat(tEnd)
                metadata.endMs = int64(posixtime(tEnd)*1000);
            end
            if ~isnat(tStart), tsMs = int64(posixtime(tStart)*1000); kind = 'cube'; end
        otherwise
            metadata.parse = 'unknown';
    end
end

function id = makeId(sourceId, filePath)
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8([sourceId '|' filePath]));
    id = lower(reshape(dec2hex(typecast(md.digest(),'uint8'))',1,[]));
end
