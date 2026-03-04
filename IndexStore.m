function varargout = IndexStore(action, varargin)
% IndexStore SQLite-backed item index for timeline queries.
    switch lower(action)
        case 'init'
            dbPath = varargin{1};
            if ~isfolder(fileparts(dbPath)), mkdir(fileparts(dbPath)); end
            if isfile(dbPath)
                conn = sqlite(dbPath, 'connect');
            else
                conn = sqlite(dbPath, 'create');
            end
            cleanupObj = onCleanup(@() close(conn)); %#ok<NASGU>
            exec(conn, ['CREATE TABLE IF NOT EXISTS sources(' ...
                'source_id TEXT PRIMARY KEY, type TEXT, label TEXT, root_path TEXT, enabled INTEGER, created_at INTEGER, updated_at INTEGER)']);
            exec(conn, ['CREATE TABLE IF NOT EXISTS directories(' ...
                'dir_path TEXT PRIMARY KEY, source_id TEXT, mtime INTEGER, last_scanned_at INTEGER)']);
            exec(conn, ['CREATE TABLE IF NOT EXISTS items(' ...
                'item_id TEXT PRIMARY KEY, source_id TEXT, file_path TEXT, timestamp_utc INTEGER, duration_ms INTEGER, kind TEXT, metadata_json TEXT, file_mtime INTEGER, file_size INTEGER)']);
            exec(conn, 'CREATE INDEX IF NOT EXISTS idx_items_source_time ON items(source_id, timestamp_utc)');
            exec(conn, 'CREATE INDEX IF NOT EXISTS idx_items_file_path ON items(file_path)');
            exec(conn, 'CREATE INDEX IF NOT EXISTS idx_dirs_source ON directories(source_id)');
        case 'replacesources'
            dbPath = varargin{1}; sources = varargin{2};
            conn = sqlite(dbPath, 'connect'); c = onCleanup(@() close(conn)); %#ok<NASGU>
            nowTs = int64(posixtime(datetime('now'))*1000);
            for i = 1:numel(sources)
                s = sources(i);
                sql = sprintf(['INSERT INTO sources(source_id,type,label,root_path,enabled,created_at,updated_at) VALUES(%s,%s,%s,%s,%d,%d,%d) ' ...
                    'ON CONFLICT(source_id) DO UPDATE SET type=excluded.type,label=excluded.label,root_path=excluded.root_path,enabled=excluded.enabled,updated_at=excluded.updated_at'], ...
                    q(s.id), q(s.type), q(s.label), q(s.rootPath), logical(s.enabled), nowTs, nowTs);
                exec(conn, sql);
            end
        case 'dirmtime'
            dbPath = varargin{1}; dirPath = varargin{2};
            conn = sqlite(dbPath, 'connect'); c = onCleanup(@() close(conn)); %#ok<NASGU>
            rows = fetch(conn, sprintf('SELECT mtime FROM directories WHERE dir_path=%s', q(dirPath)));
            varargout{1} = firstScalar(rows);
        case 'touchdir'
            dbPath = varargin{1}; sourceId = varargin{2}; dirPath = varargin{3}; mtime = varargin{4};
            conn = sqlite(dbPath, 'connect'); c = onCleanup(@() close(conn)); %#ok<NASGU>
            nowTs = int64(posixtime(datetime('now'))*1000);
            exec(conn, sprintf(['INSERT INTO directories(dir_path,source_id,mtime,last_scanned_at) VALUES(%s,%s,%d,%d) ' ...
                'ON CONFLICT(dir_path) DO UPDATE SET source_id=excluded.source_id,mtime=excluded.mtime,last_scanned_at=excluded.last_scanned_at'], ...
                q(dirPath), q(sourceId), mtime, nowTs));
        case 'upsertitem'
            dbPath = varargin{1}; item = varargin{2};
            conn = sqlite(dbPath, 'connect'); c = onCleanup(@() close(conn)); %#ok<NASGU>
            exec(conn, sprintf(['INSERT INTO items(item_id,source_id,file_path,timestamp_utc,duration_ms,kind,metadata_json,file_mtime,file_size) ' ...
                'VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s) ON CONFLICT(item_id) DO UPDATE SET timestamp_utc=excluded.timestamp_utc,duration_ms=excluded.duration_ms,kind=excluded.kind,metadata_json=excluded.metadata_json,file_mtime=excluded.file_mtime,file_size=excluded.file_size'], ...
                q(item.item_id), q(item.source_id), q(item.file_path), n(item.timestamp_utc), n(item.duration_ms), q(item.kind), q(item.metadata_json), n(item.file_mtime), n(item.file_size)));
        case 'queryrange'
            dbPath = varargin{1}; sourceIds = varargin{2}; startMs = varargin{3}; endMs = varargin{4};
            if isempty(sourceIds), varargout{1} = table(); return; end
            in = strjoin(cellfun(@q, sourceIds, 'UniformOutput', false), ',');
            conn = sqlite(dbPath, 'connect'); c = onCleanup(@() close(conn)); %#ok<NASGU>
            rows = fetch(conn, sprintf(['SELECT item_id,source_id,file_path,timestamp_utc,' ...
                'COALESCE(duration_ms,0) AS duration_ms,' ...
                'COALESCE(kind,'''') AS kind,' ...
                'COALESCE(metadata_json,''{}'') AS metadata_json,' ...
                'COALESCE(file_mtime,0) AS file_mtime,' ...
                'COALESCE(file_size,0) AS file_size FROM items ' ...
                'WHERE source_id IN (%s) AND timestamp_utc BETWEEN %d AND %d ORDER BY timestamp_utc'], in, int64(startMs), int64(endMs)));
            varargout{1} = toTable(rows, {'item_id','source_id','file_path','timestamp_utc','duration_ms','kind','metadata_json','file_mtime','file_size'});
        case 'summaries'
            dbPath = varargin{1};
            conn = sqlite(dbPath, 'connect'); c = onCleanup(@() close(conn)); %#ok<NASGU>
            rows = fetch(conn, ['SELECT s.source_id,s.label,s.type,s.root_path,s.enabled,COUNT(i.item_id),MIN(i.timestamp_utc),MAX(i.timestamp_utc) ' ...
                'FROM sources s LEFT JOIN items i ON s.source_id=i.source_id GROUP BY s.source_id,s.label,s.type,s.root_path,s.enabled']);
            varargout{1} = toTable(rows, {'source_id','label','type','root_path','enabled','item_count','min_ts','max_ts'});
        otherwise
            error('IndexStore:UnknownAction', 'Unknown action: %s', action);
    end
end

function out = q(s)
    out = ['''' strrep(char(string(s)), '''', '''''') ''''];
end
function out = n(v)
    if isempty(v) || (isnumeric(v) && isnan(v))
        out = 'NULL';
    else
        out = sprintf('%d', int64(v));
    end
end


function v = firstScalar(rows)
    if isempty(rows)
        v = [];
        return;
    end
    if istable(rows)
        if height(rows) == 0
            v = [];
            return;
        end
        v = rows{1,1};
        return;
    end
    if iscell(rows)
        v = rows{1};
    else
        v = rows(1);
    end
end

function t = toTable(rows, names)
    if isempty(rows)
        t = table();
        return;
    end
    if istable(rows)
        t = rows;
        if width(t) == numel(names)
            t.Properties.VariableNames = names;
        end
        return;
    end
    t = cell2table(rows, 'VariableNames', names);
end
