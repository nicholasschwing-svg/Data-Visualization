function test_index_store()
    tmpDb = fullfile(tempdir, ['timeline_idx_' char(java.util.UUID.randomUUID) '.sqlite']);
    c = onCleanup(@() cleanup(tmpDb)); %#ok<NASGU>

    IndexStore('init', tmpDb);
    src = struct('id','fridge_test','type','FRIDGE','label','FRIDGE','rootPath',tempdir,'enabled',true,'config',struct());
    IndexStore('replacesources', tmpDb, src);

    item = struct('item_id','abc','source_id','fridge_test','file_path','/tmp/a.hdr','timestamp_utc',int64(1000), ...
        'duration_ms',[], 'kind','cube','metadata_json','{}','file_mtime',int64(1),'file_size',int64(2));
    IndexStore('upsertitem', tmpDb, item);

    rows = IndexStore('queryrange', tmpDb, {'fridge_test'}, int64(0), int64(5000));
    assert(height(rows)==1, 'Expected one indexed row.');

    [dt, modality] = parseFastFilenameTime('2024-11-20_1817_17_point-00_LWIR_cal_hsi.hdr'); %#ok<ASGLU>
    assert(~isnat(dt), 'FAST parse should produce a datetime.');
    assert(strcmp(modality, 'LWIR'), 'Expected LWIR modality token.');
end

function cleanup(p)
    if isfile(p), delete(p); end
end
