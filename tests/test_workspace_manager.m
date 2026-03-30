function test_workspace_manager()
    ws = WorkspaceManager('default');
    assert(contains(ws.indexDbPath, 'timeline_index_'), 'Expected derived index DB name.');
    assert(~strcmp(ws.indexDbPath, fullfile(prefdir, 'timeline_app_cache', 'timeline_index.sqlite')), ...
        'Default workspace should no longer use the shared legacy DB name.');

    legacy = ws;
    legacy.campaignRoot = fullfile(tempdir, 'Drive_A', 'Campaign_1');
    legacy.indexDbPath = fullfile(prefdir, 'timeline_app_cache', 'timeline_index.sqlite');
    normA = WorkspaceManager('normalize', legacy);
    assert(~strcmp(normA.indexDbPath, legacy.indexDbPath), 'Legacy DB path should be upgraded.');

    legacy.campaignRoot = fullfile(tempdir, 'Drive_B', 'Campaign_2');
    normB = WorkspaceManager('normalize', legacy);
    assert(~strcmp(normA.indexDbPath, normB.indexDbPath), ...
        'Different campaign roots should map to different index DB files.');

    cacheDir = fullfile(prefdir, 'timeline_app_cache');
    if ~isfolder(cacheDir), mkdir(cacheDir); end
    testDb = fullfile(cacheDir, 'timeline_index_test_clear.sqlite');
    fid = fopen(testDb, 'w'); fwrite(fid, 'x'); fclose(fid);
    WorkspaceManager('setlast', fullfile(tempdir, 'dummy_workspace.json'));
    report = WorkspaceManager('clearcache');
    assert(any(strcmp(report.deletedPaths, testDb)), 'Expected test cache DB to be deleted.');
    assert(~isfile(testDb), 'Expected test cache DB file to be removed.');
    assert(~isfile(fullfile(cacheDir, 'last_workspace.json')), 'Expected last workspace pointer to be removed.');
end
