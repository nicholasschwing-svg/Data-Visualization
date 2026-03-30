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
end
