function test_multiviewer_sync_utils()
    testNearest();
    testOverlap();
    testDesync();
end

function testNearest()
    base = datetime(2024,1,1,0,0,0);
    times = base + milliseconds([0 100 200 1000]);
    out = mv_resolve_panel_sample(times, base + milliseconds(180), struct('toleranceMs', 50));
    assert(strcmp(out.status, 'OK'));
    assert(out.idx == 3);

    out2 = mv_resolve_panel_sample(times, base + milliseconds(600), struct('toleranceMs', 50));
    assert(strcmp(out2.status, 'NO_SAMPLE_NEARBY'));
end

function testOverlap()
    base = datetime(2024,1,1,0,0,0);
    m = containers.Map('KeyType','char','ValueType','any');
    m('A') = base + milliseconds([0 1000 2000]);
    m('B') = base + milliseconds([40 1030 4000]);
    m('C') = base + milliseconds([5000]);

    [t, info] = mv_find_next_overlap(base + milliseconds(10), 1, m, 60, 2, 'ANY', 'A');
    assert(info.found);
    assert(abs(milliseconds(t - (base + milliseconds(1000)))) <= 1);
end

function testDesync()
    base = datetime(2024,1,1,0,0,0);
    tf = mv_is_desynced(base, {base, base + milliseconds(50)}, 100);
    assert(~tf);
    tf2 = mv_is_desynced(base, {base, base + milliseconds(150)}, 100);
    assert(tf2);
end
