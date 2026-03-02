function test_discovery_filter()
    root = fullfile(tempdir, ['discover_' char(java.util.UUID.randomUUID)]);
    mkdir(root);
    c = onCleanup(@() cleanup(root)); %#ok<NASGU>

    mkdir(fullfile(root, 'MX-20SW', 'MX-20SW', 'cal', 'CAL0209', 'hsic'));
    fid = fopen(fullfile(root, 'MX-20SW', 'MX-20SW', 'cal', 'CAL0209', 'hsic', 'sample.hsic'), 'w');
    fwrite(fid, 'x'); fclose(fid);

    mkdir(fullfile(root, 'HSI', 'run1'));
    fid = fopen(fullfile(root, 'HSI', 'run1', 'scene.hsic'), 'w');
    fwrite(fid, 'x'); fclose(fid);

    mkdir(fullfile(root, 'FAST_radiance', 'batch'));
    fid = fopen(fullfile(root, 'FAST_radiance', 'batch', 'capture.hsic'), 'w');
    fwrite(fid, 'x'); fclose(fid);

    % Should not be discovered because it has no evidence files.
    mkdir(fullfile(root, 'LWIR_V1K'));

    src = discoverDataSources(root, {'@tmp','@eaDir'}, 4);
    labels = lower(string({src.label}));

    assert(any(labels == "hsi"), 'Expected HSI source to be discovered.');
    assert(any(contains(labels, "fast")), 'Expected FAST source to be discovered.');
    assert(any(contains(labels, "mx-20sw")), 'Expected MX-20SW source to be discovered.');
    assert(~any(labels == "hsic"), 'Leaf hsic folders should not be separate sources.');
    assert(~any(labels == "lwir_v1k"), 'Source without .hsic evidence should not be listed.');
end

function cleanup(root)
    if isfolder(root)
        rmdir(root, 's');
    end
end
