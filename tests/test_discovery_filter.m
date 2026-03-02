function test_discovery_filter()
    root = fullfile(tempdir, ['discover_' char(java.util.UUID.randomUUID)]);
    mkdir(root);
    c = onCleanup(@() cleanup(root)); %#ok<NASGU>

    mkdir(fullfile(root, 'MX-20SW', 'MX-20SW', 'cal', 'CAL0209', 'hsic'));
    mkdir(fullfile(root, 'HSI'));
    mkdir(fullfile(root, 'FAST_radiance'));

    src = discoverDataSources(root, {'@tmp','@eaDir'}, 4);
    labels = lower(string({src.label}));

    assert(any(labels == "hsi"), 'Expected HSI root source to be discovered.');
    assert(any(contains(labels, "fast")), 'Expected FAST source to be discovered.');
    assert(~any(labels == "hsic"), 'Leaf hsic folders should not be separate sources.');
end

function cleanup(root)
    if isfolder(root)
        rmdir(root, 's');
    end
end
