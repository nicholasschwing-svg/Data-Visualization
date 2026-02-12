function tests = test_fridge_read_frame_cache
    tests = functiontests(localfunctions);
end

function testReadWithCache(testCase)
    tmp = tempname;
    mkdir(tmp);
    c = onCleanup(@() rmdir(tmp, 's'));

    rawPath = fullfile(tmp, 'sample_LWIR.raw');
    hdrPath = fullfile(tmp, 'sample_LWIR.hdr');

    fid = fopen(rawPath, 'w');
    fwrite(fid, uint16([1 2 3 4]), 'uint16');
    fclose(fid);

    fid = fopen(hdrPath, 'w');
    fprintf(fid, 'samples = 2\n');
    fprintf(fid, 'lines = 2\n');
    fprintf(fid, 'bands = 1\n');
    fprintf(fid, 'data type = 12\n');
    fprintf(fid, 'header offset = 0\n');
    fprintf(fid, 'interleave = bsq\n');
    fprintf(fid, 'byte order = 0\n');
    fclose(fid);

    hdr = readENVIHeader(hdrPath);
    hdrsMap = containers.Map({'LWIR'}, {hdr});
    filesMap = containers.Map({'LWIR'}, {rawPath});

    img1 = fridge_read_frame('LWIR', 1, hdrsMap, filesMap, 'UseCache', true);
    img2 = fridge_read_frame('LWIR', 1, hdrsMap, filesMap, 'UseCache', true);

    verifyEqual(testCase, img1, img2);
    verifyEqual(testCase, img1, uint16([1 2; 3 4]));

    clear c;
end
