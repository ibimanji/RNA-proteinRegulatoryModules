function [output_imgs, dims] = test_LoadImageStacks_zf(inputPath, sub_dir, round_ids, channel_ids, expected_z_slices, input_format, output_format, useGPU)
%TEST_LOADIMAGESTACKS_ZF Reviewable extracted loader for PyStar preprocessing runtime.

    warning('off', 'all');

    if nargin < 6 || isempty(input_format)
        input_format = 'uint8';
    end
    if nargin < 7 || isempty(output_format)
        output_format = 'uint8';
    end
    if nargin < 8
        useGPU = false;
    end
    global PYSTAR_LOADER_FILENAME_PATTERN PYSTAR_LOADER_FOV_ID;
    filename_pattern = PYSTAR_LOADER_FILENAME_PATTERN;
    fov_id = PYSTAR_LOADER_FOV_ID;
    if isempty(filename_pattern)
        filename_pattern = '';
    end
    if isempty(fov_id)
        fov_id = 0;
    end

    round_ids = reshape(double(round_ids), 1, []);
    channel_ids = reshape(double(channel_ids), 1, []);
    fov_id = double(fov_id);
    Nround = numel(round_ids);
    Nchannel = numel(channel_ids);

    if Nround == 0
        error('test_LoadImageStacks_zf:InvalidInput', 'round_ids must not be empty');
    end
    if Nchannel == 0
        error('test_LoadImageStacks_zf:InvalidInput', 'channel_ids must not be empty');
    end

    resolved_paths = cell(Nchannel, Nround);
    tried_patterns = cell(Nchannel, Nround);
    first_real_img = [];
    first_real_c_idx = 0;
    first_real_r_idx = 0;

    for r_idx = 1:Nround
        round_id = round_ids(r_idx);
        for c_idx = 1:Nchannel
            channel_id = channel_ids(c_idx);
            [curr_path, curr_patterns] = local_resolve_round_channel_path(inputPath, sub_dir, filename_pattern, round_id, channel_id, fov_id);
            resolved_paths{c_idx, r_idx} = curr_path;
            tried_patterns{c_idx, r_idx} = curr_patterns;

            if isempty(first_real_img) && ~isempty(curr_path)
                first_real_img = new_LoadMultipageTiff(curr_path, input_format, output_format, useGPU);
                first_real_img = local_adjust_z_depth(first_real_img, expected_z_slices);
                first_real_c_idx = c_idx;
                first_real_r_idx = r_idx;
            end
        end
    end

    if isempty(first_real_img)
        error( ...
            'test_LoadImageStacks_zf:NoInputTiffs', ...
            'No TIFF files found for any requested round/channel under %s/%s', ...
            inputPath, sub_dir ...
        );
    end

    inferred_z_slices = size(first_real_img, 3);
    output_imgs = zeros( ...
        size(first_real_img, 1), ...
        size(first_real_img, 2), ...
        inferred_z_slices, ...
        Nchannel, ...
        Nround, ...
        output_format ...
    );

    for r_idx = 1:Nround
        round_id = round_ids(r_idx);
        fprintf('Loading round %d...\n', round_id);

        for c_idx = 1:Nchannel
            channel_id = channel_ids(c_idx);
            curr_path = resolved_paths{c_idx, r_idx};

            if isempty(curr_path)
                fprintf('Missing round %d channel %d; filling with zeros. Tried: %s\n', ...
                    round_id, channel_id, strjoin(tried_patterns{c_idx, r_idx}, ' | '));
                curr_img = zeros( ...
                    size(output_imgs, 1), ...
                    size(output_imgs, 2), ...
                    inferred_z_slices, ...
                    output_format ...
                );
            elseif c_idx == first_real_c_idx && r_idx == first_real_r_idx
                curr_img = first_real_img;
            else
                curr_img = new_LoadMultipageTiff(curr_path, input_format, output_format, useGPU);
                curr_img = local_adjust_z_depth(curr_img, expected_z_slices);
            end

            output_imgs(:,:,:,c_idx,r_idx) = curr_img;
        end
    end

    dims = size(output_imgs);
end


function [curr_path, glob_patterns] = local_resolve_round_channel_path(inputPath, sub_dir, filename_pattern, round_id, channel_id, fov_id)
    round_dir = sprintf('round%03d', round_id);
    if isempty(filename_pattern)
        glob_patterns = {
            fullfile(inputPath, round_dir, sub_dir, sprintf('*ch%02d*.tif', channel_id)), ...
            fullfile(inputPath, round_dir, sub_dir, sprintf('*ch%d*.tif', channel_id)) ...
        };
    else
        glob_patterns = {
            local_build_glob_path(inputPath, filename_pattern, round_id, channel_id, fov_id, true), ...
            local_build_glob_path(inputPath, filename_pattern, round_id, channel_id, fov_id, false) ...
        };
    end

    matches = [];
    for pattern_idx = 1:numel(glob_patterns)
        curr_matches = dir(glob_patterns{pattern_idx});
        if ~isempty(curr_matches)
            matches = curr_matches;
            break;
        end
    end

    if isempty(matches)
        curr_path = '';
        return;
    end


    if numel(matches) ~= 1
        error( ...
            'test_LoadImageStacks_zf:AmbiguousInput', ...
            'Expected exactly one TIFF for round %d channel %d under pattern rooted at %s, found %d', ...
            round_id, channel_id, inputPath, numel(matches) ...
        );
    end

    curr_path = fullfile(matches(1).folder, matches(1).name);
end


function glob_path = local_build_glob_path(inputPath, filename_pattern, round_id, channel_id, fov_id, zero_pad_channel)
    pattern = char(filename_pattern);
    if zero_pad_channel
        channel_text = sprintf('%02d', channel_id);
    else
        channel_text = sprintf('%d', channel_id);
    end

    pattern = strrep(pattern, '{round:03d}', sprintf('%03d', round_id));
    pattern = strrep(pattern, '{round}', sprintf('%d', round_id));
    pattern = strrep(pattern, '{fov:03d}', sprintf('%03d', fov_id));
    pattern = strrep(pattern, '{fov}', sprintf('%d', fov_id));
    pattern = strrep(pattern, '{ch:02d}', sprintf('%02d', channel_id));
    pattern = strrep(pattern, '{ch}', channel_text);

    if local_is_absolute_path(pattern)
        glob_path = pattern;
    else
        glob_path = fullfile(inputPath, pattern);
    end
end


function tf = local_is_absolute_path(path_value)
    path_text = char(path_value);
    tf = startsWith(path_text, filesep) || ...
        ~isempty(regexp(path_text, '^[A-Za-z]:[\\/]', 'once'));
end


function adjusted = local_adjust_z_depth(curr_img, expected_z_slices)
    if nargin < 2 || isempty(expected_z_slices)
        adjusted = curr_img;
        return;
    end

    expected_z_slices = double(expected_z_slices);
    if size(curr_img, 3) == expected_z_slices
        adjusted = curr_img;
        return;
    end

    if size(curr_img, 3) > expected_z_slices
        error( ...
            'test_LoadImageStacks_zf:UnexpectedZDepth', ...
            'Input TIFF contains %d z-slices, expected at most %d', ...
            size(curr_img, 3), expected_z_slices ...
        );
    end

    adjusted = curr_img;
    adjusted(:,:,end+1:expected_z_slices) = zeros( ...
        size(curr_img, 1), ...
        size(curr_img, 2), ...
        expected_z_slices - size(curr_img, 3), ...
        class(curr_img) ...
    );
end
