function metadata_json = pystar_if_register_entry(config_json)

cfg = jsondecode(config_json);

input_path = fullfile(cfg.user_dir, cfg.sample, cfg.source_data_dir);
output_path = fullfile(cfg.user_dir, cfg.sample, cfg.registration_dir);

position_name = cfg.position_name;
protein_folder = cfg.protein_folder;
protein_round = cfg.protein_round;
protein_stains = string(cfg.protein_stains);

round1_path = fullfile(input_path, "round001", position_name);
round1_files = dir(fullfile(round1_path, "*.tif"));

round1_imgs = cell(numel(round1_files), 1);
for c = 1:numel(round1_files)
    round1_imgs{c} = pystar_load_multipage_tiff(fullfile(round1_files(c).folder, round1_files(c).name), 'uint8', 'uint8', false);
end

ref_img = round1_imgs{1};
for c = 2:numel(round1_imgs)
    ref_img = ref_img + round1_imgs{c};
end

protein_path = fullfile(input_path, protein_folder, position_name);
protein_files = dir(fullfile(protein_path, "*.tif"));

protein_imgs = cell(numel(protein_files), 1);
for c = 1:numel(protein_files)
    protein_imgs{c} = pystar_load_multipage_tiff(fullfile(protein_files(c).folder, protein_files(c).name), 'uint8', 'uint8', false);
end

moving_channel = cfg.registration_channel;
mov_img = protein_imgs{moving_channel};

if size(mov_img, 3) == 1
    params = pystar_dft_register_2d(ref_img, mov_img, false);
    for c = 1:numel(protein_imgs)
        protein_imgs{c} = uint8(pystar_dft_apply_2d(protein_imgs{c}, params, false));
    end
else
    params = pystar_dft_register_3d(ref_img, mov_img, false);
    for c = 1:numel(protein_imgs)
        protein_imgs{c} = uint8(pystar_dft_apply_3d(protein_imgs{c}, params, false));
    end
end

protein_output_dir = fullfile(output_path, protein_round);
pystar_save_protein_images(protein_output_dir, protein_imgs, position_name, protein_stains);

metadata = struct();
metadata.position_name = position_name;
metadata.output_dir = protein_output_dir;
metadata.shifts = params.shifts;
metadata.registration_channel = moving_channel;
metadata.n_channels = numel(protein_imgs);

metadata_json = jsonencode(metadata);
end
