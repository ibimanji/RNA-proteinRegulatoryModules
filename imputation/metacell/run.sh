unset R_HOME

Rscript run_rds_wnn_nospatial.R \
  --input states_with_plaque_info.rds \
  --output_rds states_nospatial_metacells.rds \
  --membership_csv states_nospatial_membership.csv \
  --gamma 30 \
  --k_knn 30 \
  --total_dims 1:30 \
  --rb_dims 1:30 \
  --nfeatures 3000 \
  --min_cells_one_metacell 150 \
  --split_by type \
  --split_values all