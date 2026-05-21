cd /dfs/share/zenghuLab/2301920002/ADdecon/pystar

bash scripts/run_pystar.sh config/experiment_config_mix.yaml main
bash scripts/run_pystar.sh config/experiment_config_mix.yaml if

cd ../clustermap
bash run_cm_sweep_multiPos_final.sh params.tsv "Position078"
