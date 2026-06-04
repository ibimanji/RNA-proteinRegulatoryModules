import sys
import time
import argparse
import logging
import traceback
import json
import subprocess
from pathlib import Path

# 确保能找到 pystar 包
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import pystar
from pystar.infrastructure import load_config
from pystar.preprocessing import DataSanitizer
from pystar.registration import RegistrationEngine
from pystar.spot_finding import SpotFinder
from pystar.mining import SignalMiner
from pystar.decoding import Decoder


def setup_logger(fov_id, output_cfg, log_name, logger_name):
    export_base = output_cfg.export_directory or output_cfg.directory
    position_name = f"Position{int(fov_id):0{output_cfg.export_fov_digits}d}"
    log_dir = Path(export_base) / position_name
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / log_name

    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")

    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger, log_file


def run_main_pipeline(cfg, current_fov, logger):
    # === Stage 1: Preprocessing ===
    t0 = time.time()
    logger.info(">>> Stage 1: Preprocessing (Sanitization)...")
    sanitizer = DataSanitizer(cfg)
    sanitizer.sanitize_fov(current_fov)
    logger.info(f"    Done in {time.time() - t0:.2f}s")

    # === Stage 2: Registration ===
    t0 = time.time()
    logger.info(">>> Stage 2: Registration...")
    from pystar.io import ImageLoader

    loader = ImageLoader(cfg)
    data_xr = loader.load_fov(current_fov)

    reg_engine = RegistrationEngine(cfg)
    reg_engine.register_fov(data_xr, current_fov)
    logger.info(f"    Done in {time.time() - t0:.2f}s")

    del data_xr
    del loader

    # === Stage 3: Spot Finding ===
    t0 = time.time()
    logger.info(">>> Stage 3: Spot Finding...")
    finder = SpotFinder(cfg)
    finder.find_spots_in_fov(current_fov)
    logger.info(f"    Done in {time.time() - t0:.2f}s")

    # === Stage 4: Mining ===
    t0 = time.time()
    logger.info(">>> Stage 4: Signal Extraction...")
    miner = SignalMiner(cfg)
    miner.mine_fov(current_fov)
    logger.info(f"    Done in {time.time() - t0:.2f}s")

    # === Stage 5: Decoding ===
    t0 = time.time()
    logger.info(">>> Stage 5: Decoding...")
    decoder = Decoder(cfg)
    decoder.decode_fov(current_fov)
    logger.info(f"    Done in {time.time() - t0:.2f}s")


def run_if_registration(cfg, current_fov, logger, log_file):
    if not cfg.pipeline.if_registration.enable:
        logger.info(">>> Stage 6: IF / Protein Registration skipped: disabled in config.")
        return

    t0 = time.time()
    logger.info(">>> Stage 6: IF / Protein Registration...")

    if_cfg = cfg.pipeline.if_registration
    position_name = f"Position{int(current_fov):0{cfg.pipeline.output.export_fov_digits}d}"
    runtime_dir = Path(__file__).resolve().parents[1] / if_cfg.runtime_path

    payload = {
        "sample": if_cfg.sample,
        "user_dir": if_cfg.user_dir,
        "source_data_dir": if_cfg.source_data_dir,
        "registration_dir": if_cfg.registration_dir,
        "position_name": position_name,
        "protein_folder": if_cfg.protein_folder,
        "protein_round": if_cfg.protein_round,
        "protein_stains": if_cfg.protein_stains,
        "registration_channel": if_cfg.registration_channel,
    }

    config_json = json.dumps(payload).replace("'", "''")
    matlab_cmd = (
        f"addpath('{runtime_dir}'); "
        f"{if_cfg.entrypoint}('{config_json}')"
    )

    logger.info(f"Runtime dir: {runtime_dir}")
    logger.info(f"Position: {position_name}")

    with open(log_file, "a") as if_log_handle:
        subprocess.run(
            ["matlab", "-batch", matlab_cmd],
            stdout=if_log_handle,
            stderr=subprocess.STDOUT,
            check=True,
        )

    logger.info(f"    Done in {time.time() - t0:.2f}s")


def main():
    parser = argparse.ArgumentParser(description="PyStar Worker Node")
    parser.add_argument("--config", required=True, help="Path to experiment_config.yaml")
    parser.add_argument("--task_id", required=True, type=int, help="Slurm Array Task ID (1-based)")
    parser.add_argument(
        "--mode",
        choices=["main", "if", "all"],
        default="all",
        help="Which stages to run: main, if, or all",
    )
    args = parser.parse_args()

    try:
        cfg = load_config(args.config)
    except Exception as e:
        print(f"FATAL: Config load failed: {e}")
        sys.exit(1)

    fov_list = cfg.dataset.parsed_fovs
    total_jobs = len(fov_list)

    if args.task_id < 1 or args.task_id > total_jobs:
        print(f"FATAL: Task ID {args.task_id} is out of range [1, {total_jobs}]")
        sys.exit(1)

    current_fov = fov_list[args.task_id - 1]

    decode_logger, decode_log_file = setup_logger(
        current_fov,
        cfg.pipeline.output,
        "decode.log",
        f"pystar.decode.{current_fov}",
    )

    if_logger, if_log_file = setup_logger(
        current_fov,
        cfg.pipeline.output,
        "if.log",
        f"pystar.if.{current_fov}",
    )

    active_loggers = []
    if args.mode in {"main", "all"}:
        active_loggers.append(decode_logger)
    if args.mode in {"if", "all"}:
        active_loggers.append(if_logger)

    for logger in active_loggers:
        logger.info(f"{'=' * 40}")
        logger.info(" PyStar Worker Started")
        logger.info(f" Mode: {args.mode}")
        logger.info(f" Task ID: {args.task_id} / {total_jobs}")
        logger.info(f" Target FOV: {current_fov}")
        logger.info(f" Config source: {cfg.config_source_path}")
        logger.info(f" PyStar package: {Path(pystar.__file__).resolve()}")
        logger.info(f" Repo root: {REPO_ROOT}")
        logger.info(f" Raw data path: {cfg.dataset.raw_data_path}")
        logger.info(f" Filename pattern: {cfg.dataset.filename_pattern}")
        logger.info(f" Output directory: {cfg.pipeline.output.directory}")
        logger.info(f" Export directory: {cfg.pipeline.output.export_directory}")
        logger.info(f" Gene list: {cfg.codebook.gene_list}")
        logger.info(f"{'=' * 40}")

    start_time_global = time.time()

    try:
        if args.mode in {"main", "all"}:
            run_main_pipeline(cfg, current_fov, decode_logger)

        if args.mode in {"if", "all"}:
            run_if_registration(cfg, current_fov, if_logger, if_log_file)

        total_time = time.time() - start_time_global
        for logger in active_loggers:
            logger.info(f"{'=' * 40}")
            logger.info(f" SUCCESS: FOV {current_fov} mode={args.mode} complete.")
            logger.info(f" Total Time: {total_time / 60:.2f} minutes")
            logger.info(f"{'=' * 40}")

    except Exception:
        for logger in active_loggers:
            logger.error(f"CRITICAL FAILURE on FOV {current_fov} mode={args.mode}")
            logger.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
