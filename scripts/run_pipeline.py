#!/usr/bin/env python3
"""
Orquestador del pipeline completo.

Ejecuta los cinco pasos en orden, saltando los que ya tienen su salida
generada salvo que se pase --force. Registra el tiempo de cada etapa,
dato que el formulario de metodos pide explicitamente.

Uso:
    python run_pipeline.py                    # completo
    python run_pipeline.py --from 3           # desde el filtrado
    python run_pipeline.py --only 4 5         # solo dos pasos
    python run_pipeline.py --force            # reejecuta todo
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import config
from common import get_logger

SCRIPT_DIR = Path(__file__).resolve().parent

STEPS = [
    (1, "01_normalize.py", "Normalizacion del VCF", config.NORM_VCF),
    (2, "02_annotate.py", "Anotacion con VEP", config.VEP_TSV),
    (3, "03_filter_variants.py", "Filtrado de variantes raras", config.RARE_TSV),
    (4, "04_recessive_model.py", "Modelo recesivo", config.CANDIDATES_TSV),
    (5, "05_rank_submission.py", "Ranking y entrega", config.SUBMISSION_CSV),
]


def run_step(script: str, logger, extra_args: list[str] | None = None) -> float:
    cmd = [sys.executable, str(SCRIPT_DIR / script)] + (extra_args or [])
    start = time.time()
    result = subprocess.run(cmd)
    elapsed = time.time() - start
    if result.returncode != 0:
        logger.error("El paso %s fallo tras %.1f s", script, elapsed)
        raise SystemExit(result.returncode)
    return elapsed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--from", dest="start_at", type=int, default=1)
    parser.add_argument("--only", type=int, nargs="+")
    parser.add_argument("--force", action="store_true",
                        help="reejecuta aunque la salida exista")
    parser.add_argument("--proband-id", default=None)
    args = parser.parse_args()

    config.ensure_directories()
    logger = get_logger("pipeline", config.LOGS_DIR / "pipeline.log")

    selected = set(args.only) if args.only else {
        n for n, *_ in STEPS if n >= args.start_at}

    timings: list[tuple[str, float]] = []
    total_start = time.time()

    for number, script, description, output in STEPS:
        if number not in selected:
            continue
        if output.exists() and not args.force:
            logger.info("Paso %d (%s) omitido: %s ya existe",
                        number, description, output.name)
            continue

        logger.info("=" * 62)
        logger.info("Paso %d: %s", number, description)
        logger.info("=" * 62)

        extra: list[str] = []
        if number == 5 and args.proband_id:
            extra = ["--proband-id", args.proband_id]

        elapsed = run_step(script, logger, extra)
        timings.append((description, elapsed))
        logger.info("Completado en %.1f s (%.1f min)", elapsed, elapsed / 60)

    total = time.time() - total_start
    logger.info("=" * 62)
    logger.info("Resumen de tiempos")
    for description, elapsed in timings:
        logger.info("  %-34s %8.1f s", description, elapsed)
    logger.info("  %-34s %8.1f s (%.1f min)", "TOTAL", total, total / 60)


if __name__ == "__main__":
    main()
