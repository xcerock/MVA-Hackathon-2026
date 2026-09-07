#!/usr/bin/env python3
"""
Paso 2: anotacion funcional con Ensembl VEP.

Convierte coordenadas en biologia interpretable: consecuencia sobre el
transcrito, frecuencias poblacionales de gnomAD (exomas y genomas),
predictores SIFT y PolyPhen, clasificaciones de ClinVar y nomenclatura
HGVS.

Notas sobre los flags:

  --everything        activa el conjunto completo de anotaciones, entre
                      ellas MAX_AF, que es la frecuencia maxima observada
                      en cualquier poblacion y el campo mas util para
                      filtrar por rareza.
  --pick_allele_gene  deja una consecuencia por gen y alelo en lugar de
                      una linea por transcrito. Sin esto la salida crece
                      un orden de magnitud.
  --hgvs --fasta      necesarios para obtener HGVSc y HGVSp. VEP desactiva
                      HGVS en modo offline si no se le pasa la referencia.

El entorno de VEP suele estar separado del principal por conflictos de
dependencias de Perl. Usar --vep-env para invocarlo via micromamba run.

Uso:
    python 02_annotate.py --vep-env vep
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import config
from common import get_logger, run_command


def build_vep_command(input_vcf: Path, output_tsv: Path, stats_file: Path,
                      cache_dir: Path, reference: Path, threads: int,
                      buffer_size: int, vep_env: str | None) -> list[str]:
    cmd: list[str] = []
    if vep_env:
        cmd += ["micromamba", "run", "-n", vep_env]
    cmd += [
        "vep",
        "--input_file", str(input_vcf),
        "--output_file", str(output_tsv),
        "--tab",
        "--force_overwrite",
        "--cache",
        "--offline",
        "--dir_cache", str(cache_dir),
        "--fasta", str(reference),
        "--species", "homo_sapiens",
        "--assembly", "GRCh38",
        "--everything",
        "--hgvs",
        "--pick_allele_gene",
        "--fork", str(threads),
        "--buffer_size", str(buffer_size),
        "--stats_file", str(stats_file),
    ]
    return cmd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, default=config.NORM_VCF)
    parser.add_argument("--output", type=Path, default=config.VEP_TSV)
    parser.add_argument("--stats", type=Path, default=config.VEP_STATS)
    parser.add_argument("--cache-dir", type=Path, default=config.VEP_CACHE)
    parser.add_argument("--reference", type=Path, default=config.REFERENCE_FASTA)
    parser.add_argument("--threads", type=int, default=config.THREADS)
    parser.add_argument("--buffer-size", type=int, default=config.VEP_BUFFER_SIZE)
    parser.add_argument("--vep-env", default="vep",
                        help="entorno de micromamba con VEP; vacio para usar el PATH actual")
    args = parser.parse_args()

    config.ensure_directories()
    logger = get_logger("02_annotate", config.LOGS_DIR / "02_annotate.log")

    if not args.input.exists():
        logger.error("Falta el VCF normalizado: %s", args.input)
        raise SystemExit(1)
    if not args.cache_dir.exists():
        logger.error("Falta el cache de VEP: %s", args.cache_dir)
        raise SystemExit(1)

    vep_env = args.vep_env or None
    if vep_env and shutil.which("micromamba") is None:
        logger.warning("micromamba no esta en el PATH; se invoca vep directamente")
        vep_env = None

    cmd = build_vep_command(args.input, args.output, args.stats,
                            args.cache_dir, args.reference, args.threads,
                            args.buffer_size, vep_env)
    logger.info("La anotacion de un WGS completo tarda entre 1 y 2 horas")
    run_command(cmd, logger)
    logger.info("Anotacion escrita en %s", args.output)


if __name__ == "__main__":
    main()
