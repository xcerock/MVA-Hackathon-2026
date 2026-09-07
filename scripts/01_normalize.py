#!/usr/bin/env python3
"""
Paso 1: preparacion del VCF.

Tres operaciones, todas necesarias antes de anotar:

  1. Conservar solo variantes con FILTER=PASS. El VCF original fue
     filtrado con GATK VariantFiltration usando filtros duros
     (QD<2, MQ<40, ReadPosRankSum<-8, FS>60, MQRankSum<-12.5).

  2. Restringir a cromosomas primarios. El VCF incluye contigs decoy
     (Un_JTFH01*) y el mitocondrial nombrado como "M", que no existen
     en la referencia de Ensembl y rompen la normalizacion.

  3. Separar sitios multialelicos en lineas independientes y realinear
     los indels a la izquierda. Sin esto, el cruce con bases de datos
     poblacionales y clinicas falla en los sitios multialelicos.

Uso:
    python 01_normalize.py [--input VCF] [--output VCF]
"""
from __future__ import annotations

import argparse
import shlex
from pathlib import Path

import config
from common import get_logger, run_command


def build_region_string(chromosomes: list[str]) -> str:
    return ",".join(chromosomes)


def normalize(input_vcf: Path, output_vcf: Path, reference: Path,
              chromosomes: list[str], threads: int, logger) -> None:
    regions = build_region_string(chromosomes)
    logger.info("Restringiendo a %d cromosomas primarios", len(chromosomes))

    cmd = (
        f"bcftools view -f PASS -r {shlex.quote(regions)} "
        f"{shlex.quote(str(input_vcf))} -Ou "
        f"| bcftools norm -m -any -f {shlex.quote(str(reference))} "
        f"-Oz -o {shlex.quote(str(output_vcf))} --threads {threads}"
    )
    run_command(cmd, logger, shell=True)
    run_command(["bcftools", "index", "-t", str(output_vcf)], logger)


def report_stats(vcf: Path, logger) -> None:
    import subprocess
    out = subprocess.run(["bcftools", "stats", str(vcf)],
                         capture_output=True, text=True)
    for line in out.stdout.splitlines():
        if line.startswith("SN"):
            logger.info(line.replace("\t", " "))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, default=config.RAW_VCF)
    parser.add_argument("--output", type=Path, default=config.NORM_VCF)
    parser.add_argument("--reference", type=Path, default=config.REFERENCE_FASTA)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    config.ensure_directories()
    logger = get_logger("01_normalize", config.LOGS_DIR / "01_normalize.log")

    for path in (args.input, args.reference):
        if not path.exists():
            logger.error("No existe: %s", path)
            raise SystemExit(1)

    normalize(args.input, args.output, args.reference,
              config.PRIMARY_CHROMOSOMES, args.threads, logger)
    report_stats(args.output, logger)
    logger.info("VCF normalizado en %s", args.output)


if __name__ == "__main__":
    main()
