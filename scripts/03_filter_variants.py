#!/usr/bin/env python3
"""
Paso 3: filtrado de variantes raras y funcionalmente relevantes.

Reduce los millones de variantes anotadas a unos cientos aplicando tres
criterios independientes:

  Frecuencia: MAX_AF por debajo del umbral. Se conservan tambien las
  variantes ausentes de gnomAD, que son las mas interesantes en una
  enfermedad ultra rara. Nota: no se usa el campo AF del INFO del VCF,
  que en este dataset es la fraccion alelica dentro de la muestra con
  AN=2 y no una frecuencia poblacional; filtrar por el no elimina nada.

  Consecuencia: impacto HIGH o MODERATE, mas cualquier consecuencia que
  mencione splicing (VEP las clasifica como LOW pero pueden alterar el
  procesamiento del transcrito).

  Calidad del genotipo: profundidad, calidad y balance alelico. Este
  ultimo descarta clusters artefactuales, donde el caller produce muchas
  variantes consecutivas con fraccion alelica muy por debajo de 0,5.

Identificacion de variantes: el VCF fue llamado con dbSNP como sitios
conocidos, asi que VEP escribe el rsID en Uploaded_variation en lugar de
las coordenadas. La localizacion se hace por el campo Location mas el
alelo alterno, y solo se recurre a Uploaded_variation como respaldo.

Uso:
    python 03_filter_variants.py
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

import config
from common import (VepTable, get_logger, load_genotypes, parse_af,
                    parse_location, parse_uploaded_variation)

GENOTYPE_COLUMNS = ["gt", "dp", "gq", "ad_ref", "ad_alt", "vaf", "pid", "pgt"]


def passes_frequency(row: dict[str, str], max_af: float) -> bool:
    af = parse_af(row.get("MAX_AF", "-"))
    return af is None or af < max_af


def passes_consequence(row: dict[str, str]) -> bool:
    if row.get("IMPACT", "") in config.ACCEPTED_IMPACTS:
        return True
    consequence = row.get("Consequence", "")
    return any(k in consequence for k in config.SPLICE_KEYWORDS)


def locate_genotype(row: dict[str, str], index, reasons: Counter):
    """
    Resuelve el genotipo correspondiente a una fila de VEP.

    Estrategia principal: campo Location mas el alelo alterno.
    Respaldo: campo Uploaded_variation cuando trae coordenadas.
    """
    location = parse_location(row.get("Location", ""))
    allele = row.get("Allele", "")
    if location is not None:
        chrom, start, end = location
        genotype = index.lookup_by_location(chrom, start, end, allele)
        if genotype is not None:
            return genotype

    parsed = parse_uploaded_variation(row.get("Uploaded_variation", ""))
    if parsed is not None:
        chrom, pos, ref, alt = parsed
        genotype = index.lookup_exact(chrom, pos, ref, alt)
        if genotype is not None:
            return genotype

    reasons["genotipo_no_localizado"] += 1
    return None


def passes_quality(gt, reasons: Counter) -> bool:
    """Control de calidad del genotipo. Registra el motivo de descarte."""
    if gt.gt not in ("0/1", "1/0", "1/1"):
        reasons["genotipo_no_variante"] += 1
        return False
    if gt.dp is not None and gt.dp < config.MIN_DEPTH:
        reasons["profundidad_baja"] += 1
        return False
    if gt.gq is not None and gt.gq < config.MIN_GQ:
        reasons["calidad_genotipo_baja"] += 1
        return False
    if gt.ad_alt is not None and gt.ad_alt < config.MIN_ALT_READS:
        reasons["pocas_lecturas_alternas"] += 1
        return False

    vaf = gt.vaf
    if vaf is not None:
        if gt.is_het and not (config.HET_VAF_MIN <= vaf <= config.HET_VAF_MAX):
            reasons["balance_alelico_heterocigoto"] += 1
            return False
        if gt.is_hom_alt and vaf < config.HOM_VAF_MIN:
            reasons["balance_alelico_homocigoto"] += 1
            return False
    return True


def filter_variants(vep_tsv: Path, norm_vcf: Path, output_tsv: Path,
                    max_af: float, logger) -> None:
    logger.info("Cargando genotipos desde %s", norm_vcf)
    index = load_genotypes(norm_vcf)
    logger.info("Genotipos cargados: %d", len(index))

    table = VepTable(vep_tsv)
    reasons: Counter = Counter()
    seen: set[str] = set()
    kept = 0
    total = 0
    writer = None

    with open(output_tsv, "w", newline="") as out:
        for row in table:
            total += 1
            if total % 1_000_000 == 0:
                logger.info("Procesadas %d lineas, conservadas %d", total, kept)

            if not passes_frequency(row, max_af):
                reasons["frecuencia_alta"] += 1
                continue
            if not passes_consequence(row):
                reasons["consecuencia_no_relevante"] += 1
                continue

            gt = locate_genotype(row, index, reasons)
            if gt is None:
                continue
            if not passes_quality(gt, reasons):
                continue

            # Una variante puede aparecer varias veces si afecta a mas de
            # un gen. Se conserva cada combinacion variante-gen una vez.
            dedup_key = f"{gt.chrom}:{gt.pos}:{gt.ref}:{gt.alt}:{row.get('SYMBOL','')}"
            if dedup_key in seen:
                reasons["duplicado_variante_gen"] += 1
                continue
            seen.add(dedup_key)

            record = dict(row)
            record.update({
                "chrom": gt.chrom, "pos": gt.pos, "ref": gt.ref, "alt": gt.alt,
                "gt": gt.gt,
                "dp": gt.dp if gt.dp is not None else "",
                "gq": gt.gq if gt.gq is not None else "",
                "ad_ref": gt.ad_ref if gt.ad_ref is not None else "",
                "ad_alt": gt.ad_alt if gt.ad_alt is not None else "",
                "vaf": f"{gt.vaf:.3f}" if gt.vaf is not None else "",
                "pid": gt.pid, "pgt": gt.pgt,
            })

            if writer is None:
                fieldnames = (table.header()
                              + ["chrom", "pos", "ref", "alt"]
                              + GENOTYPE_COLUMNS)
                writer = csv.DictWriter(out, fieldnames=fieldnames,
                                        delimiter="\t", extrasaction="ignore")
                writer.writeheader()
            writer.writerow(record)
            kept += 1

    logger.info("Lineas totales analizadas: %d", total)
    logger.info("Variantes conservadas: %d", kept)
    logger.info("Motivos de descarte:")
    for reason, count in reasons.most_common():
        logger.info("  %-32s %d", reason, count)

    unlocated = reasons.get("genotipo_no_localizado", 0)
    if unlocated:
        logger.warning("%d variantes no pudieron localizarse en el VCF. "
                       "Revisar si la anotacion y el VCF proceden de la "
                       "misma normalizacion.", unlocated)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vep-tsv", type=Path, default=config.VEP_TSV)
    parser.add_argument("--vcf", type=Path, default=config.NORM_VCF)
    parser.add_argument("--output", type=Path, default=config.RARE_TSV)
    parser.add_argument("--max-af", type=float, default=config.MAX_POPULATION_AF)
    args = parser.parse_args()

    config.ensure_directories()
    logger = get_logger("03_filter", config.LOGS_DIR / "03_filter.log")

    for path in (args.vep_tsv, args.vcf):
        if not path.exists():
            logger.error("No existe: %s", path)
            raise SystemExit(1)

    filter_variants(args.vep_tsv, args.vcf, args.output, args.max_af, logger)
    logger.info("Resultado en %s", args.output)


if __name__ == "__main__":
    main()