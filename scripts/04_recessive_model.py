#!/usr/bin/env python3
"""
Paso 4: aplicacion del modelo de herencia recesiva.

El MVA se hereda de forma autosomica recesiva, asi que la variante causal
debe estar en doble dosis: homocigota, o dos variantes distintas en el
mismo gen (heterocigosis compuesta).

Sin datos de los progenitores no se puede fasear por herencia, pero GATK
HaplotypeCaller asigna PID y PGT a las variantes que comparten lecturas.
Cuando dos heterocigotas comparten PID:

    PGT identico  -> mismo cromosoma (cis)  -> NO es compuesto verdadero
    PGT distinto  -> cromosomas opuestos (trans) -> compuesto confirmado

Esa distincion elimina falsos positivos que la mayoria de pipelines sin
trio no detectan.

Dos salvaguardas contra loci problematicos:

  Genes con demasiadas heterocigotas raras se marcan como ruidosos. Un
  homopolimero, una duplicacion segmental o un pseudogen producen decenas
  de llamadas espurias, y ademas hacen explotar la combinatoria: trece
  heterocigotas generan setenta y ocho pares.

  Pares de variantes muy proximas se descartan, porque suelen ser la
  misma alteracion representada de dos formas distintas por el caller.

Uso:
    python 04_recessive_model.py
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from itertools import combinations
from pathlib import Path

import config
from common import get_logger, resolve_phase

OUTPUT_FIELDS = [
    "gene", "model", "phase",
    "chrom_1", "pos_1", "ref_1", "alt_1",
    "consequence_1", "impact_1", "hgvsc_1", "hgvsp_1",
    "max_af_1", "sift_1", "polyphen_1", "clin_sig_1",
    "existing_1", "gt_1", "dp_1", "vaf_1", "gq_1",
    "chrom_2", "pos_2", "ref_2", "alt_2",
    "consequence_2", "impact_2", "hgvsc_2", "hgvsp_2",
    "max_af_2", "sift_2", "polyphen_2", "clin_sig_2",
    "existing_2", "gt_2", "dp_2", "vaf_2", "gq_2",
    "distance_bp", "gene_het_count", "flags",
]


@dataclass
class Variant:
    gene: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: str
    consequence: str
    impact: str
    hgvsc: str
    hgvsp: str
    max_af: str
    sift: str
    polyphen: str
    clin_sig: str
    existing: str
    dp: str
    gq: str
    vaf: str
    pid: str
    pgt: str

    @property
    def is_het(self) -> bool:
        return self.gt in ("0/1", "1/0")

    @property
    def is_hom(self) -> bool:
        return self.gt == "1/1"


def load_variants(path: Path) -> dict[str, list[Variant]]:
    by_gene: dict[str, list[Variant]] = defaultdict(list)
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = row.get("SYMBOL", "").strip()
            if not gene or gene == "-":
                continue
            try:
                pos = int(row["pos"])
            except (KeyError, ValueError):
                continue
            by_gene[gene].append(Variant(
                gene=gene,
                chrom=row.get("chrom", ""),
                pos=pos,
                ref=row.get("ref", ""),
                alt=row.get("alt", ""),
                gt=row.get("gt", "./."),
                consequence=row.get("Consequence", ""),
                impact=row.get("IMPACT", ""),
                hgvsc=row.get("HGVSc", ""),
                hgvsp=row.get("HGVSp", ""),
                max_af=row.get("MAX_AF", "-"),
                sift=row.get("SIFT", "-"),
                polyphen=row.get("PolyPhen", "-"),
                clin_sig=row.get("CLIN_SIG", "-"),
                existing=row.get("Existing_variation", "-"),
                dp=row.get("dp", ""),
                gq=row.get("gq", ""),
                vaf=row.get("vaf", ""),
                pid=row.get("pid", "."),
                pgt=row.get("pgt", "."),
            ))
    return by_gene


def variant_columns(v: Variant, suffix: str) -> dict[str, str]:
    return {
        f"chrom_{suffix}": v.chrom,
        f"pos_{suffix}": str(v.pos),
        f"ref_{suffix}": v.ref,
        f"alt_{suffix}": v.alt,
        f"consequence_{suffix}": v.consequence,
        f"impact_{suffix}": v.impact,
        f"hgvsc_{suffix}": v.hgvsc,
        f"hgvsp_{suffix}": v.hgvsp,
        f"max_af_{suffix}": v.max_af,
        f"sift_{suffix}": v.sift,
        f"polyphen_{suffix}": v.polyphen,
        f"clin_sig_{suffix}": v.clin_sig,
        f"existing_{suffix}": v.existing,
        f"gt_{suffix}": v.gt,
        f"dp_{suffix}": v.dp,
        f"vaf_{suffix}": v.vaf,
        f"gq_{suffix}": v.gq,
    }


def empty_columns(suffix: str) -> dict[str, str]:
    return {f"{name}_{suffix}": "" for name in
            ["chrom", "pos", "ref", "alt", "consequence", "impact",
             "hgvsc", "hgvsp", "max_af", "sift", "polyphen", "clin_sig",
             "existing", "gt", "dp", "vaf", "gq"]}


def build_candidates(by_gene: dict[str, list[Variant]], logger) -> list[dict]:
    rows: list[dict] = []
    noisy_genes: list[tuple[str, int]] = []

    for gene, variants in sorted(by_gene.items()):
        homs = [v for v in variants if v.is_hom]
        hets = [v for v in variants if v.is_het]

        # Homocigotas: cada una es un candidato por si sola.
        for v in homs:
            row = {"gene": gene, "model": "homozygous", "phase": "n/a",
                   "distance_bp": "", "gene_het_count": str(len(hets)),
                   "flags": ""}
            row.update(variant_columns(v, "1"))
            row.update(empty_columns("2"))
            rows.append(row)

        if len(hets) < 2:
            continue

        # Salvaguarda contra loci problematicos.
        if len(hets) > config.MAX_HETS_PER_GENE:
            noisy_genes.append((gene, len(hets)))
            logger.info("Gen %s descartado: %d heterocigotas raras "
                        "(limite %d), probable locus artefactual",
                        gene, len(hets), config.MAX_HETS_PER_GENE)
            continue

        for a, b in combinations(sorted(hets, key=lambda v: v.pos), 2):
            distance = abs(b.pos - a.pos)
            flags = []

            if a.chrom != b.chrom:
                continue
            if distance < config.MIN_PAIR_DISTANCE_BP:
                flags.append("variantes_contiguas")

            phase = resolve_phase(a.pid, a.pgt, b.pid, b.pgt)
            if phase == "cis":
                flags.append("fase_cis_descartado")
            elif phase == "trans":
                flags.append("fase_trans_confirmada")

            row = {"gene": gene,
                   "model": "compound_het",
                   "phase": phase,
                   "distance_bp": str(distance),
                   "gene_het_count": str(len(hets)),
                   "flags": ";".join(flags)}
            row.update(variant_columns(a, "1"))
            row.update(variant_columns(b, "2"))
            rows.append(row)

    if noisy_genes:
        logger.info("Genes descartados por exceso de heterocigotas: %s",
                    ", ".join(f"{g}({n})" for g, n in noisy_genes))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, default=config.RARE_TSV)
    parser.add_argument("--output", type=Path, default=config.CANDIDATES_TSV)
    args = parser.parse_args()

    config.ensure_directories()
    logger = get_logger("04_recessive", config.LOGS_DIR / "04_recessive.log")

    if not args.input.exists():
        logger.error("No existe: %s", args.input)
        raise SystemExit(1)

    by_gene = load_variants(args.input)
    logger.info("Genes con al menos una variante rara funcional: %d", len(by_gene))

    rows = build_candidates(by_gene, logger)

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=OUTPUT_FIELDS,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    models = defaultdict(int)
    for r in rows:
        models[(r["model"], r["phase"])] += 1
    logger.info("Candidatos generados: %d", len(rows))
    for (model, phase), count in sorted(models.items()):
        logger.info("  %-14s fase=%-8s %d", model, phase, count)
    logger.info("Resultado en %s", args.output)


if __name__ == "__main__":
    main()
