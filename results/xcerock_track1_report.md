# Track 1: Variant Prediction

## Candidato principal

BUB1B (NM_001211.6), heterocigosis compuesta:

- chr15:40209701 T>G, c.2210T>G, p.Leu737Ter, stop_gained,
  gnomAD MAX_AF 9.98e-05, ClinVar pathogenic/likely_pathogenic (rs759242053)
- chr15:40220612 T>G, c.3006T>G, p.Asn1002Lys, missense,
  ausente de gnomAD, SIFT deleterious (0.01), PolyPhen probably_damaging (0.997)

Ambas en el transcrito MANE Select. La parada prematura en el residuo 737
elimina el dominio quinasa; la missense en el residuo 1002 cae dentro de
ese dominio. Patron compatible con MVA1: alelo nulo mas alelo hipomorfico.

## Metodo

VCF filtrado a PASS y cromosomas primarios, normalizado contra GRCh38 y
anotado con VEP 116.2 (cache offline, --everything --hgvs).
Filtrado por MAX_AF < 0,01, impacto HIGH o MODERATE, y control de calidad
de genotipo (DP, GQ, balance alelico). Modelo recesivo automatico:
homocigotos y pares de heterocigotos por gen, con resolucion de fase via
PID/PGT. Ranking por evidencia acumulada.

El pipeline es agnostico respecto al gen. Con los paneles de conocimiento
previo desactivados, BUB1B sigue en primera posicion (100 puntos frente a
52 del siguiente candidato).

## Divulgacion de uso de IA

Claude Opus 5 - Plan Pro 

## Repositorio

[[URL del repositorio](https://github.com/xcerock/MVA-Hackathon-2026.git)]
