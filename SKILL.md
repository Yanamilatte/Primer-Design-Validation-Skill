---
name: primer-design-validation-skill
description: Automate RT-PCR and RT-qPCR primer validation and redesign with NCBI Primer-BLAST. Use when Codex needs to check existing forward/reverse primer pairs, validate transcript or gene specificity, interpret Primer-BLAST HTML results, review off-target amplicons, redesign failed primer pairs, handle NCBI hit-review pages, cross-check species background specificity, or produce a concise primer validity/redesign report with Tm, GC, product length, and in-silico caveats.
---

# RT Primer-BLAST NCBI

## Workflow

1. Collect the primer pair as ordered from the vendor: forward primer 5'->3' and reverse primer 5'->3'. Do not reverse-complement the reverse primer unless the user explicitly says their sequence is not the primer oligo sequence.
2. Get the intended organism or taxonomy ID. Prefer a target RefSeq accession or FASTA sequence when the intended RT-PCR product is known; Primer-BLAST can run with only a primer pair, but a strict "valid for this target" conclusion requires an intended template.
3. Decide whether the assay is isoform-specific or gene-level. For routine RT-qPCR gene-expression assays, use `--allow-transcript-variants` and state that same-gene transcript variants are allowed targets; for isoform-specific assays, do not allow variants unless the user explicitly requests gene-level detection.
4. Prefer the smallest relevant database. For RT-qPCR against mature transcripts, start with `refseq_mrna` plus organism. Use `refseq_rna` for noncoding RefSeq RNA, a genome database for genomic DNA risk, or `core_nt`/`nt` only when broad coverage is necessary.
5. Run the bundled script from the skill directory:

```bash
python3 scripts/primer_blast_check.py \
  --forward AGAGCTACGAGCTGCCTGAC \
  --reverse AGCACTGTGTTGGCGTACAG \
  --template NM_001101.5 \
  --organism "Homo sapiens" \
  --database refseq_mrna \
  --allow-transcript-variants \
  --out-dir primer_blast_actb
```

6. If NCBI has already been run manually, parse the saved result page instead of submitting again:

```bash
python3 scripts/primer_blast_check.py \
  --from-html primer_blast_result.html \
  --out-dir primer_blast_parsed
```

7. If the pair fails, classify the failure before redesign:
   - High-risk: unintended products in the qPCR window, same/near-same product length, wrong locus/species, or target/background species cross-amplification when species-specific detection is required.
   - Lower-risk but not clean: only long off-target products (for example >500 bp) or predicted transcripts outside the qPCR window.
   - Blocked: NCBI hit-review page, timeout, no completed result, or missing target accession.
8. Redesign failed high-risk or blocked primer pairs when the user asks for a usable order list. Use NCBI Primer-BLAST design with the intended RefSeq/FASTA template, `refseq_mrna`, organism, expected product size 80-180 bp for qPCR, Tm near 60 C, and gene-level transcript variants allowed when appropriate. Keep redesigned candidates separate from original primers until they pass validation.
9. If NCBI returns a hit-review page, select only intended/allowed targets:
   - Allow same-gene transcript variants for gene-level RT-qPCR.
   - Do not allow unrelated homologs, pseudogenes, "like" genes, or wrong-species targets merely to force a pass.
   - After selecting allowed targets, resubmit and parse the final Primer-BLAST result.
10. For species-specific assays, add a background check in the species that must not amplify. The intended target species should pass, and the host/background species should show no target templates in the selected database.
11. Output a concise report and, for order-list tasks, a clean final table/workbook with only passable primer pairs marked as orderable. Mark pairs that cannot be cleanly redesigned as "do not order" or "needs vendor/manual redesign"; do not silently include failed candidate sequences.
12. Clean temporary raw HTML, job directories, helper scripts, and intermediate workbooks when the user asks for cleanup. Preserve only user-facing final reports/order lists and the original source file unless the user asks otherwise.
13. Report the verdict as in-silico evidence, not wet-lab proof. RT-qPCR efficiency still requires standard curve, melt curve or gel confirmation, and negative controls.

## Interpretation Rules

- Mark `pass` only when Primer-BLAST reports the pair is specific to the input template or shows an intended target product with no potentially unintended templates in the selected database.
- Mark `warning` when specificity is clean but the run has important limitations: no template accession, product size outside the requested RT-qPCR window, local primer metrics are weak, or the database does not cover the organism well.
- Mark `fail` when there are unintended target products, no intended product, multiple products inconsistent with the assay goal, or the hit is the wrong locus/species.
- For transcript-specific assays, treat same-gene splice variants as off-targets unless the user requested gene-level detection. For gene-level RT assays, use `--allow-transcript-variants` and state that the result is gene-specific rather than isoform-specific.
- Treat same-size or qPCR-window off-target products as strong redesign signals. Long off-target products outside normal qPCR conditions may be reported as lower practical risk, but they are still not strict Primer-BLAST passes.
- For cross-species expression validation, a primer pair that passes the intended species but amplifies the host/background species is not usable for species-specific detection.
- For housekeeping genes with many pseudogene-like or paralog hits, prefer a clean redesigned pair or suggest alternative reference genes rather than forcing a failed design into the order list.
- Do not run large batches against NCBI Primer-BLAST. For many primer pairs, advise local BLAST/Primer3-style screening or space requests out manually.

## Script Options

- `--forward` and `--reverse`: required for new NCBI submissions.
- `--template`: RefSeq/GenBank accession or FASTA sequence for the intended template. Strongly recommended for an actionable RT primer verdict.
- `--organism`: organism name or taxid; default is `Homo sapiens`. Use `--allow-no-organism` only when broad off-target discovery is intentional.
- `--database`: one of `refseq_mrna`, `refseq_rna`, `refseq_representative_genomes`, `PRIMERDB/genome_selected_species`, `core_nt`, `nt`, or `Custom`.
- `--product-min` and `--product-max`: expected product-size range, default `70-250` for RT-qPCR.
- `--max-target-size`: largest amplicon Primer-BLAST should consider for specificity; default `4000`.
- `--dry-run`: validate inputs and show the NCBI payload without submitting.
- `--from-html`: parse a saved Primer-BLAST HTML result and generate `summary.json` plus `report.md`.

## Redesign Playbook

Use this when a user asks to "redesign", "重新设计", "fix failed primers", or wants a final order list after validation.

1. Start from target RefSeq accessions, not the failed primer sequences.
2. Design with Primer-BLAST using product 80-180 bp for SYBR RT-qPCR unless the user has a different assay window.
3. Return multiple candidates internally, but mark only validated candidates as orderable.
4. If a target repeatedly fails because of homologs/pseudogenes, mark it as not orderable and explain the biological reason.
5. For mixed-species systems, validate both the intended target species and the host/background species.
6. Keep final outputs clean: one orderable table plus a short report with blocked/failed items and wet-lab validation reminders.

## References

Read `references/ncbi-primer-blast-notes.md` when changing database choices, result parsing, or specificity thresholds.
