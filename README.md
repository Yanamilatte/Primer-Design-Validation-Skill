# Primer Design Validation Skill

A Codex skill for validating and redesigning RT-PCR and RT-qPCR primers with NCBI Primer-BLAST.

This skill is intended for primer order-list review: it checks existing forward/reverse primer pairs, flags specificity risks, redesigns failed pairs when needed, re-validates replacement candidates, and keeps the final orderable list clean.

## What It Does

- Validates RT-PCR and RT-qPCR primer pairs against NCBI Primer-BLAST results.
- Reviews intended products, unintended templates, transcript variants, homologs, and pseudogene-like hits.
- Supports gene-level or isoform-specific interpretation depending on the assay goal.
- Handles cross-species background checks, such as human transgene primers in rat or mouse samples.
- Redesigns failed primer pairs using target RefSeq/FASTA templates and re-checks replacement candidates before marking them as orderable.
- Produces concise reports with pass/warning/fail decisions, product size, Tm/GC notes, and practical caveats.
- Cleans temporary Primer-BLAST HTML, job folders, helper files, and intermediate workbooks when cleanup is requested.

## Skill Name

Use this skill locally as:

```text
$primer-design-validation-skill
```

Example prompt:

```text
Use $primer-design-validation-skill to check whether these RT-qPCR primers are specific, then redesign and re-check any failed primer pairs.
```

## Recommended Workflow

1. Start from the vendor primer order list and keep the original forward/reverse oligo sequences in 5' to 3' orientation.
2. Identify the intended organism, target gene, assay goal, and RefSeq accession or target FASTA whenever possible.
3. Validate each primer pair with NCBI Primer-BLAST, preferably using `refseq_mrna` for RT-qPCR gene-expression assays.
4. Classify failed or risky primers before redesigning them.
5. Redesign only failed, blocked, or high-risk pairs from the intended RefSeq/FASTA target.
6. Re-check redesigned candidates before adding them to the final order list.
7. Mark unresolved targets as not orderable rather than silently including risky primer pairs.
8. Preserve only final reports/order workbooks and remove temporary files when cleanup is requested.

## Repository Contents

- `SKILL.md` - Core workflow and interpretation rules for the Codex skill.
- `agents/openai.yaml` - UI-facing metadata and default prompt.
- `scripts/primer_blast_check.py` - Helper script for submitting or parsing Primer-BLAST specificity checks.
- `references/ncbi-primer-blast-notes.md` - Notes on NCBI Primer-BLAST behavior, result interpretation, and specificity thresholds.

## Important Notes

Primer-BLAST validation is in-silico evidence, not wet-lab proof. Final RT-qPCR assays still need experimental confirmation, including efficiency testing, melt curve or gel checks, no-template controls, and appropriate biological controls.
