# NCBI Primer-BLAST Notes

## Official behavior to preserve

- NCBI Primer-BLAST accepts an existing forward/reverse primer pair and can perform a specificity check. If only one primer is supplied, a template is required.
- NCBI recommends selecting the organism and the smallest database likely to contain the target. Broad `nt` searches without organism limits are slower and less precise for routine RT-qPCR validation.
- NCBI says a RefSeq accession is preferred as PCR template input because it helps Primer-BLAST identify the intended template and check specificity.
- The specificity check evaluates possible PCR products from forward-reverse, forward-forward, and reverse-reverse primer combinations, not just ordinary single-primer BLAST hits.
- Primer-BLAST result pages include a summary phrase such as "Primer pairs are specific to input template as no other targets were found..." and detailed sections named "Products on intended targets" and "Products on potentially unintended templates".

## Current form fields verified on 2026-06-06

Endpoint: `https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi`

Core fields:

- `CMD=request`
- `INPUT_SEQUENCE`: optional template accession or FASTA sequence
- `PRIMER_LEFT_INPUT`: forward primer sequence, 5'->3'
- `PRIMER_RIGHT_INPUT`: reverse primer sequence, 5'->3' on minus strand
- `SEARCH_SPECIFIC_PRIMER=on`
- `SEARCHMODE=0`
- `PRIMER_SPECIFICITY_DATABASE`: database value
- `ORGANISM`: organism name or taxid
- `TOTAL_PRIMER_SPECIFICITY_MISMATCH`: zero-based form value; displayed mismatch count minus 1
- `PRIMER_3END_SPECIFICITY_MISMATCH`: zero-based form value; displayed mismatch count minus 1
- `MISMATCH_REGION_LENGTH`: displayed 3' region length
- `TOTAL_MISMATCH_IGNORE`: displayed total mismatch ignore threshold
- `MAX_TARGET_SIZE`: max target amplicon size for specificity checking
- `HITSIZE`, `EVALUE`, `WORD_SIZE`: BLAST sensitivity controls
- `NUM_TARGETS_WITH_PRIMERS`: max targets to show for pre-designed primers
- `MAX_TARGET_PER_TEMPLATE`: max targets per sequence

Useful optional checkboxes:

- `ALLOW_TRANSCRIPT_VARIANTS=on`: gene-level rather than transcript-specific acceptance of same-gene splice variants
- `EXCLUDE_XM=on`: exclude predicted RefSeq transcripts
- `EXCLUDE_ENV=on`: exclude uncultured/environmental sequences
- `LOW_COMPLEXITY_FILTER=on`: avoid low-complexity regions

## RT primer verdict policy

Use these labels in reports:

- `pass`: clean Primer-BLAST specificity in the selected database, intended product present if a template was supplied, and no major local primer warnings.
- `warning`: specificity appears clean, but a design or evidence limitation remains.
- `fail`: off-target product(s), wrong target, no intended product, or Primer-BLAST reports an error that invalidates the check.
- `pending`: NCBI accepted the job but it has not finished.
- `blocked`: network, NCBI, or input errors prevented a valid check.

Always state the selected database and organism. A pass in `refseq_mrna` is not the same as a pass against genome assemblies or `nt`.
