#!/usr/bin/env python3
"""Submit or parse NCBI Primer-BLAST specificity checks for RT/qPCR primers."""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import uuid
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from html import unescape
from html.parser import HTMLParser
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin
from urllib.request import Request, urlopen


ENDPOINT = "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi"
ALLOWED_BASES = set("ACGT")
DATABASES = {
    "refseq_mrna",
    "refseq_rna",
    "refseq_representative_genomes",
    "PRIMERDB/genome_selected_species",
    "core_nt",
    "nt",
    "Custom",
}


class PrimerBlastParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.inputs: List[Dict[str, str]] = []
        self.text_parts: List[str] = []
        self._skip_depth = 0

    def handle_starttag(self, tag: str, attrs: List[Tuple[str, Optional[str]]]) -> None:
        if tag.lower() in {"script", "style"}:
            self._skip_depth += 1
            return
        if tag.lower() == "input":
            self.inputs.append({k.lower(): (v or "") for k, v in attrs})

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() in {"script", "style"} and self._skip_depth:
            self._skip_depth -= 1

    def handle_data(self, data: str) -> None:
        if self._skip_depth:
            return
        stripped = data.strip()
        if stripped:
            self.text_parts.append(stripped)


@dataclass
class PrimerMetrics:
    sequence: str
    length: int
    gc_percent: float
    wallace_tm: float
    max_homopolymer: int
    gc_3prime_last5: int


def normalize_primer(seq: str) -> str:
    seq = seq.strip().upper()
    seq = re.sub(r"^5['`-]*", "", seq)
    seq = re.sub(r"[-`']*3['`]*$", "", seq)
    seq = re.sub(r"[^A-Z]", "", seq)
    return seq


def validate_primer(seq: str, label: str) -> None:
    if not seq:
        raise ValueError(f"{label} primer is empty.")
    bad = sorted(set(seq) - ALLOWED_BASES)
    if bad:
        raise ValueError(
            f"{label} primer contains unsupported bases: {''.join(bad)}. "
            "Use only A/C/G/T for NCBI Primer-BLAST."
        )


def primer_metrics(seq: str) -> PrimerMetrics:
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    max_run = max((len(m.group(0)) for m in re.finditer(r"(A+|C+|G+|T+)", seq)), default=0)
    return PrimerMetrics(
        sequence=seq,
        length=len(seq),
        gc_percent=round(100.0 * gc / len(seq), 2),
        wallace_tm=float(2 * at + 4 * gc),
        max_homopolymer=max_run,
        gc_3prime_last5=seq[-5:].count("G") + seq[-5:].count("C"),
    )


def local_flags(fw: PrimerMetrics, rv: PrimerMetrics, product_min: int, product_max: int) -> List[str]:
    flags: List[str] = []
    for label, metrics in (("forward", fw), ("reverse", rv)):
        if metrics.length < 18 or metrics.length > 25:
            flags.append(f"{label} length {metrics.length} is outside the common 18-25 nt RT-qPCR range")
        if metrics.gc_percent < 40 or metrics.gc_percent > 60:
            flags.append(f"{label} GC {metrics.gc_percent:.1f}% is outside the common 40-60% range")
        if metrics.max_homopolymer > 4:
            flags.append(f"{label} has a homopolymer run of {metrics.max_homopolymer}")
        if metrics.gc_3prime_last5 == 0:
            flags.append(f"{label} has no G/C in the last 5 bases")
        if metrics.gc_3prime_last5 > 4:
            flags.append(f"{label} has {metrics.gc_3prime_last5} G/C bases in the last 5 bases")
    if abs(fw.wallace_tm - rv.wallace_tm) > 3:
        flags.append(
            f"rough Wallace Tm difference is {abs(fw.wallace_tm - rv.wallace_tm):.1f} C; verify with Primer-BLAST/Primer3 Tm"
        )
    if product_min < 50 or product_max > 500:
        flags.append("requested product range is unusual for RT-qPCR; verify assay intent")
    return flags


def parse_html(html: str) -> Tuple[PrimerBlastParser, str, Dict[str, str]]:
    parser = PrimerBlastParser()
    parser.feed(html)
    text = " ".join(parser.text_parts)
    text = re.sub(r"\s+", " ", text).strip()
    values: Dict[str, str] = {}
    for item in parser.inputs:
        name = item.get("name")
        if name and name not in values:
            values[name] = item.get("value", "")
    return parser, text, values


def strip_tags(fragment: str) -> str:
    fragment = re.sub(r"<(script|style)\b.*?</\1>", " ", fragment, flags=re.I | re.S)
    fragment = re.sub(r"<[^>]+>", " ", fragment)
    fragment = unescape(fragment)
    return re.sub(r"[ \t\r\f\v]+", " ", fragment).strip()


def section_text(html: str, title: str) -> str:
    pattern = re.compile(
        r'<div\s+class="prPairTl">\s*'
        + re.escape(title)
        + r"\s*</div>\s*<hr\s*/?>\s*<div\s+class=\"prPairDtl\">(.*?)</div>",
        flags=re.I | re.S,
    )
    match = pattern.search(html)
    return strip_tags(match.group(1)) if match else ""


def extract_refresh_url(html: str) -> Optional[str]:
    match = re.search(
        r"<meta[^>]+http-equiv\s*=\s*[\"']?refresh[\"']?[^>]+content\s*=\s*[\"'][^;]+;\s*URL=([^\"']+)",
        html,
        flags=re.I,
    )
    if not match:
        return None
    return urljoin(ENDPOINT, unescape(match.group(1).strip()))


def extract_job_key(html: str, url: Optional[str] = None) -> Optional[str]:
    for source in (url or "", html):
        match = re.search(r"job_key=([A-Za-z0-9_.-]+)", source)
        if match:
            return match.group(1)
        match = re.search(r"JOB ID:([A-Za-z0-9_.-]+)", source)
        if match:
            return match.group(1)
    return None


def extract_summary_text(html: str, values: Dict[str, str]) -> str:
    if values.get("PRIMER_RESULTS_INFO"):
        return unescape(values["PRIMER_RESULTS_INFO"])
    match = re.search(r"<dt>\s*Specificity of primers\s*</dt>\s*<dd>(.*?)</dd>", html, re.I | re.S)
    if match:
        return strip_tags(match.group(1))
    return ""


def numbered_values(values: Dict[str, str], prefix: str, count: int) -> List[str]:
    return [values.get(f"{prefix}_{idx}", "") for idx in range(count)]


def summarize_accessions(text: str, limit: int = 12) -> List[str]:
    hits: List[str] = []
    for line in re.split(r"\s{2,}|\n", text):
        line = line.strip()
        if not line or line.lower().startswith(("product length", "forward primer", "reverse primer", "template")):
            continue
        if re.search(r"\b[A-Z]{1,4}_?\d+(?:\.\d+)?\b", line) or " " in line:
            hits.append(line)
        if len(hits) >= limit:
            break
    return hits


def parse_result(
    html: str,
    *,
    forward: Optional[str],
    reverse: Optional[str],
    template: Optional[str],
    organism: Optional[str],
    database: Optional[str],
    product_min: int,
    product_max: int,
    result_url: Optional[str],
) -> Dict[str, object]:
    _, text, values = parse_html(html)
    refresh_url = extract_refresh_url(html)
    job_key = extract_job_key(html, result_url)
    result_info = extract_summary_text(html, values)
    exception = ""
    match = re.search(r"Exception error:\s*([^<.]+(?:\.[^<]*)?)", text, re.I)
    if match:
        exception = match.group(1).strip()

    pairs_count = int(values.get("PRIMER_PAIRS_NUMBER") or 0)
    intended = section_text(html, "Products on intended targets")
    unintended = section_text(html, "Products on potentially unintended templates")
    allowed = section_text(html, "Products on allowed targets")
    allowed_variants = section_text(html, "Products on allowed transcript variants")

    primer_pairs: List[Dict[str, object]] = []
    for idx in range(pairs_count):
        pair = {
            "index": idx + 1,
            "forward_sequence": values.get(f"FW_PRIMER_SEQ_{idx}", ""),
            "reverse_sequence": values.get(f"RV_PRIMER_SEQ_{idx}", ""),
            "forward_tm": values.get(f"FW_PRIMER_TM_{idx}", ""),
            "reverse_tm": values.get(f"RV_PRIMER_TM_{idx}", ""),
            "forward_gc": values.get(f"FW_PRIMER_GC_{idx}", ""),
            "reverse_gc": values.get(f"RV_PRIMER_GC_{idx}", ""),
            "forward_self_any": values.get(f"FW_PRIMER_SFCM_{idx}", ""),
            "forward_self_3prime": values.get(f"FW_PRIMER_SF3CM_{idx}", ""),
            "reverse_self_any": values.get(f"RV_PRIMER_SFCM_{idx}", ""),
            "reverse_self_3prime": values.get(f"RV_PRIMER_SF3CM_{idx}", ""),
            "product_length": values.get(f"PRODUCT_LENGTH_{idx}", ""),
        }
        primer_pairs.append(pair)

    supplied_fw = normalize_primer(forward or values.get("PRIMER_LEFT_INPUT", ""))
    supplied_rv = normalize_primer(reverse or values.get("PRIMER_RIGHT_INPUT", ""))
    metrics: Dict[str, object] = {}
    flags: List[str] = []
    if supplied_fw and supplied_rv:
        fw_metrics = primer_metrics(supplied_fw)
        rv_metrics = primer_metrics(supplied_rv)
        metrics = {"forward": asdict(fw_metrics), "reverse": asdict(rv_metrics)}
        flags = local_flags(fw_metrics, rv_metrics, product_min, product_max)

    for pair in primer_pairs:
        product_length = pair.get("product_length")
        if product_length:
            try:
                length = int(str(product_length))
            except ValueError:
                continue
            if length < product_min or length > product_max:
                flags.append(f"Primer-BLAST product length {length} is outside requested {product_min}-{product_max} bp")

    if exception:
        verdict = "blocked"
        reasons = [exception]
    elif refresh_url and "Primer-BLAST Results" not in text:
        verdict = "pending"
        reasons = ["NCBI job is still pending"]
    elif unintended.strip():
        verdict = "fail"
        reasons = ["Potentially unintended templates were reported"]
    elif pairs_count == 0 and not result_info:
        verdict = "blocked"
        reasons = ["Could not find a completed Primer-BLAST result in the HTML"]
    elif pairs_count == 0:
        verdict = "fail"
        reasons = ["No primer pair/product was reported"]
    elif re.search(r"specific to input template.*no other targets", result_info, re.I):
        verdict = "warning" if flags else "pass"
        reasons = flags[:] if flags else ["No unintended targets were found in the selected database"]
    elif not (template or values.get("INPUT_SEQUENCE")):
        verdict = "warning"
        reasons = ["No intended template was supplied; review listed products against the assay goal"]
        reasons.extend(flags)
    elif not unintended.strip() and intended.strip():
        verdict = "warning" if flags else "pass"
        reasons = flags[:] if flags else ["Intended product found and unintended-target section is empty"]
    else:
        verdict = "warning"
        reasons = ["Primer-BLAST result needs manual review"]
        reasons.extend(flags)

    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "verdict": verdict,
        "reasons": reasons,
        "job_key": job_key,
        "result_url": result_url,
        "refresh_url": refresh_url,
        "inputs": {
            "forward": supplied_fw,
            "reverse": supplied_rv,
            "template": template or values.get("INPUT_SEQUENCE", ""),
            "organism": organism or values.get("ORGANISM", ""),
            "database": database or values.get("PRIMER_SPECIFICITY_DATABASE", ""),
            "product_min": product_min,
            "product_max": product_max,
        },
        "ncbi_specificity_summary": result_info,
        "primer_pairs": primer_pairs,
        "local_metrics": metrics,
        "local_flags": flags,
        "sections": {
            "intended_targets": intended,
            "allowed_targets": allowed,
            "allowed_transcript_variants": allowed_variants,
            "potentially_unintended_templates": unintended,
        },
        "unintended_target_preview": summarize_accessions(unintended),
    }


def build_payload(args: argparse.Namespace) -> Dict[str, str]:
    forward = normalize_primer(args.forward or "")
    reverse = normalize_primer(args.reverse or "")
    validate_primer(forward, "forward")
    validate_primer(reverse, "reverse")
    if args.database not in DATABASES:
        raise ValueError(f"Unsupported database {args.database!r}. Use one of: {', '.join(sorted(DATABASES))}")
    if not args.organism and not args.allow_no_organism and args.database != "Custom":
        raise ValueError("Provide --organism or explicitly use --allow-no-organism.")

    payload = {
        "CMD": "request",
        "INPUT_SEQUENCE": args.template or "",
        "PRIMER_LEFT_INPUT": forward,
        "PRIMER_RIGHT_INPUT": reverse,
        "PRIMER_PRODUCT_MIN": str(args.product_min),
        "PRIMER_PRODUCT_MAX": str(args.product_max),
        "PRIMER_NUM_RETURN": "1",
        "PRIMER_MIN_TM": "57.0",
        "PRIMER_OPT_TM": "60.0",
        "PRIMER_MAX_TM": "63.0",
        "PRIMER_MAX_DIFF_TM": "3",
        "PRIMER_ON_SPLICE_SITE": str(args.junction_span),
        "SEARCH_SPECIFIC_PRIMER": "on",
        "SEARCHMODE": "0",
        "PRIMER_SPECIFICITY_DATABASE": args.database,
        "CUSTOM_DB": args.custom_db or "",
        "ORGANISM": args.organism or "",
        "ALLOW_NO_ORGANISM": "on" if args.allow_no_organism else "",
        "ENTREZ_QUERY": args.entrez_query or "",
        "TOTAL_PRIMER_SPECIFICITY_MISMATCH": str(args.total_mismatch - 1),
        "PRIMER_3END_SPECIFICITY_MISMATCH": str(args.three_prime_mismatch - 1),
        "MISMATCH_REGION_LENGTH": str(args.mismatch_region_length),
        "TOTAL_MISMATCH_IGNORE": str(args.total_mismatch_ignore),
        "MAX_TARGET_SIZE": str(args.max_target_size),
        "HITSIZE": str(args.hitsize),
        "UNGAPPED_BLAST": "on",
        "EVALUE": str(args.evalue),
        "WORD_SIZE": str(args.word_size),
        "MAX_CANDIDATE_PRIMER": "500",
        "NUM_TARGETS_WITH_PRIMERS": str(args.num_targets_with_primers),
        "MAX_TARGET_PER_TEMPLATE": str(args.max_target_per_template),
        "PRIMER_MIN_SIZE": "15",
        "PRIMER_OPT_SIZE": "20",
        "PRIMER_MAX_SIZE": "25",
        "PRIMER_MIN_GC": "20.0",
        "PRIMER_MAX_GC": "80.0",
        "GC_CLAMP": "0",
        "POLYX": "5",
        "PRIMER_MAX_END_STABILITY": "9",
        "PRIMER_MAX_END_GC": "5",
        "MONO_CATIONS": "50.0",
        "DIVA_CATIONS": "1.5",
        "CON_DNTPS": "0.6",
        "SALT_FORMULAR": "1",
        "TM_METHOD": "1",
        "CON_ANEAL_OLIGO": "50.0",
        "PRIMER_MISPRIMING_LIBRARY": "AUTO",
        "LOW_COMPLEXITY_FILTER": "on",
        "SHOW_SVIEWER": "on",
        "LINK_LOC": "",
        "NUM_DIFFS": "0",
        "NUM_OPTS_DIFFS": "0",
    }
    if args.span_intron:
        payload["SPAN_INTRON"] = "on"
        payload["MIN_INTRON_SIZE"] = str(args.min_intron_size)
        payload["MAX_INTRON_SIZE"] = str(args.max_intron_size)
    if args.allow_transcript_variants:
        payload["ALLOW_TRANSCRIPT_VARIANTS"] = "on"
    if args.exclude_predicted:
        payload["EXCLUDE_XM"] = "on"
    if args.exclude_env:
        payload["EXCLUDE_ENV"] = "on"
    if args.no_snp:
        payload["NO_SNP"] = "on"
    return payload


def encode_multipart(fields: Dict[str, str]) -> Tuple[bytes, str]:
    boundary = "----CodexPrimerBlast" + uuid.uuid4().hex
    chunks: List[bytes] = []
    for key, value in fields.items():
        chunks.append(f"--{boundary}\r\n".encode())
        chunks.append(f'Content-Disposition: form-data; name="{key}"\r\n\r\n'.encode())
        chunks.append(str(value).encode("utf-8"))
        chunks.append(b"\r\n")
    chunks.append(f"--{boundary}--\r\n".encode())
    return b"".join(chunks), boundary


def http_request(url: str, *, fields: Optional[Dict[str, str]], user_agent: str, timeout: int = 120) -> str:
    data = None
    headers = {"User-Agent": user_agent}
    if fields is not None:
        data, boundary = encode_multipart(fields)
        headers["Content-Type"] = f"multipart/form-data; boundary={boundary}"
    request = Request(url, data=data, headers=headers)
    try:
        with urlopen(request, timeout=timeout) as response:
            raw = response.read()
    except HTTPError as exc:
        raise RuntimeError(f"NCBI returned HTTP {exc.code}: {exc.reason}") from exc
    except URLError as exc:
        raise RuntimeError(f"Could not reach NCBI: {exc.reason}") from exc
    return raw.decode("utf-8", errors="replace")


def submit_and_poll(args: argparse.Namespace, payload: Dict[str, str], out_dir: Path) -> Tuple[str, Optional[str]]:
    user_agent = "rt-primer-blast-ncbi-skill/1.0"
    if args.email:
        user_agent += f" ({args.email})"

    html = http_request(ENDPOINT, fields=payload, user_agent=user_agent)
    (out_dir / "submit.html").write_text(html, encoding="utf-8")
    result_url = extract_refresh_url(html) or ENDPOINT
    deadline = time.time() + args.timeout_minutes * 60

    while time.time() < deadline:
        summary = parse_result(
            html,
            forward=args.forward,
            reverse=args.reverse,
            template=args.template,
            organism=args.organism,
            database=args.database,
            product_min=args.product_min,
            product_max=args.product_max,
            result_url=result_url,
        )
        if summary["verdict"] not in {"pending"}:
            return html, result_url
        wait_url = summary.get("refresh_url") or result_url
        if not isinstance(wait_url, str) or not wait_url:
            return html, result_url
        time.sleep(args.poll_seconds)
        html = http_request(wait_url, fields=None, user_agent=user_agent)
        result_url = wait_url

    raise TimeoutError(f"Primer-BLAST job did not finish within {args.timeout_minutes} minutes.")


def write_report(summary: Dict[str, object], out_dir: Path) -> None:
    lines: List[str] = []
    lines.append("# Primer-BLAST Specificity Report")
    lines.append("")
    lines.append(f"Verdict: {str(summary['verdict']).upper()}")
    lines.append("")
    lines.append("## Inputs")
    inputs = summary.get("inputs", {})
    if isinstance(inputs, dict):
        for key in ("forward", "reverse", "template", "organism", "database"):
            lines.append(f"- {key}: {inputs.get(key, '')}")
        lines.append(f"- expected_product_size: {inputs.get('product_min')} - {inputs.get('product_max')} bp")
    if summary.get("job_key"):
        lines.append(f"- ncbi_job_key: {summary['job_key']}")
    if summary.get("result_url"):
        lines.append(f"- result_url: {summary['result_url']}")
    lines.append("")
    lines.append("## Reasons")
    for reason in summary.get("reasons", []):
        lines.append(f"- {reason}")
    lines.append("")
    lines.append("## NCBI Specificity")
    lines.append(str(summary.get("ncbi_specificity_summary") or "No Primer-BLAST specificity summary was found."))
    lines.append("")
    pairs = summary.get("primer_pairs") or []
    lines.append("## Primer Pair Metrics")
    if pairs:
        for pair in pairs:  # type: ignore[assignment]
            if not isinstance(pair, dict):
                continue
            lines.append(
                "- Pair {index}: product {product_length} bp; F Tm {forward_tm}, GC {forward_gc}%; "
                "R Tm {reverse_tm}, GC {reverse_gc}%".format(**pair)
            )
    else:
        lines.append("- No primer pair metrics were parsed.")
    lines.append("")
    lines.append("## Potentially Unintended Targets")
    preview = summary.get("unintended_target_preview") or []
    if preview:
        for item in preview:
            lines.append(f"- {item}")
    else:
        lines.append("- None parsed from the Primer-BLAST unintended-template section.")
    lines.append("")
    lines.append("## Caveat")
    lines.append(
        "This is an in-silico Primer-BLAST specificity check. Primer efficiency and assay validity still need wet-lab confirmation."
    )
    (out_dir / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_out_dir(path: Optional[str]) -> Path:
    if path:
        out_dir = Path(path)
    else:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = Path(f"primer_blast_{stamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--forward", help="Forward primer sequence, 5'->3'.")
    parser.add_argument("--reverse", help="Reverse primer sequence as ordered, 5'->3' on minus strand.")
    parser.add_argument("--template", help="Target accession or FASTA sequence. Strongly recommended.")
    parser.add_argument("--organism", default="Homo sapiens", help="Organism name or taxid.")
    parser.add_argument("--allow-no-organism", action="store_true", help="Allow broad search without organism limit.")
    parser.add_argument("--database", default="refseq_mrna", help="Primer-BLAST specificity database.")
    parser.add_argument("--custom-db", help="Custom accession, assembly accession, or FASTA database input.")
    parser.add_argument("--entrez-query", help="Optional Entrez query restriction.")
    parser.add_argument("--product-min", type=positive_int, default=70, help="Expected minimum product size.")
    parser.add_argument("--product-max", type=positive_int, default=250, help="Expected maximum product size.")
    parser.add_argument("--max-target-size", type=positive_int, default=4000, help="Maximum target amplicon size.")
    parser.add_argument("--total-mismatch", type=positive_int, default=2, choices=range(1, 7))
    parser.add_argument("--three-prime-mismatch", type=positive_int, default=2, choices=range(1, 7))
    parser.add_argument("--mismatch-region-length", type=positive_int, default=5)
    parser.add_argument("--total-mismatch-ignore", type=positive_int, default=6, choices=range(1, 10))
    parser.add_argument("--hitsize", type=positive_int, default=50000)
    parser.add_argument("--evalue", type=positive_int, default=30000)
    parser.add_argument("--word-size", type=positive_int, default=7)
    parser.add_argument("--num-targets-with-primers", type=positive_int, default=1000)
    parser.add_argument("--max-target-per-template", type=positive_int, default=100)
    parser.add_argument("--allow-transcript-variants", action="store_true")
    parser.add_argument("--exclude-predicted", action="store_true")
    parser.add_argument("--exclude-env", action="store_true")
    parser.add_argument("--span-intron", action="store_true")
    parser.add_argument("--min-intron-size", type=positive_int, default=1000)
    parser.add_argument("--max-intron-size", type=positive_int, default=1000000)
    parser.add_argument("--junction-span", type=int, default=0, choices=(0, 1, 2), help="0 no preference, 1 must span exon junction, 2 may not span.")
    parser.add_argument("--no-snp", action="store_true", help="Avoid known SNPs when a suitable accession is supplied.")
    parser.add_argument("--poll-seconds", type=positive_int, default=20)
    parser.add_argument("--timeout-minutes", type=positive_int, default=30)
    parser.add_argument("--email", help="Email/contact string included in User-Agent.")
    parser.add_argument("--out-dir", help="Output directory for raw HTML, summary.json, and report.md.")
    parser.add_argument("--from-html", help="Parse an existing Primer-BLAST result HTML file instead of submitting.")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs and print payload without contacting NCBI.")
    return parser


def main(argv: Optional[Iterable[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    out_dir = make_out_dir(args.out_dir)

    try:
        if args.from_html:
            html = Path(args.from_html).read_text(encoding="utf-8", errors="replace")
            result_url = None
        else:
            payload = build_payload(args)
            if args.dry_run:
                fw = primer_metrics(payload["PRIMER_LEFT_INPUT"])
                rv = primer_metrics(payload["PRIMER_RIGHT_INPUT"])
                dry = {
                    "payload": payload,
                    "local_metrics": {"forward": asdict(fw), "reverse": asdict(rv)},
                    "local_flags": local_flags(fw, rv, args.product_min, args.product_max),
                }
                print(json.dumps(dry, indent=2, sort_keys=True))
                return 0
            html, result_url = submit_and_poll(args, payload, out_dir)
            (out_dir / "primer_blast_result.html").write_text(html, encoding="utf-8")

        summary = parse_result(
            html,
            forward=args.forward,
            reverse=args.reverse,
            template=args.template,
            organism=args.organism,
            database=args.database,
            product_min=args.product_min,
            product_max=args.product_max,
            result_url=result_url,
        )
        (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
        write_report(summary, out_dir)
        print(f"{summary['verdict'].upper()}: wrote {out_dir / 'summary.json'} and {out_dir / 'report.md'}")
        return 0 if summary["verdict"] in {"pass", "warning", "pending"} else 2
    except Exception as exc:
        error = {"generated_at": datetime.now(timezone.utc).isoformat(), "verdict": "blocked", "error": str(exc)}
        (out_dir / "summary.json").write_text(json.dumps(error, indent=2, sort_keys=True), encoding="utf-8")
        print(f"BLOCKED: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
