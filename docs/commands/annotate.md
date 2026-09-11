# karyoscope annotate

Annotate sequences against a KaryoScope database, producing per-feature-set BED tracks.

## Synopsis

```
karyoscope annotate -i INPUT [OPTIONS]
```

## Description

`karyoscope annotate` assigns every k-mer in the input to a feature in a single alignment-free pass, by querying the database's index. The backend follows the database: an HKS-index database (e.g. `HKS_human_CHM13_v2`) is queried with `hks`, and a KMC-index database with the bundled `get_featureIDs` helper. It produces one BED per feature set, and by default writes BOTH a "presmoothed" (raw) and a "smoothed" (hierarchy-smoothed) BED for each feature set. k-mers that are not present in the index render as `novel`. Input may be FASTA, FASTQ, BAM or CRAM. The expensive k-mer-query step is resumable across reruns, so an interrupted (for example, OOM-killed) run can resume straight into smoothing.

## Options

| Option | Description |
| --- | --- |
| `-i`, `--input FILE` | Input sequence file. Accepts FASTA (`.fasta`/`.fa`/`.fna`, plain or `.gz`), FASTQ (`.fastq`/`.fq`, plain or `.gz`), BAM (`.bam`) or CRAM (`.cram`). BAM/CRAM inputs are converted with `samtools fasta` (requires `samtools` on PATH); CRAM also requires `--reference`. **[required]** |
| `--reference FILE` | Reference FASTA a CRAM input was aligned against. **Required for `.cram`**, ignored otherwise. See [CRAM input](#cram-input). |
| `--query-names` / `--no-query-names` | Identify output sequences by name rather than ordinal rank. Assemblies default to names; **refused for read-level input**, which always uses ranks. See [Paired-end reads](#paired-end-reads). |
| `--query-names-sidecar` | Also write `<outdir>/<input>.query_names.txt.gz` for each input: every record's name, one per line, in the order the query assigns ranks (line N+1 is rank N). This is how rank-identified read-level output is joined back to read names. For BAM/CRAM on the HKS backend it is teed off the decode `annotate` performs anyway, so it is nearly free; for FASTA/FASTQ input (either backend) it is one `awk` pass over the file, and for BAM/CRAM on the KMC backend a separate decode. See [Query names](#query-names). |
| `-o`, `--outdir DIRECTORY` | Directory to write output BEDs into. Default: same directory as `--input`. |
| `--db TEXT` | Database id to use (e.g., `KS_human_CHM13_v2`). Default: the unique installed database if there's exactly one. |
| `--db-root DIRECTORY` | Override the database root directory (default: `$KARYOSCOPE_DB` or `~/.karyoscope/db/`). |
| `--feature-set TEXT` | Restrict output to this feature set. Repeatable. Default: all feature sets declared in the database's manifest. |
| `-t`, `--threads INTEGER` | Threads for both k-mer querying and smoothing. `0` means auto-detect. [default: `0`] |
| `--k INTEGER` | Query k-mer length. Defaults to the database's k. Only a variable-k HKS index (built with `karyoscope build --variable-k`) accepts a value other than its k — use it for a k-sweep; on a fixed-k index any other value is an error. A database built with priorities is always fixed-k. Outputs are tagged `.k<k>` so runs into one directory don't collide. |
| `--smooth` / `--no-smooth` | Produce the hierarchy-smoothed BED in addition to the presmoothed BED. [default: `smooth`] |
| `--keep-presmoothed` / `--no-keep-presmoothed` | Keep the presmoothed BED. Pass `--no-keep-presmoothed` to write only the smoothed output. [default: `keep-presmoothed`] |
| `--keep-intermediates` | Keep the combined `.featureIDs.bed` from the C++ step (useful for debugging). |
| `--force` | Regenerate the combined intermediate even if a complete one already exists. By default a rerun reuses a verified combined BED left by a previous (e.g. OOM-killed) run and skips the `get_featureIDs` step, resuming straight into smoothing. A partial file from a killed run is never reused regardless of this flag. |
| `--bgzip` / `--no-bgzip` | bgzip the per-feature-set output BEDs. Pass `--no-bgzip` to write plain `.bed` files. Note this shrinks the final output but not peak disk usage: compression runs after every BED has been written in full. [default: `bgzip`] |
| `--no-resource-check` | Skip the up-front disk **and memory** checks. Disk: the output footprint on `--outdir`, see [Disk space](#disk-space). Memory: HKS databases must hold their index resident (~10 GB for a human database), so an under-allocated run is killed by the OS mid-query; the check reads the requirement from the index files themselves. |
| `--preserve-order` / `--no-preserve-order` | Write output BEDs with sequences in the same order as the input. Pass `--no-preserve-order` for the fastest path when order doesn't matter downstream (typically read data). [default: `preserve-order`] |
| `-h`, `--help` | Show this message and exit. |

## Examples

```bash
# Annotate an assembly with all feature sets (default), plain-text BEDs
karyoscope annotate --input hg002v1.1.fasta.gz --outdir results/ --threads 16 --no-bgzip

# Only the chromosome feature set
karyoscope annotate -i asm.fa --feature-set chromosome -o results/

# Read-level input (FASTQ): faster writes when output order doesn't matter
karyoscope annotate -i reads.fastq.gz -o results/ --no-preserve-order

# Read-level input, keeping the rank -> read-name mapping next to the output
karyoscope annotate -i tumor.cram --reference GRCh38.fasta -o results/ --query-names-sidecar
```

## CRAM input

CRAM stores bases as a diff against the reference used for alignment, so it cannot be
decoded without that same FASTA. Pass it with `--reference`:

```bash
karyoscope annotate -i tumor.cram --reference GRCh38.fasta \
                    --db HKS_human_CHM13_cytoband --feature-set cytoband -o results/
```

`--reference` is the **alignment** reference and has nothing to do with `--db`. `annotate`
is alignment-free — it only ever sees read sequences — so annotating GRCh38-aligned reads
against a CHM13 database is normal, not a mismatch.

Omitting `--reference` is refused up front rather than left to htslib, which would
otherwise resolve each contig's M5 checksum through `$REF_PATH`/`$REF_CACHE` and could
decode against a *different* build of the same genome, yielding plausible but incorrect
sequence.

**What gets converted.** `samtools fasta -F 0x900 -N` — primary records only, with the
`/1`,`/2` mate suffix forced on. That is exactly one full-length copy of every read:
each read has one primary record carrying full-length SEQ, whereas supplementary records
are hard-clipped slices of reads the primary already supplied. Unmapped reads and
duplicates are **kept** — they are genuine distinct reads. Minus-strand records are
reverse-complemented back to sequencing orientation.

**Scratch space.** The KMC backend streams the conversion with no temp file. The HKS
backend cannot — `hks lookup` needs a seekable path — so it materialises a temp FASTA
next to the output first (the system tempdir on cluster nodes is often small or
RAM-backed). **That file is full size**: a 64x human WGS CRAM measured 28.9 GB in and
254 GB out, so budget that space on the output filesystem on top of the BEDs.

## Paired-end reads

`annotate` has no concept of pairing: each mate is an independent query sequence.
What matters is whether the output can still be joined back into fragments afterwards.

Read-level output is identified by **ordinal rank**, and `--query-names` is *refused*
for read-level input. This is not a style preference — `hks` collects every sequence
name into memory before querying (`load_seq_names` builds a `Vec<String>` of the lot),
which is fine for an assembly's few thousand contigs and fatal for reads. Measured on a
64x human WGS CRAM: 1.315 G names at ~68 bytes each is ~90 GB on top of the ~12.3 GB
index, and the run is OOM-killed roughly 20 minutes in, having already spent 15 of them
decoding. Raising the memory limit buys nothing but longer identifiers.

Nothing is lost. Rank N is the Nth record of the query file, so the mapping back to read
names is a property of the input, and `--query-names-sidecar` records it next to the output
(see [Query names](#query-names)). Mates group by stripping the `/1`,`/2` suffix — which is
why `-N` is forced during the alignment decode. In a coordinate-sorted BAM/CRAM the two
mates of a fragment are nowhere near each other in rank, so pairing needs the names; in
split R1/R2 FASTQ files mate *i* of one file is mate *i* of the other, and the ranks pair
by themselves.

Note the coordinates are k-mer offsets **within each read**, not genomic positions: a
151 bp read queried at k=31 spans 0..121.

## Query names

`--query-names-sidecar` writes `<outdir>/<input>.query_names.txt.gz` for each input: one
record name per line, in the order the query assigns ranks, so line N+1 is rank N. The name
is the record header up to its first whitespace (`SRR123.7` from `@SRR123.7 1 length=151`),
with the `/1`,`/2` mate suffix for BAM/CRAM input. The file is plain `gzip`, not bgzip. It is
named from the input alone — not the database — because the mapping does not depend on
what the reads were annotated against, so one CRAM run against two databases yields one
sidecar.

What it costs depends on where the names come from:

| Input | Backend | How the names are produced |
|---|---|---|
| BAM / CRAM | HKS | Teed off the `samtools fasta` decode `annotate` runs anyway. Nearly free; regenerating the mapping afterwards is a second full decode (~25 min on a 56 GB CRAM). |
| FASTA / FASTQ | either | One `awk` pass over the file (`gzip -dc` in front of it for `.gz`), picking header lines by `>` for FASTA and by position — every fourth line — for FASTQ. One extra read of the input, small next to a lookup that reads it once per feature set; decompression is the whole cost. |
| BAM / CRAM | KMC | A separate `samtools fasta` pass — a second decode, because `get_featureIDs` is fed by a streaming pipe that leaves no seekable copy behind. |

Only `hks` emits ranks; `get_featureIDs` (the KMC backend) always writes names into the BED, so there the sidecar is written for consistency rather than out of necessity.

The FASTQ scan takes four lines per record without parsing them, which is safe because it is exactly how `hks` reads FASTQ: a record that spans more lines is a parse error there, not a rank. Headers are picked by position rather than by a leading `@`, since a quality line may begin with `@`.

Without the flag the mapping is reproducible by hand, since the query file is a
deterministic function of the input:

```bash
# BAM / CRAM: the same decode annotate performs
samtools fasta -F 0x900 -N --reference GRCh38.fasta tumor.cram | grep '^>'
# FASTA
grep '^>' asm.fa
# FASTQ (four lines per record, as hks reads it)
gzip -dc reads.fastq.gz | awk 'NR % 4 == 1'
```

## Output

For each feature set, `annotate` writes up to two BEDs:

- Presmoothed (raw): `<input>.<dbid>.<feature_set>.presmoothed.bed[.gz]`
- Smoothed (hierarchy-smoothed): `<input>.<dbid>.<feature_set>.smoothed.bed[.gz]`

With `--query-names-sidecar`, one more file per input, independent of the database: `<input>.query_names.txt.gz` (see [Query names](#query-names)).

Each BED's 4th column is the human-readable feature name. k-mers absent from the index render as `novel`. The `.gz` suffix is present unless `--no-bgzip` is passed; `--no-bgzip` keeps BEDs as plain text.

Smoothing promotes short noisy intervals (especially short `novel` runs flanked by specific features) to the lowest common ancestor of their flankers in the feature set's hierarchy.

For human-scale inputs, use at least 16 threads. Memory to request depends on the backend and input shape. With the **HKS backend**, the peak is set by the index rather than the input — the shared base index plus one feature set's labeling at a time — so it is ~10 GB whether you annotate a single haplotype or a combined diploid assembly such as HG002 v1.1 (request ≥ 16 GB). The **KMC backend**'s largest process peaks at ~36 GB on a haplotype and ~46 GB on a combined diploid, and its smoothing workers can collectively need as much again (request ≥ 96 GB). At `-t 16`, a haplotype takes ~7–9 min (HKS) or ~13–14 min (KMC); HG002 v1.1 combined takes ~14–16 min (HKS) or ~23–27 min (KMC).

## Progress output

`annotate` reports what it is doing as it goes, so a long run is distinguishable from a hung one:

```
Annotating hg002v1.1.fasta.gz against HKS_human_CHM13_v2
  6 feature set(s), 16 thread(s), ~34 GB estimated output
  [1/6] chromosome    4m05s
  [2/6] region        3m45s
  ...
  bgzip (12 file(s))  1m31s
Wrote:
  ...
```

The milestone shape follows the backend. HKS queries each feature set in sequence, so it reports `[i/N]` per set. KMC runs one combined query and then smooths every set in a single streaming pass, so it reports named phases instead (`k-mer query`, `smoothing 6 feature set(s)`).

Pass `-q` (before the subcommand: `karyoscope -q annotate ...`) to suppress the narration; the closing `Wrote:` block still prints, since that is the command's result. Detailed per-step timings remain on the logging channel — `karyoscope -v annotate ...` for those.

## Disk space

Output is large. Budget roughly **0.8 GB per feature set per Gbp of input**, plus one intermediate the size of the largest single feature set's presmoothed BED:

| Input | Feature sets | Output BEDs (uncompressed) | Peak including intermediate |
|---|---|---|---|
| Single human haplotype (~3.1 Gbp) | 6 | ~15 GB | ~18 GB |
| Diploid assembly, e.g. HG002 v1.1 (~6.0 Gbp) | 6 | ~29 GB | ~34 GB |

Measured on HG002 v1.1 against `HKS_human_CHM13_v2`: 21.7 GB of presmoothed BED, 7.3 GB of smoothed BED, and a 5.4 GB lookup TSV live at the peak. The BED content is the same for either backend, so the figures hold for both.

`--bgzip` (the default) does **not** lower the peak — every BED is written in full before the compression pass starts. What does lower it:

- `--feature-set` to annotate a subset at a time
- `--no-keep-presmoothed` (drops the larger of the two variants, ~0.6 of the 0.8 GB) or `--no-smooth`
- splitting a diploid assembly into per-haplotype runs

`annotate` estimates this footprint before starting and refuses if `--outdir` can't hold it, so a full disk is caught in a second rather than after the k-mer query. The estimate uses an exact base count from a sibling `.fai` when one exists, and otherwise scales the input file size. Pass `--no-resource-check` if the estimate is wrong for your input.

## See also

- [`karyoscope karyotype`](karyotype.md) — render karyotype SVGs (runs annotate under the hood)
- [`karyoscope bin`](bin.md) — aggregate base-pair BEDs into larger bins
