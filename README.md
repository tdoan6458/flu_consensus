# flu_consensus

## Building Reference Database

1. Download latest Genbank NCBI Influenza virus sequences.
2. Filtering by minimum length requirement for each segment.
3. Filter out sequences that can’t be assigned segment.
4. Assign segment# to each sequence using keyword search.

```
awk '
BEGIN{
  # minimum length by segment
  minlen[1]=2200; minlen[2]=2200; minlen[3]=2200;
  minlen[4]=1700; minlen[5]=1700;
  minlen[6]=1300;
  minlen[7]=900;
  minlen[8]=800;

  processed=0; unknown=0; tooshort=0;

  OFS="\t";
}

function seg_from_header(h,    t){
  t = tolower(h);

  # --- Gene-name based mapping (case-insensitive) ---
  if (t ~ /(^|[^a-z0-9])pb2([^a-z0-9]|$)/) return 1;
  if (t ~ /(^|[^a-z0-9])pb1([^a-z0-9]|$)/) return 2;

  # PA: try to avoid matching random "pa" inside other words
  if (t ~ /(^|[^a-z0-9])pa([^a-z0-9]|$)/ || t ~ /polymerase[[:space:]]+acidic/) return 3;

  # HA: HA1, hemagglutinin, or standalone HA
  if (t ~ /ha1/ || t ~ /hemagglutinin/ || t ~ /(^|[^a-z0-9])ha([^a-z0-9]|$)/) return 4;

  if (t ~ /(^|[^a-z0-9])np([^a-z0-9]|$)/) return 5;

  # NA / NB / neuraminidase
  if (t ~ /(^|[^a-z0-9])na([^a-z0-9]|$)/ || t ~ /(^|[^a-z0-9])nb([^a-z0-9]|$)/ || t ~ /neuraminidase/) return 6;

  # M / BM2 (prefer bm2/matrix terms; standalone m as fallback)
  if (t ~ /bm2/ || t ~ /matrix[[:space:]]+protein/ || t ~ /(^|[^a-z0-9])m([^a-z0-9]|$)/) return 7;

  # NS / NS1 / NS2 / nonstructural
  if (t ~ /ns1/ || t ~ /ns2/ || t ~ /nonstructural/ || t ~ /(^|[^a-z0-9])ns([^a-z0-9]|$)/) return 8;

  # --- Fallback: "segment: N" parsing ---
  if (match(t, /segment[: ]+[0-9]+/)) {
    s = substr(t, RSTART, RLENGTH);
    gsub(/[^0-9]/, "", s);
    if (s+0 >= 1 && s+0 <= 8) return (s+0);
  }

  return 0; # unknown
}

function flush_record(    seg, need){
  if (header == "") return;

  seg = seg_from_header(header);

  if (seg == 0) {
    print header > "unprocessed_headers.txt";
    unknown++;
    return;
  }

  need = minlen[seg];
  if (seqlen < need) {
    print header "\tLEN=" seqlen "\tSEG=" seg "\tMIN=" need > "too_short_headers.txt";
    tooshort++;
    return;
  }

  # passed
  print acc, seg;
  print header > "processed_headers.txt";
  processed++;
}

# FASTA parsing
/^>/{
  flush_record();

  header = $0;
  seqlen = 0;

  # accession = first token after ">"
  if (match($0, /^>[^ ]+/, m)) acc = substr(m[0], 2);
  else acc = "";
  next;
}

{
  gsub(/[[:space:]]/, "", $0);
  seqlen += length($0);
}

END{
  flush_record();

  print "Processed records (passed length): " processed > "/dev/stderr";
  print "Unknown segment (see unprocessed_headers.txt): " unknown > "/dev/stderr";
  print "Too short (see too_short_headers.txt): " tooshort > "/dev/stderr";
}
' Influenza_all_ref.fasta > acc2segment.tsv
```
In this example, I downloaded all Influenza B sequences from [NCBI Genbank](https://www.ncbi.nlm.nih.gov/nuccore/?term=txid11520[organism:exp]%20AND%20biomol_genomic[prop]). Sequences that are too short (usually are partial cds) are filtered out and printed to [too_short_headers.txt]. Sequences that are assigned with segment# are printed to [processed_headers.txt]. The [acc2segment.tsv] have two columns, first column being the accessionID and second column being segment# (from 1 to 8) that the each sequence are assigned to.
> \[!NOTE]
>
> Here I was lenient on setting the minimum length requirements for each segments. Change them if you need to. Check the unprocessed_headers.txt and acc2segment.tsv to see if you could add or remove keywords for assigning sequences to their segments. You might want to blast the sequence to NCBI to check if there are any mismatches.
```
awk '
NR==FNR {
    ids[$0];          # store full header line (including >)
    next
}
/^>/ {
    keep = ($0 in ids)   # check if header exactly matches
}
keep
' processed_headers.txt InfluenzaB_all_ref.fasta > InfluenzaB_filtered_ref.fasta
```

## Align samples to reference database

### 1. Align sample reads to assigned, filtered reference database
```
bwa mem -t 12 \
    InfluenzaB_filtered_ref.fasta \
    ${library}_1.fastq.gz \
    ${library}_2.fastq.gz | \
    samtools view -bS -F 4 - | \
    samtools sort -o /bwa_output/${library}.sorted.bam -

samtools index /bwa_output/${library}.sorted.bam
```
### 2. Compute max meandepth by each segment
```
samtools coverage \
    -o /coverage/${library}_coverage.txt \
    /bwa_output/${library}.sorted.bam

awk 'FNR==NR {seg[$1]=$2; next} {print $0, seg[$1]}' acc2segment.tsv /coverage/${library}_coverage.txt > ${library}_coverage_with_segment.txt

sed -i '1s/$/	segment/' coverage/${library}_coverage_with_segment.txt
```
### 3. Build individual reference fasta file for each sample
```
# 1) Get the IDs from coverage file (skip header, take first column)
tail -n +2 coverage/${library}_highest_meandepth_by_segment.txt | awk '{print $1}' > /coverage/${library}_best_aligned_ids.txt

# 2) Extract matching sequences from FASTA
seqkit grep -f /coverage/${library}_best_aligned_ids.txt InfluenzaB_filtered_ref.fasta > /references/${library}_ref.fasta
```
### 4. Align each sample to its own new reference and generate consensus
```
bwa index /references/${library}_ref.fasta
bwa mem -t 12 \
    /references/${library}_ref.fasta \
    ${library}_1.fastq.gz \
    ${library}_2.fastq.gz | \
    samtools view -bS -F 4 - | \
    samtools sort -o /bwa_output2/${library}.sorted.bam -

samtools index /bwa_output2/${library}.sorted.bam
```
### 6. Assess coverage and depth
```
samtools coverage \
    -o coverage2/${library}_coverage.txt \
    bwa_output2/${library}.sorted.bam

awk 'FNR==NR {seg[$1]=$2; next} {print $0, seg[$1]}' acc2segment.tsv coverage2/${library}_coverage.txt > coverage2/${library}_coverage_with_segment.txt

sed -i '1s/$/	segment/' coverage2/${library}_coverage_with_segment.txt

samtools depth \
    -a \
    -o coverage2/${library}_depth.txt \
    bwa_output2/${library}.sorted.bam

awk 'FNR==NR {seg[$1]=$2; next} {print $0, seg[$1]}' acc2segment.tsv coverage2/${library}_depth.txt > coverage2/${library}_depth_with_segment.txt
```
### 7. Generate consensus sequence
```
samtools consensus /bwa_output2/${library}.sorted.bam > /consensus/${library}_consensus.fasta
```
