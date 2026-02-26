# flu_consensus

## Building Reference Database

1. Download latest Genbank NCBI Influenza virus sequences.
3. Assign segment# to each sequence using keyword search.
4. Filtering by minimum length requirement for each segment.
5. Filter out sequences that can’t be assigned segment.

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
' InfluenzaB_all_ref.fasta > acc2segment.tsv

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

## Align samples to reference

1. Align sample reads to assigned, filtered reference database.
2. Compute max meandepth by each segment.
3. Build individual reference fasta file for each sample.
4. Align each sample to its own new reference.
5. Assess coverage and depth
6. Generate consensus sequence
