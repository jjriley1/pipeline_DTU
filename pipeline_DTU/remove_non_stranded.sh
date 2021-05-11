gzip -d export/agg-agg-agg.gtf.gz && grep -v -P '\t\.\t\.\t' export/agg-agg-agg.gtf > export/agg-agg-agg-strandedOnly.gtf && gzip export/agg-agg-agg.gtf
