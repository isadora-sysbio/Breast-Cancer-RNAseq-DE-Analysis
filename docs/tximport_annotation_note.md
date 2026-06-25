# tximport Annotation Note

During Salmon-to-DESeq2 preparation, transcript-level Salmon quantifications were imported using tximport.

A transcript-to-gene mapping table (tx2gene) was generated from the Ensembl GRCh38.112 GTF annotation file.

A total of 13,448 transcript IDs present in the Salmon quant.sf files were not found in the tx2gene table.

This represented approximately 7% of the transcript IDs in the Salmon output. Around 93% of transcript IDs matched successfully.

The unmatched IDs were ENST transcript IDs. This likely reflects small differences between the transcript FASTA used to build the Salmon index and the GTF annotation used to generate tx2gene.

tximport completed successfully and summarized abundance, counts, and effective length values for the matched transcript-to-gene mappings.

This was treated as an annotation/reference compatibility note, not as a pipeline failure.
