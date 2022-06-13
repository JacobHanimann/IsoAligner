# IsoAligner

Aligning protein isoform sequences is often performed in cancer diagnostics to homogenise mutation
annotations from different diagnostic assays. However, most alignment tools are fitted for homologous
sequences, leading often to alignments of non-identical exonic regions. Here, we present the interactive
alignment webservice IsoAligner for exact mapping of exonic protein subsequences. The tool uses a
customized Needleman-Wunsch algorithm including an open gap penalty combined with a gene-specific
minimal exon length function and dynamically adjustable parameters. As an input, IsoAligner accepts
either various gene/transcript/protein IDs from different databases (Ensembl, UniProt, RefSeq) or
raw amino acid sequences. The output of IsoAligner consists of pairwise alignments and a table
of mapped amino acid positions between the canonical or supplied isoform IDs and all alternative
isoforms. IsoAligner’s human isoform library comprises of over 1.3 million IDs mapped on over
120,000 protein sequences. IsoAligner, is a fast and interactive alignment tool for retrieving amino
acids positions between different protein isoforms. Its application will allow diagnostic and precision
medicine labs to detect inconsistent variant annotations between different assays and databases.
Availability: This tool is available as a Webservice on www.isoaligner.org. A REST API is available
for programmatic access. The source code for both services can be found on GitHub.
