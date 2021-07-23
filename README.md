# IsoAligner

Abstract
Summary: 
Aligning protein sequences of isoforms is necessary in order to transfer annotated mutational data from one splice variant to another. However, most alignment tools available online are fitted for homologue sequences, while simple mapping of amino acids across isoforms comes short. Here, we present a interactive alignment webservice which is tailored for the positional mapping of exons. The tool uses a customized Needleman-Wunsch algorithm including an open gap penalty combined with a gene-specific minimal exon length function and dynamically adjustable parameters. As an input, either gene/transcript/protein ID’s of various sequence databases (Ensembl, Uniprot, Refseq) or raw amino acid sequences, can be provided. The output of IsoAligner consists of pairwise alignments and a table of mapped amino acid positions between the canonical or supplied isoform ID and all alternative isoforms. IsoAligner’s human isoform library comprises in total of +1.3 million IDs mapped on 120k protein sequences. IsoAligner, is a fast and interactive alignment tool helping scientists to retrieve amino acids positions between different protein isoforms.

Availability: This tool is available as a Webservice on www.IsoAligner.org. A REST-API is available for programmatic access. The source code for both services can be found on https://github.com/JacobHanimann/IsoAligner

https://share.streamlit.io/jacobhanimann/isoaligner/main/Webinterface.py

