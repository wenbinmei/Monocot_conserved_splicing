# This repo contains the main code used in paper **Evolutionarily Conserved Alternative Splicing Across Monocots** 
## You can find current manuscript on the [bioRxiv](http://www.biorxiv.org/content/early/2017/03/25/120469)
## This repo also contains the GTF file of alternative splicing isoiforms for each species generated in the paper
## If you are using the code or file in this repo, we are welcoming you to cite our paper [Mei et al. 2017](http://www.biorxiv.org/content/early/2017/03/25/120469)

##  This paper has two major conponents of the analyses: identify alternative splicing events and identify conserved alternative splicing across species
### first step in identify alternative splicing in each species
1. First of all you can use your favourite softwares to assembly the transcripts, either via denovo or reference-based assembly. In this analysis, in order to capture more transcripts diversity we utilized three different software to assemble the transcripts: cufflinks, Trinity and StringTie.
2. After you assemble the transcripts and you can filter them based on some criterials (details see our paper) and doing additional clustering via PASA to build the overall isoforms.
3. Next is the multiple-steps to clean the isoforms built from PASA, start to remove the loci built from PASA only have single transcript.
4. Use Spanki to calculate entropy score and to filter the isoforms based on splice junction support.
5. remove intron retention isoforms with low intron coverage support.
6. filtering the isoforms based on expression level and isoform ratio.

### second part is about identify the conserved alternative splicing across species, the Orthologroup is identified via OrthoFinder.
1. To extract the splice junction sequence around the alternative splicing events.
2. combine the splice junction of the same alternative splicing events (such as intron retention, exon skip) across species and pooled these splice junction.
3. tblastx the pool of splice junction tag with each other, and require to be considered conserved alternative splicing if both two sides of splice junction are conserved and also these two genes are in the same OrthoGroup identified by OrthoFinder, and also alternative splicing event type needs to be the same. Then, we can clustering these conserved alternative splicing evnets across different species if they are in the same OrthoGroup and share with each other.

## In the folder AS_GTF contains all the GTF files for alternative splicing isoforms identified in this study in each species

## License
This repo is free and open source for research usage, licensed under [GPLv2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html). 
