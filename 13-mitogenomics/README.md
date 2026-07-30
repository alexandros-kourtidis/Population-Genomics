Analysis of a multi-sample FASTA file containing the entire mitochondrial genome called from WGS data (see 04-variantcalling directory). The desired output is a haplotype network and/or a simple phylogenetic tree.

# VCF --> FASTA --> NEXUS --> Haplotype networks

**1**. Aligned single-sample FASTA files were extracted from each singular vcf using the **bcftools consensus** command and then concantenated in a multi-sample FASTA file. The bcftools consensus specifications were:

  (a) make the absence of data visible (noted with _, but our dataset with WGS only had the first base missing in some samples in the mitogenome).

  (b) annotate deletions with -, which is comprehended as such by MEGA

  (c) annotate and remove insertions, which are messing up downstream analyses

  Code snipset:
  bcftools consensus \
    -f $REF \
    -a '_' \ #absence of data
    --mark-del '-' \ #deletions
    --mark-ins lc \ #insertions
    "$vcf" | \
   tr -d '[a-z]' | \
   sed "1s/.*/>$sample_name/" > "$OUT_DIR/${sample_name}_indel.fasta"
  done

**2**. Curation was done manually in **MEGA12** and exported as FASTA. Manually removed ambiguous sites (like Y, R, W, X, ...).
   
**3**. Imported in **DnaSP**, calculated haplotypes and exported as nexus. Manually added a traits block for the creation of population-informative haplotype network in PopART. An example of the file to be imported to popart is included here.

**4**. The haplotype network was constructed in **popart-1.7** under minimum spanning / median joining, and was manually edited for colouration and position of the nodes. This worklfow makes it possible to visualise the haplotypes in map version as well, with a piechart for each population.
