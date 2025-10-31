Analysis of a multi-sample FASTA file containing the entire mitochondrial genome called from WGS data (see 04-variantcalling directory). The desired output is a phylogenetic tree and a haplotype network.

1. The FASTA files were created with the command, in order to:
(a) make the absence of data visible (noted with _, but our dataset with WGS only had the first base missing in some samples in the mitogenome).
(b) annotate deletions with -, which is comprehended as such by MEGA
(c) annotate and remove insertions, which are messing up downstream analyses

	# Generate consensus FASTA for variant sites only
 bcftools consensus \
  -f $REF \
  -a '_' \ #absence of data
  --mark-del '-' \ #deletions
  --mark-ins lc \ #insertions
  "$vcf" | \
 tr -d '[a-z]' | \
 sed "1s/.*/>$sample_name/" \
 > "$OUT_DIR/${sample_name}_indel.fasta"
done

2. Alignment was done manually in MEGA12 and exported as nexus. Added a traits block for the creation of population-informative haplotype network in PopART.
   
 # Example snipset from the traits block in nexus
Begin Traits;
Dimensions NTraits=20;
Format labels=yes missing=? separator=Comma;
TraitLatitude 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
TraitLongitude 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
TraitLabels A1_AnOudin A2_BuSN A3_Mech A4_Pl15_YEL A5_GenM B1_BKN1 B2_OM2 B3_ZW B4_OHZ B5_DA2 C1_BlfN C2_MO C3_Ter1 C4_BW_48630 C5_BW_36962 D1_CBOO6 D2_LRV D3_BKLE5 D4_BW_62256 D5_BW_22050;
Matrix
A1_AnOudin_C02_combined 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
[...]
END
;

3. The haplotype network was constructed in popart-1.7 under minimum spanning / median joining, and was manually edited for colouration and position of the nodes.
