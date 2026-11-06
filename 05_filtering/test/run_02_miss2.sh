#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=filter_test
#SBATCH --output=logs/filtertest.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem


# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load VCFtools/0.1.16-GCC-10.3.0

# ── Directories & Files ─────────────────────────────────────────────────────────
cd /lustre1/scratch/XXX/vscXXXXX/filtering

IN_DIR="/lustre1/scratch/XXX/vscXXXXX/variantcalling/VCFs"
VCFDIR="test/vcfs"
SUMDIR="test/summaries"
mkdir -p $VCFDIR $SUMDIR

base="batch1234_scaffold_10"
IN_VCF="${IN_DIR}/${base}.vcf.gz"
VCF1="${VCFDIR}/${base}.DPonly_noMISS.vcf.gz"
VCF2="${VCFDIR}/${base}.snpQUALDPMQ_DPGQ_filtered.vcf.gz"
VCF3="${VCFDIR}/${base}.snpQUALDPMQ_filtered.vcf.gz"

# Creating imiss files of the  vcf files
#echo "creating imiss files for $IN_VCF ..."
#vcftools --gzvcf $IN_VCF --missing-site --out "${IN_VCF}_miss"
#vcftools --gzvcf $IN_VCF --missing-indv --out "${IN_VCF}_miss"

echo "creating imiss files for $VCF1 ..."
#vcftools --gzvcf $VCF1 --missing-site --out "${VCF1}_miss"
vcftools --gzvcf $VCF1 --missing-indv --out "${VCF1}_miss"

#echo "creating imiss files for $VCF2 ..."
#vcftools --gzvcf $VCF2 --missing-site --out "${VCF2}_miss"
#vcftools --gzvcf $VCF2 --missing-indv --out "${VCF2}_miss"

#echo "creating imiss files for $VCF3 ..."
#vcftools --gzvcf $VCF3 --missing-site --out "${VCF3}_miss"
#vcftools --gzvcf $VCF3 --missing-indv --out "${VCF3}_miss"

echo "done"