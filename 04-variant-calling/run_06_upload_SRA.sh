#!/bin/bash -l
 
# -------------------------------------------------------------------------------
#   upload FASTQ files on SRA 
# -------------------------------------------------------------------------------

#SBATCH --cluster=wice
#SBATCH --job-name=upload_SRA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --output=logs/upload_SRA.%j.log
#SBATCH -A lp_svbelleghem

# Environment setup
module load cluster/wice/batch

cd /lustre1/scratch/363/vsc36396/all_fqfiles

# FTP variables
FTP_HOST="ftp-private.ncbi.nlm.nih.gov"
FTP_USER="subftp"
FTP_PASS="B1ckyadixdaycsuD"
REMOTE_DIR="uploads/alexandros.kourtidis_kuleuven.be_mp7fUE8X/SUB16419844_raw_fq"

# Automate standard FTP
ftp -inv $FTP_HOST <<EOF
user $FTP_USER $FTP_PASS
cd $REMOTE_DIR
binary
mput *.fq.gz
quit
EOF

echo "done"
