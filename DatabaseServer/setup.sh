#! /bin/bash
# N.B. RUN WITH SUDO PRIVILEGES

len=${echo $1 | wc -c}
if [ -z $1 || $1 = ^[0-9]+$ || len -ne 2 ]; then
    echo "A two digit data fragment larger than zero must be specified!"
    echo "Ex: sudo ./setup.sh 01"
fi
# Gets necessary packages
yum install -y git, wget

pip3 install -r "requirements.txt"

# Sets up blast directories
cd $HOME 
sudo mkdir bin blastdb queries fasta results blastdb_custom

# Pulls data fragment
cd blastdb
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt.00.tar.gz
wget "https://ftp.ncbi.nlm.nih.gov/blast/db/nt.{$1}.tar.gz"
tar -xvf *.tar.gz
rm *.tar.gz

echo "Installation complete"