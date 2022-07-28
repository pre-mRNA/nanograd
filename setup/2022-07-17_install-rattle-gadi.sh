


# download
cd /g/data/lf10/as7425/apps/
git clone --recurse-submodules https://github.com/comprna/RATTLE
cd RATTLE

# install 
module load gcc/system
./build.sh

# rattle step 1: cluster reads
./rattle cluster -i reads.fq,reads2.fq,reads3.fq -l lib1,lib2,lib3 -t 48 --iso --rna -o ${out_dir}

# step 2: cluster extraction
