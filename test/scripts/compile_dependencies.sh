cd test
git clone https://github.com/lh3/bwa.git
cd bwa
git checkout tags/0.7.9a
make
cd ..
PATH=`pwd`/bwa:$PATH
git clone https://github.com/arq5x/bedtools2.git
cd bedtools2
git checkout tags/v2.25.0
make
cd ..
PATH=`pwd`/bedtools2/bin:$PATH
cd ..
echo $PATH