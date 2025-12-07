#/bin/bash -e

rm -rf LKH-3.0.13 LKH-3.0.13.tgz
wget http://webhotel4.ruc.dk/~keld/research/LKH-3/LKH-3.0.13.tgz
tar xvfz LKH-3.0.13.tgz
pushd LKH-3.0.13
make
cp LKH ${HOME}/.local/bin/
popd