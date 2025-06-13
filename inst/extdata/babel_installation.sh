# 1. Install dependencies
sudo apt-get update
sudo apt-get install -y \
build-essential \
cmake \
libeigen3-dev \
libxml2-dev \
zlib1g-dev \
libboost-all-dev \
git

# 2. Download Open Babel 2.4.1
wget https://github.com/openbabel/openbabel/archive/refs/tags/openbabel-2-4-1.tar.gz
tar -xzf openbabel-2-4-1.tar.gz
cd openbabel-openbabel-2-4-1

# 3. Build and install it
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local/openbabel2 -DOPENBABEL_INSTALL_LIB_DIR=lib
make -j$(nproc)
sudo make install