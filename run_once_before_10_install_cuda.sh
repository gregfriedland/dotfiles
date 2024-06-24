#!/bin/bash

if test -f /opt/google/cuda-installer
then
  exit
fi

sudo mkdir -p /opt/google/cuda-installer
cd /opt/google/cuda-installer/ || exit

sudo curl -fSsL -O https://github.com/GoogleCloudPlatform/compute-gpu-installation/releases/download/cuda-installer-v1.1.0/cuda_installer.pyz
sudo python3 cuda_installer.pyz install_cuda
