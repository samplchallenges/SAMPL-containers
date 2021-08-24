# Description
Details how to setup and build Autodock Vina Base docker container for Autodock Vina Docking container to inherit from. 

# Setup:
1. `mkdir SAMPL-league/attic/adv-base/dependencies`
2. Download Autodock Tools linux x86 .tgz file (`autodock_vina_1_1_2_linux_x86.tgz`) from http://vina.scripps.edu/download.html
3. `mv autodock_vina_1_1_2_linux_x86.tgz SAMPL-league/attic/adv-base/dependencies`
4. `cd SAMPL-league/attic/adv-base/dependencies`
5. `tar -xvf autodock_vina_1_1_2_linux_x86.tgz`
6. `rm autodock_vina_1_1_2_linux_x86.tgz`
7. `mv autodock_vina_1_1_2_linux_x86 adv`
8. Download MGL Tools linux x86 .tar.gz (`mgltools_x86_64Linux2_1.5.6.tar.gz`) from http://mgltools.scripps.edu/downloads
9. `mv mgltools_x86_64Linux2_1.5.6.tar.gz SAMPL-league/attic/adv-base/dependencies`
10. `cd SAMPL-league/attic/adv-base/dependencies`
11. `tar -xvf mgltools_x86_64Linux2_1.5.6.tar.gz`
12. `rm mgltools_x86_64Linux2_1.5.6.tar.gz`
13. `mv mgltools_x86_64Linux2_1.5.6 mgl`
14. `vi mgl/install.sh`
15. Change line 6 `TarDir=` to `TarDir="/opt/app/dependencies/mgl/"`
16. Change line 7 `export MGL_ROOT=""` to `export MGL_ROOT="/opt/app/dependencies/mgl/"`
17. Save install.sh

# Build
1. `cd adv-base/`
2. `docker build -t adv-base .`
