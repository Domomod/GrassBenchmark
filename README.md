Example of setting up benchmark data and running GrassSV on `fastadna_all` dataset

```
sudo apt-get install git-lfs
git clone https://github.com/Domomod/GrassBenchmark.git
git clone https://github.com/Domomod/GrassSV.git
git clone https://github.com/swacisko/ALGA.git
PATH=$PATH:<path to GrassBenchmark/src>:<path to GrassSV/>:<path to ALGA/build/>
cd GrassBenchmark
git lfs fetch --all
cd fastadna_all
run_grass.sh -g ../ref/ref.fa -r reads1.fq -R reads2.fq -l 200 --allow-overwrite
```
