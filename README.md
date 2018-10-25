# alex_loops_the_loop
An example loop over ATLAS' data format. Does some neat things too!

# Installation and Running

```bash
mkdir analysis_dir/; cd analysis_dir/
mkdir source/; cd source/
git clone git@github.com:dantrim/alex_loops_the_loop.git
lsetup "asetup AnalysisBase,21.2.45,here"
cd ..
mkdir build; cd build;
make -j
source x86*/setup.sh
run_loop -h
```

# Example File
Right now I have an example file on UCI brick: `/data/uclhc/uci/user/dantrim/susyNtProduction/n0303/samples/data16_13TeV.00303291.physics_Main.deriv.DAOD_SUSY2.r9264_p3083_p3372/DAOD_SUSY2.12594827._000052.pool.root.1`

# Run with file
run_loop -i <file>
  
# Output
Creates a file `loop_fun.root` with some histograms in them with random things being histogrammed.
