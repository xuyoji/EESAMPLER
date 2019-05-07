#BUILDING
Linux environment is Required.
Install dependencies
```
sudo apt install git g++ make python-minimal
```

Clone repos
```
git clone https://github.com/xuyoji/EESAMPLER.git
git clone https://github.com/Z3Prover/z3.git
```

Build z3
```
cd z3
python scripts/mk_make.py
cd build
make
sudo make install
cd ../..
```

Build EESampler

```
cd EESAMPLER
make
```

# Running
'''
./EESampler -n 1000000 -t 7200.0 -i -r formula.cnf
'''
Use '-i' to enable the independent support information.
Use '-n' to set the maximum samples number, default value is 10000000.
Use '-t' to set the time limit, default value is 7200.
Use '-r' to enable the heuristic soft constraint choose method, which will random choose the soft constraint variable by the probility of the repeat rate of samples.(beta function)
Use '-d' to enable debug mode, which will print debug information on the terminal.

#Benchmarks
We use Benchmarks from UniGen. You can find them in 'unigen-benchmarks' folder which 370 DIMACS files included.



