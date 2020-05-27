# QuantumCircuits

Here we implement some arithmetic circuits required by the floating point operations featured in https://arxiv.org/pdf/1807.02023.pdf.
The addition circuit from https://arxiv.org/pdf/0910.2530.pdf is also implemented.

All implementations make use of ProjectQ: An open source software framework for quantum computing https://projectq.ch.

I recommend using the a conda environment when using Project Q and download using

``` 
conda install -c rpmuller projectq
```

To see an example of running the addition circuit, run the following:

```
python main.py
```