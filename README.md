# Fast Community Search

This repository contains the implementation of fast community search algorithm presented in the [paper](https://dl.acm.org/doi/pdf/10.1145/3200863) "Searching for a Single Community in a Graph" by Ray et al. (2018). 

Please refer our paper if you find this code/algorithm useful in your research.

```
@article{ray2018searching,
  title={Searching for a Single Community in a Graph},
  author={Ray, Avik and Sanghavi, Sujay and Shakkottai, Sanjay},
  journal={ACM Transactions on Modeling and Performance Evaluation of Computing Systems (TOMPECS)},
  volume={3},
  number={3},
  pages={1--17},
  year={2018},
  publisher={ACM New York, NY, USA}
}
``` 

## Environment

The code was implemented and tested using `Matlab (R2016a)`.

## Basic Usage

```
% define community membership file
CommMembershipFile = 'NOComm_Neq_n1000_K5_v4';

% define edge probabilities
p = 0.1;
q = 0.01;

% define side information, number of labeled nodes/community
numLabel = 10;

% perform community search
communitySearchAll(CommMembershipFile, p, q, numLabel);

% example output (should have average error ~ 1%)
Generating graph ...
Running community search algo for k = 1
Running community search algo for k = 2
Running community search algo for k = 3
Running community search algo for k = 4
Running community search algo for k = 5
Percentage error = 0.2%
Experiment complete !
```
