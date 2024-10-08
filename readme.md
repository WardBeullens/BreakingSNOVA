# Demonstration of weak key attack against SNOVA

This repo has a demonstration of a weak-key forgery attack against SNOVA with parameter set (v,o,q,l) = (37, 17, 16, 2)

The attack runs in a few minutes and works against one out of roughly every 140000 public keys. 

This repo contains a weak key pk.txt on which the attack is performed. You can run the attack on your own weak key by replacing the key in pk.txt. 
You can generate your own weak keys with our weak key finding script. 

This project builds on the implementation of SNOVA by the SNOVA team. The SNOVA readme can be found below.

## Running the attack

Before running the attack make sure m4gb is installed in the m4gb subfolder.

If necessary, run 

```
git submodule update --init
```

to download m4gb. Follow the instructions in the 4mgb project to install m4gb. 
If m4gb is installed, you can run the sage script: 

```
sage Attack_demo.sage
```

## Searching for weak keys 

You can compile the key search program using the makefile: 

```
make find_vulnerable_key
```

Then run the program with two integers as arguments. Running the program with arguments X and Y generates Y public keys and checks if they are weak. The keys are generated deterministically as in the NIST script for generating KATs. 
The first X keys are skipped (This can be useful if you want to parallelize the search).
The first weak public key that is found is printed on the screen. You can copy-paste it into the pk.txt file to run the attack on the weak key you found.

For example, running

```
./find_vulnerable_key 0 1000000
``` 

will print a weak public key, as well as: 

```
first key found after 734846 attempts
tries: 1000000, tries/second = 2941.167821
lowest ranks:
 0: 0, 0.000000
 1: 4, 0.000004
 2: 20775, 0.020775
 3: 979221, 0.979221
 4: 0, 0.000000
```

Indicating that 4 out of the first 1.000.000 keys is vulnerable to the weak key attack with rank 1. 

README from SNOVA implementation
=======
This repository contains the latest official reference & AVX implementation of the SNOVA signature scheme in C language.

Please refer to this [document](https://github.com/PQCLAB-SNOVA/SNOVA/blob/main/doc/NIST_Submission_Source_Code_Differences_Document.md) for the main differences between the previously submitted code to NIST and the current version.

SNOVA parameter
-------
| SL |         Name  |  V |  O |  L | sk size (esk) | sk size (ssk) |    pk size   | sign size  |
|----| --------------|----|----|----|---------------|---------------|--------------|------------|
|  1 | SNOVA_37_17_2 | 37 | 17 |  2 |    60008(+48) |            48 |    9826(+16) |   108(+16) |
|  1 |  SNOVA_25_8_3 | 25 |  8 |  3 |    37962(+48) |            48 |    2304(+16) | 148.5(+16) |
|  1 |  SNOVA_24_5_4 | 24 |  5 |  4 |    34112(+48) |            48 |    1000(+16) |   232(+16) |
|  3 | SNOVA_56_25_2 | 56 | 25 |  2 |   202132(+48) |            48 |   31250(+16) |   162(+16) |
|  3 | SNOVA_49_11_3 | 49 | 11 |  3 |   174798(+48) |            48 |  5989.5(+16) |   270(+16) |
|  3 |  SNOVA_37_8_4 | 37 |  8 |  4 |   128384(+48) |            48 |    4096(+16) |   360(+16) |
|  5 | SNOVA_75_33_2 | 75 | 33 |  2 |   515360(+48) |            48 |   71874(+16) |   216(+16) |
|  5 | SNOVA_66_15_3 | 66 | 15 |  3 |   432297(+48) |            48 | 15187.5(+16) | 364.5(+16) |
|  5 | SNOVA_60_10_4 | 60 | 10 |  4 |   389312(+48) |            48 |    8000(+16) |   560(+16) |

Build instructions
-------
These implementations include multiple test programs and an easy-to-compile Makefile.

Prerequisites for compiling this repository
-------
None.
This repository contains the symmetric primitive libraries as source code. See the respective files for credits and copyright notices.

Test programs
-------
The steps to compile and run the SNOVA signature scheme's test program on Linux are as follows:

1. test keypair, sign and verify
```bash
make clean
make test
```
2. test nist api
```bash
make clean
make test_api
```
3. test speed
```bash
make clean
make test_speed
```
```bash
make clean
make test_speed  TEST_SPEED_N=4096
```

Set the parameters
-------
Only need to modify the SNOVA_V, SNOVA_O, SNOVA_L, SK_IS_SEED in the Makefile file.

Example: (makefile line 37)
```make
SNOVA_V ?= 24
SNOVA_O ?= 5
SNOVA_L ?= 4

SK_IS_SEED ?= 0 # 0: sk = ssk; 1: sk = esk
```

Genarate KAT
-------
```bash
make clean
make PQCgenKAT
```

Tip: The following parameters can all be input during MAKE execution, such as
-------
```bash
make test_speed SNOVA_V=24 SNOVA_O=5 SNOVA_L=4 TEST_SPEED_N=4096
```

The following are the latest reference parameters.
-------

SL 1: 
| SNOVA_V | SNOVA_O | SNOVA_L |
|---------|---------|---------|
|      37 |      17 |       2 |
|      25 |       8 |       3 |
|      24 |       5 |       4 |

SL 3: 
| SNOVA_V | SNOVA_O | SNOVA_L |
|---------|---------|---------|
|      56 |      25 |       2 |
|      49 |      11 |       3 |
|      37 |       8 |       4 |

SL 5: 
| SNOVA_V | SNOVA_O | SNOVA_L |
|---------|---------|---------|
|      75 |      33 |       2 |
|      66 |      15 |       3 |
|      60 |      10 |       4 |


Example:
SL 5 (60, 10, 4)
```
make test_speed SNOVA_V=60 SNOVA_O=10 SNOVA_L=4 TEST_SPEED_N=512
```
```
make PQCgenKAT SNOVA_V=60 SNOVA_O=10 SNOVA_L=4 SK_IS_SEED=1
```
```
make test_api SNOVA_V=60 SNOVA_O=10 SNOVA_L=4 SK_IS_SEED=1
```


## Team Members

Thank you to our team members from SNOVA:

- Lih-Chung Wang
- Chun-Yen Chou
- Jintai Ding
- Yen-Liang Kuan
- Jan Adriaan Leegwater
- Ming-Siou Li
- Bo-Shu Tseng
- Po-En Tseng
- Chia-Chun Wang

## Team Contribution Disclaimer
All commits in this project are the result of team collaboration. The contributions made by the **pqclab-zero** account represent collective efforts of the team, **not individual contributions**.
