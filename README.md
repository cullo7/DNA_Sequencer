# DNA_Sequencer

![alt text][four_branch]

###1. Description
Compiler that finds the melting temperature for dissociation of a DNA duplex with one or two consecutive mismatched base pairs.

###2. Parameters
Input can be a double-stranded DNA with a up to two mismatched base pairs, G-C pairs and and A-T pairs.

###3. Purpose
If we can find the melting temperature for any DNA strand, we can accurately denature and anneal any DNA strand within a very small margin of error.

###4. Formula
Our formula for calculating melting temperature is based on the Gibb's free energy, enthalpy, and entropy contribution of each base pair using the nearest neighbor calculation. With the base melting temperature, we then take into account numerous other factors such as:
* Salinity percentage
* G-C content
* Oligonucleotide molarity

Our formula is derived from the Van 't Hoff equation where we solve for the dissociation temperature.

Our Equation:

![alt text][Our_Equation]

Van 't Hoff Equation:

![alt text][Van_Hoff]
 
###5. Setup
DNA_Sequencer is configured for Python3. All the necessary modules can be downloaded with pip install 
and our requirements.txt file.
```
pip install -r requirements.txt
```
###6. Usage

####Running Program

```
$ make run
or
$ make irun (python3 main.py < input.txt)
or
$ python3 main.py
```

####Cleaning, Editing

```
$ make clean --cleans build and execution artifacts
$ make lint --check style with flake8
$ make isort --sort import statements
$ make help --information on makefile
```

####Running Multiple Sequences

```
[DNA_compiler]: multiple (or 'm')
...
1: R3 --random sequence of length 3
2: 3M3 --three matches, one mismatch, three matches
3: R5MR2 --five random matches, one mismatch, two random matches
4: quit
```

####Viewing Test Cases

```
[DNA_compiler]: test [1-15] (or 't')
...
energy:#
enthalpy:#
entropy:#
temperature:#
length:#
```

####Running Test Cases

```
[DNA_compiler]: show [1-15] (or 's')
T A
G C
C G
...
Expected:
energy:#
enthalpy:#
entropy:#
temperature:#
length:#
```

###7. References
1. "Nucleic acid thermodynamics." Wikipedia. Wikimedia Foundation, n.d. Web. 31 Dec. 2016.

2. Allawi, Hatim T., and John Santalucia. "Thermodynamics and NMR of Internal G·T Mismatches in DNA." Biochemistry 36.34 (1997): 10581-0594. Web.

3. Allawi, H. "Thermodynamics of internal C.T mismatches in DNA." Nucleic Acids Research 26.11 (1998): 2694-701. Web.

4. Santalucia, John, Hatim T. Allawi, and P. Ananda Seneviratne. "Improved Nearest-Neighbor Parameters for Predicting DNA Duplex Stability †." Biochemistry 35.11 (1996): 3555-562. Web.

5. "Calculating Tm of Oligonucleotide Duplexes - IDT Biophysics." N.p., n.d. Web. 31 Dec. 2016.

6. Allawi, Hatim T., and John Santalucia. "Nearest-Neighbor Thermodynamics of Internal A·C Mismatches in DNA:  Sequence Dependence and pH Effects." Biochemistry 37.26 (1998): 9435-444. Web.

7. Peyret, Nicolas, P. Ananda Seneviratne, Hatim T. Allawi, and John Santalucia. " Nearest-Neighbor Thermodynamics and NMR of DNA Sequences with Internal A‚A, C‚C, G‚G, and T‚T Mismatches†." Biochemistry 38 (1998): 3468-477. Web. 31 Dec. 2016.

8. Allawi, Hatim T., and John Santalucia. "Nearest Neighbor Thermodynamic Parameters for Internal G·A Mismatches in DNA." Biochemistry 37.8 (1998): 2170-179. Web.

9. Lomzov, Alexander A., Yury N. Vorobjev, and Dmitrii V. Pyshnyi. "Evaluation of the Gibbs Free Energy Changes and Melting Temperatures of DNA/DNA Duplexes Using Hybridization Enthalpy Calculated by Molecular Dynamics Simulation." The Journal of Physical Chemistry B 119.49 (2015): 15221-5234. Web.

10. Santalucia, J. "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." Proceedings of the National Academy of Sciences 95.4 (1998): 1460-465. Web.

[four_branch]:https://upload.wikimedia.org/wikipedia/commons/thumb/9/92/Holliday_junction_coloured.png/400px-Holliday_junction_coloured.png
[Van_Hoff]:https://encrypted-tbn3.gstatic.com/images?q=tbn:ANd9GcRpOoWiRMsxcG6t11c1Tj4NLbHZQT_n9cu-h8c0Q9vV58TUh7Xy0Q
[Our_Equation]:http://biotools.nubic.northwestern.edu/images/thermoeq5.gif
