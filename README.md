# DNA_compiler

![alt text][four_branch]

###1. Description
Compiler that finds the melting temperature for dissociation of a DNA duplex.

###2. Parameters
Input can be a double-stranded DNA with a variable number of base pair mismatches as well as any amount of salt content, G-C pairs and other determinants of melting temperature (section 4)

###3. Purpose
If we can find the melting temperature for any DNA strand, we can accurately dissociate and associate any DNA strand within a very small margin of error.

###4. Formula
Our formula(Fig 1) for melting temperature is based on the Gibb's free energy of each base pair using the nearest neighbor calculation. With the base melting temperature, we then take into account numerous other factors such as:
* Salinity percentage
* G-C content
* 

###5. References
1. Seeman, Nadrian C. Structural DNA nanotechnology. Cambridge, United Kingdom: Cambridge U Press is part of the U of Cambridge, 2015. Print.
2. Lomzov, Alexander A., Yury N. Vorobjev, and Dmitrii V. Pyshnyi. "Evaluation of the Gibbs Free Energy Changes and Melting Temperatures of DNA/DNA Duplexes Using Hybridization Enthalpy Calculated by Molecular Dynamics Simulation." The Journal of Physical Chemistry B119.49 (2015): 15221-5234. Web.
3. "Nucleic acid thermodynamics." Wikipedia. Wikimedia Foundation, n.d. Web. 23 Dec. 2016.
4. Santalucia, J. "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." Proceedings of the National Academy of Sciences 95.4 (1998): 1460-465. Web.
5. Santalucia, John, Hatim T. Allawi, and P. Ananda Seneviratne. "Improved Nearest-Neighbor Parameters for Predicting DNA Duplex Stability †." Biochemistry 35.11 (1996): 3555-562. Web.
6. "Calculation of Tm for Oligonucleotide Duplexes ..." Www.idtdna.com. N.p., n.d. Web. 23 Dec. 2016.

[four_branch]:https://upload.wikimedia.org/wikipedia/commons/thumb/9/92/Holliday_junction_coloured.png/400px-Holliday_junction_coloured.png
