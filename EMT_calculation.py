import os
import time
import glob
import numpy as np
from math import log
from math import sqrt
from math import exp
import cmath
from scipy.interpolate import InterpolatedUnivariateSpline


frequency_or_wavelength_m = []
frequency_or_wavelength_in = []
n_m = []
k_m = []
n_in = []
k_in = []
constant_n = []
constant_k = []


print ("")
print ("|______________________________________________________________________________|")
print ("|                                                                              |")
print ("|                            Effective medium thoery                           |")
print ("|______________________________________________________________________________|")
print ("|                                                                              |")
print ("| The effective-medium approximation is a method of treating a macroscopically |")
print ("| inhomogeneous medium,i.e., a medium in which quantities such as              |")
print ("| the conductivity σ, dielectric function ϵ,or elastic modulus vary in space.  |")
print ("|                                                                              |")
print ("|______________________________________________________________________________|")
print ("")
print ("")
print ("|--------------------------- Press enter to continue --------------------------|")
print ("")
input ("")

print (" Several methods to predict the mixing material's effective permittivity.")
print (" If you want to use the method by Maxwell Garnett model -> please press <1>")
print (" If you want to use the method by Bruggeman model -> please press <2>")
print (" If you want to use the method by Lichtenecker model -> please press <3>")
print (" If you want to use the method by Looyengat model -> please press <4>")
print (" If you wanto use the method by Monecke model -> please press <5>")
print (" If you wanto use the method by Hollow sphere equivalent model-> please press <6>")
print (" If you want to see the references, -> please press <0>")



answer = input ("\n")

if answer == '0' or answer == '1' or answer == '2' or answer == '3' or answer == '4' or answer == '5'  or answer == '6':
    answer = int(answer)

    if answer == 1:
        print ("-> You use the EMT method by Maxwell Garnett model")
        print ("")
    elif answer == 2:
        print ("-> You use the EMT method by Bruggeman model")
        print ("")
    elif answer == 3:
        print ("-> You use the EMT method by Lichtenecker model")
        print ("")
    elif answer == 4:
        print ("-> You use the EMT method by Looyengat model")
        print ("")
    elif answer == 5:
        print ("-> You use the EMT method by Monecke model")
        print ("")
    elif answer == 6:
        print ("-> You use the EMT method by Hollow sphere equivalen model")
        print ("")


while answer == 0:
    print (" References:")
    print (" 1. Maxwell Garnett model: J.C.M. Garnett, Colours in metal glasses and in metallic films, Phil. Trans. R. Soc. Lond. 203, 385-420 [1904]")
    print (" 2. Bruggeman model: D.A.G. Bruggeman, Berechnung verschiedener physikalischer Konstanten von heterogenen Substanzen, Ann. Phys. (Leipzig) 24, 636-679 [1935]")
    print (" 3. Lichtenecker model: K.Lichtenecker,Phys. Z. 27 (1926) 115.")
    print (" 4. Looyenga model: H. Looyenga, Dielectric constants of heterogeneous mixtures, Physica 31, 401-406 [1965]")
    print (" 5. Monecke mode: J. Monecke, Bergman spectral representation of a simple expression for the dielectric response of a symmetric two-component composite, J. Phys.: Cond. Mat. 6, 907-912 [1994]")
    print (" 6. Hollow sphere equivalent mode: C.F. Bohren and D.R. Huffman, Absorption and Scattering of Light by Small Particles, Wiley, New York, p. 149 [1983]")
    print ("")
    print ("")    
    print (" Several methods to predict the mixing material's effective permittivity.")
    print (" If you want to use the method by Maxwell Garnett model -> please press <1>")
    print (" If you want to use the method by Bruggeman model -> please press <2>")
    print (" If you want to use the method by Lichtenecker model -> please press <3>")
    print (" If you want to use the method by Looyengat model -> please press <4>")
    print (" If you wanto use the method by Monecke model -> please press <5>")
    print (" If you wanto use the method by Hollow sphere equivalent model-> please press <6>")

    
    answer = int(input ("\n"))
    break

while (answer != 1) and (answer != 2) and (answer != 3) and (answer != 4) and (answer != 5) and (answer != 6):

    print (" Please type appropriate number.")
    print ("")
    print (" Several methods to predict the mixing material's effective permittivity.")
    print (" If you want to use the method by Maxwell Garnett model -> please press <1>")
    print (" If you want to use the method by Bruggeman model -> please press <2>")
    print (" If you want to use the method by Lichtenecker model -> please press <3>")
    print (" If you want to use the method by Looyengat model -> please press <4>")
    print (" If you wanto use the method by Monecke model -> please press <5>")
    print (" If you wanto use the method by Hollow sphere equivalent model-> please press <6>")
    print (" If you want to see the references, -> please press <0>")
    answer = input ("\n")
    answer= int(answer)
    
    while answer == 0:
        print (" References:")
        print (" 1. Maxwell Garnett model: J.C.M. Garnett, Colours in metal glasses and in metallic films, Phil. Trans. R. Soc. Lond. 203, 385-420 [1904]")
        print (" 2. Bruggeman model: D.A.G. Bruggeman, Berechnung verschiedener physikalischer Konstanten von heterogenen Substanzen, Ann. Phys. (Leipzig) 24, 636-679 [1935]")
        print (" 3. Lichtenecker model: K.Lichtenecker,Phys. Z. 27 (1926) 115.")
        print (" 4. Looyenga model: H. Looyenga, Dielectric constants of heterogeneous mixtures, Physica 31, 401-406 [1965]")
        print (" 5. Monecke model: J. Monecke, Bergman spectral representation of a simple expression for the dielectric response of a symmetric two-component composite, J. Phys.: Cond. Mat. 6, 907-912 [1994]")
        print (" 6. Hollow sphere equivalent model: C.F. Bohren and D.R. Huffman, Absorption and Scattering of Light by Small Particles, Wiley, New York, p. 149 [1983]")
        print ("")
        print ("")    
        print (" Several methods to predict the mixing material's effective permittivity.")
        print (" If you want to use the method by Maxwell Garnett model -> please press <1>")
        print (" If you want to use the method by Bruggeman model -> please press <2>")
        print (" If you want to use the method by Lichtenecker model -> please press <3>")
        print (" If you want to use the method by Looyengat model -> please press <4>")
        print (" If you wanto use the method by Monecke model -> please press <5>")
        print (" If you wanto use the method by Hollow sphere equivalent model-> please press <6>")
        
        answer = int(input ("\n"))
        break
    break

if (answer == 1 or 2 or 3 or 4 or 5 or 6):
    print (" I suggest utilizing two sets of optical data with the same wavelength (or frequency).")
    print (" Alternatively, this code can facilitate the initial calculation/interpolation of your ")
    print (" optical data based on the chosen wavelength (or frequency) from the matrix data or inclusion data")
    print (" In this case, you will need to specify the desired wavelength for utilization.")
    print (" Note that interpolation is proceeding by the algorithm of Python InterpolatedUnivariateSpline.")
    print (" If you are using already two sets of optical data with the same wavelength (or frequency) -> please press <0> ")
    print (" If you want to use the wavelength (or frequency) of matrix -> please press <1> ")
    print (" If you want to use the wavelength (or frequency) of inclusion -> please press <2> ")
    answer_interpolation = int(input ("\n"))

while (answer_interpolation != 0) and (answer_interpolation != 1) and (answer_interpolation != 2):

    print (" Please type appropriate number.")
    print ("")
    print (" If you are using already two sets of optical data with the same wavelength (or frequency) -> please press <0> ")
    print (" If you want to use the wavelength (or frequency) of matrix -> please press <1> ")
    print (" If you want to use the wavelength (or frequency) of inclusion -> please press <2> ")
    answer_interpolation = int(input ("\n"))
    if (answer_interpolation == 0 or 1 or 2):
        break


if (answer == 1 or 2 or 3 or 4 or 5 or 6):
    print (" Please provide the name of the file for the refractive index of the matrix.")
    print (" File must contain three columns, e.g., wavelength (or frequency), n, and k")
    file_matrix = input ("")
    print ("")

    print (" Please  provide the names of the files containing the refractive index for the inclusions.")
    print (" File must contain three columns, e.g., wavelength (or frequency), n, and k")
    file_inclusion = input ("")
    print ("")
    
    if (answer != 6):

        print (" Please provide the volume fraction of the inclusions.")
        print (" Caution: 0 < volume fraction < 1") 
        vol_fr = float(input (""))
        print ("")

    if (answer == 6):
        
        print (" Hollow sphere equivalent model needs a specific parameter")
        print (" Thickness: 1 - ((r_in)**3/(r_out)**3), where")
        print (" r_in is inner and r_out is outer radius of model")
        print (" In this model, the volume fraction is determined by the thickness.")
        print (" Please provide inner radius in unit same to outer radius")
        r_in = float(input (""))
        print ("")
        print (" Please provide inner radius in unit same to inner radius")
        r_out = float(input (""))
        vol_fr = 1 - (r_in**3/r_out**3)
        print ("")

    print (" Please  provide the name of the output file.")
    file_output = input ("")
    print ("")

        
# Interpolation first

    for filename in glob.glob(file_matrix):
        print ("|______________________________________________________________________________|")
        print ("- Summary of your inputs")
        print ("")
        print ("1. Name of the file for the refractive index of the matrix: ", file_matrix)
        f1 = open (file_matrix, "r")

        for line in f1:
            line1 = line.split()
            frequency_or_wavelength_m.append(float(line1[0]))
            n_m.append(float(line1[1]))
            k_m.append(float(line1[2]))

    for filename in glob.glob(file_inclusion):
        print ("2. Name of the file for the refractive index of the matrix: ", file_inclusion)
        f2 = open (file_inclusion, "r")

        for line in f2:
            line1 = line.split()
            frequency_or_wavelength_in.append(float(line1[0]))
            n_in.append(float(line1[1]))
            k_in.append(float(line1[2]))

    if (answer_interpolation == 1):

        interpolation_function = InterpolatedUnivariateSpline(frequency_or_wavelength_in, n_in)
        n_in = interpolation_function(frequency_or_wavelength_m)

        interpolation_function = InterpolatedUnivariateSpline(frequency_or_wavelength_in, k_in)
        k_in = interpolation_function(frequency_or_wavelength_m)

    elif (answer_interpolation == 2):

        interpolation_function = InterpolatedUnivariateSpline(frequency_or_wavelength_m, n_m)
        n_m = interpolation_function(frequency_or_wavelength_in)

        interpolation_function = InterpolatedUnivariateSpline(frequency_or_wavelength_m, k_m)
        k_m = interpolation_function(frequency_or_wavelength_in)


#optical constants
    print ("3. Volume fraction of inclusion: ", vol_fr)
    print ("|______________________________________________________________________________|")
    print ("")
    print ("Being calculated...\n")

    for l in range(len(frequency_or_wavelength_m)):
        e_m1 = n_m[l]*n_m[l] - k_m[l]*k_m[l]
        e_m2 = 2*n_m[l]*k_m[l]
        e_m = complex(e_m1, e_m2)

        # print (e_m)
        e_in1 = n_in[l]*n_in[l] - k_in[l]*k_in[l]
        e_in2 = 2*n_in[l]*k_in[l]
        e_in = complex(e_in1, e_in2)

        if answer == 1:
            e_eff = (e_m * (1 + 2*vol_fr*((e_in - e_m)/(e_in + 2*e_m)))/(1 - vol_fr *((e_in - e_m)/(e_in + 2*e_m)))) #Maxwell Garnett

        if answer == 2:
            e_eff = (1/4*((3*vol_fr(e_in - e_m)) + 2*e_m - e_in + sqrt((((1-3*vol_fr)**2)*(e_in**2)) +
                    (2*(2+9*vol_fr - 9*(vol_fr**2))*e_in*e_m) + (((3*vol_fr - 2)**2) *e_m**2)))) #Bruggeman

        if answer == 3:
            e_eff = e_in**vol_fr + e_m**(1-vol_fr) #Lichtenecker

        if answer == 4:
            e_eff = (vol_fr*e_in**1/3+(1-vol_fr)*e_m**1/3)**3 # Looyenga

        if answer == 5:
           e_eff = ((2*(vol_fr*e_i + (1-vol_fr)*e_m)**2 + e_in*e_m))/((1+vol_fr)*e_in + (2-vol_fr)*e_m) #    Monecke model 
        
        if answer == 6:
           e_eff = e_in * ((3-2*vol_fr)*e_m + 2*f*e_in)/(vol_fr*e_m + (3-vol_fr)*e_in) #    Hollow sphere equivalent model

        constant_n.append(sqrt((sqrt(e_eff.real*e_eff.real + e_eff.imag*e_eff.imag) + e_eff.real)/2))
        constant_k.append(sqrt((sqrt(e_eff.real*e_eff.real + e_eff.imag*e_eff.imag) - e_eff.real)/2))

#writing into files k and n values
    f = open(file_output, "w")

    for l in range(len(frequency_or_wavelength_m)):
    
        f.write(str(frequency_or_wavelength_m[l]))
        f.write(" ")
        f.write(str(constant_n[l]))
        f.write(" ")
        f.write(str(constant_k[l]))
        f.write("\n")
    
    f.close()