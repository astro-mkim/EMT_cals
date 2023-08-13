***EMT calculation***

**Effective Medium Theory (EMT)**

The effective-medium approximation is a method used to analyze macroscopically inhomogeneous mediums, where properties such as conductivity (σ), dielectric function (ϵ), or elastic modulus vary within the space. Several models (Maxwell Garnett model, Bruggeman model, Lichtenecker model, Looyengat model, etc.) are available to predict the effective permittivity of mixed materials. For more details, please visit: https://en.wikipedia.org/wiki/Effective_medium_approximations.

This straightforward code calculates refractive indices (n and k values) based on the refractive indices of both the matrix and the inclusion. In most scenarios, I recommend using two sets of optical data with the same wavelength (or frequency). Alternatively, this code can facilitate initial calculations/interpolation of optical data based on the selected wavelength (or frequency) from either the matrix data or the inclusion data. In this case, you will need to specify your desired wavelength for calculation. Note that interpolation is carried out using Python's InterpolatedUnivariateSpline algorithm.

The file must contain three columns, e.g., wavelength (or frequency), n, and k.

**References:**

- Maxwell-Garnett model (J.C.M. Garnett, "Colours in metal glasses and in metallic films", Phil. Trans. R. Soc. Lond. 203, 385-420 [1904])
- Bruggeman model (D.A.G. Bruggeman, "Berechnung verschiedener physikalischer Konstanten von heterogenen Substanzen", Ann. Phys. (Leipzig) 24, 636-679 [1935])
- Looyenga model (H. Looyenga, "Dielectric constants of heterogeneous mixtures", Physica 31, 401-406 [1965])
- Monecke model (J. Monecke, Bergman spectral representation of a simple expression for the dielectric response of a symmetric two-component composite, J. Phys.: Cond. Mat. 6, 907-912 [1994])
- Hollow sphere equivalent model (C.F. Bohren and D.R. Huffman, Absorption and Scattering of Light by Small Particles, Wiley, New York, p. 149 [1983])


**How to use:**

Prerequisites: Python packages numpy, math (cmath), and scipy

In your terminal: python EMT_calculation.py

You will be guided through an interactive process.
