Copyright  2016, The *ALF* Project.
This is the ALF Project Documentation by the ALF contributors. It is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License. You are free to share and benefit from this documentation as long as this license is preserved and proper attribution to the authors is given. For details see the ALF project homepage [alf.physik.uni-wuerzburg.de](alf.physik.uni-wuerzburg.de).

Introduction
============

The auxiliary field quantum Monte Carlo (QMC) approach is the algorithm of choice to simulate a variety of correlated electron systems in the solid state and beyond . The phenomena one can investigate in detail include correlation effects in in the bulk and surfaces of topological insulators, quantum phase transitions between semimetals (Dirac fermions) and insulators, deconfined quantum critical points, topologically ordered phases, heavy fermion systems, nematic and magnetic quantum phase transitions in metals, superconductivity in spin orbit split bands, SU(N) symmetric models, etc. This ever growing list of phenomena is based on recent symmetry related insights that enable one to find sign free formulations of the problem thus allowing solutions in polynomial time . The aim of this project is to introduce a general formulation of the finite temperature auxiliary field method so as to quickly be able to play with different model Hamiltonians at minimal programming cost. The reader is expected to be familiar with the auxiliary field QMC approach. A detailed review containing all the prerequisites for understanding the code can be found in . In this documentation, we will briefly list the most important equations of the auxiliary field QMC and then show in all details how to implement a variety of models, run the code, and produce results for equal time and time displaced correlation functions. The program code is written in Fortran according to the 2003 standard.

The first and most important part is to define a general Hamiltonian that can accommodate a large class of models (see Sec. \[sec:def\]). Our approach is to express the model as a sum of one-body terms, a sum of two-body terms each written as a perfect square of a one one body term, as well as one-body term coupled to an Ising field with dynamics to be specified by the user. The form of the interaction in terms of sums of perfect squares allows us to use generic forms of discrete approximations to the Hubbard-Stratonovich (HS) transformation. Symmetry considerations are imperative to enhance the speed of the code. We thereby include a <span>*color*</span> index reflecting an underlying SU(N) color symmetry as well as a flavor index reflecting the fact that after the HS transformation, the fermionic determinant is block diagonal in this index. To use the code, one will require a minimal understanding of the algorithm. In Section \[sec:def\], we very briefly show how to go through the steps required to formulated the many body imaginary time propagation in terms of a sum over HS and Ising fields of one body imaginary time propagator. The user will have to provide this one body imaginary time propagator for a given configuration of HS and Ising fields.

Section \[sec:imp\] is devoted to the data structures which need to implement the model. This includes an `Operator` type to optimally work with sparse Hermitian matrices, a `Lattice` type to define one and two dimensional Bravais lattices, and two `Observable` types to handle equal time, time displaced and scalar observables.

After a description of the file structure in Sec \[sec:files\], we will give explicit examples on how to use the code for the Hubbard model on square and honeycomb lattices for different choices of the Hubbard Stratonovich transformation (See Secs. \[sec:walk1\], \[sec:walk1.1\] and  \[sec:walk1.2\]) as well as the Hubbard model on a square lattice coupled to a transverse Ising field (see Sec. \[sec:walk2\] ). Our implementation is rather general such that a variety of other models can be simulated. In Sec. \[sec:other\_models\] we provide some information on how to simulate the Kondo lattice as well as the SU(N) symmetric Hubbard-Heisenberg model.

The Monte Carlo run and the data analysis are separate: the QMC run dumps the results of <span>*bins*</span> sequentially into files which are then analyzed by analysis programs. In Sec. \[sec:analysis\], we provide a brief description of the analysis programs for our three observable types. The analysis program allow for omitting a given number of initial bins so as to allow for warmup and to rebin so as to a posteriori take into account long autocorrelation times. Finally, Sec. \[sec:running\] will provide all details required to compile and run the code.

Definition of the model Hamiltonian
===================================

The class of solvable models includes Hamiltonians $\\hat{\\mathcal{H}}$ that have the following general form:
$$\\begin{aligned}
\\hat{\\mathcal{H}}&=&\\hat{\\mathcal{H}}\_{T}+\\hat{\\mathcal{H}}\_{V} +  \\hat{\\mathcal{H}}\_{I} +   \\hat{\\mathcal{H}}\_{0,I}\\;,\\mathrm{where}
\\label{eqn:general\_ham}\\\\
\\hat{\\mathcal{H}}\_{T}
&=&
\\sum\\limits\_{k=1}^{M\_T}
\\sum\\limits\_{\\sigma=1}^{N\_{\\mathrm{col}}}
\\sum\\limits\_{s=1}^{N\_{\\mathrm{fl}}}
\\sum\\limits\_{x,y}^{N\_{\\mathrm{dim}}}
\\hat{c}^{\\dagger}\_{x \\sigma   s}T\_{xy}^{(k s)} \\hat{c}^{\\phantom\\dagger}\_{y \\sigma s}  \\equiv  \\sum\\limits\_{k=1}^{M\_T} \\hat{T}^{(k)}
\\label{eqn:general\_ham\_t}\\\\
\\hat{\\mathcal{H}}\_{V}
&=&
\\sum\\limits\_{k=1}^{M\_V}U\_{k}
\\left\\{
\\sum\\limits\_{\\sigma=1}^{N\_{\\mathrm{col}}}
\\sum\\limits\_{s=1}^{N\_{\\mathrm{fl}}}
\\left\[
\\left(
\\sum\\limits\_{x,y}^{N\_{\\mathrm{dim}}}
\\hat{c}^{\\dagger}\_{x \\sigma s}V\_{xy}^{(k s)}\\hat{c}^{\\phantom\\dagger}\_{y \\sigma s}
\\right)
+\\alpha\_{k s} 
\\right\]
\\right\\}^{2}  \\equiv   
\\sum\\limits\_{k=1}^{M\_V}U\_{k}   \\left(\\hat{V}^{(k)} \\right)^2
\\label{eqn:general\_ham\_v}\\\\
\\hat{\\mathcal{H}}\_{I}  & = &
\\sum\\limits\_{k=1}^{M\_I} \\hat{Z}\_{k} 
\\left(
\\sum\\limits\_{\\sigma=1}^{N\_{\\mathrm{col}}}
\\sum\\limits\_{s=1}^{N\_{\\mathrm{fl}}}
\\sum\\limits\_{x,y}^{N\_{\\mathrm{dim}}}
\\hat{c}^{\\dagger}\_{x \\sigma s} I\_{xy}^{(k s)}\\hat{c}^{\\phantom\\dagger}\_{y \\sigma s}
\\right) \\equiv \\sum\\limits\_{k=1}^{M\_I} \\hat{Z}\_{k}    \\hat{I}^{(k)} 
\\;.\\label{eqn:general\_ham\_i}\\end{aligned}$$
 The indices have the following meaning:

-   The number of fermion *flavors* is set by *N*<sub>*f**l*</sub>. After the Hubbard-Stratonovich transformation, the action will be block diagonal in the flavor index.

-   The number of fermion *colors* is set by *N*<sub>*c**o**l*</sub>. The Hamiltonian is invariant under SU(*N*<sub>*c**o**l*</sub>) rotations. Note that in the code *N*<sub>*c**o**l*</sub> ≡ `N_``SUN`.

-   The indices *x*, *y* label lattice sites where *x*, *y* = 1, ⋯, *N*<sub>*d**i**m*</sub>.

    *N*<sub>*d**i**m*</sub> is the total number of spacial vertices: *N*<sub>*d**i**m*</sub> = *N*<sub>unit cell</sub>*N*<sub>*o**r**b**i**t**a**l*</sub>, where *N*<sub>unit cell</sub> is the number of unit cells of the underlying Bravais lattice and *N*<sub>*o**r**b**i**t**a**l*</sub> is the number of (spacial) orbitals per unit cell.

-   Therefore, the matrices **T**<sup>(*k**s*)</sup>, **V**<sup>(*k**s*)</sup> and **I**<sup>(*k**s*)</sup> are of dimension *N*<sub>*d**i**m*</sub> × *N*<sub>*d**i**m*</sub>

-   The number of interaction terms is labelled by *M*<sub>*V*</sub> and *M*<sub>*I*</sub>. *M*<sub>*T*</sub> &gt; 1 would allow for a checkerboard decomposition.

The Ising part of the general Hamiltonian (\[eqn:general\_ham\]) is $\\hat{\\mathcal{H}}\_{0,I}+ \\hat{\\mathcal{H}}\_{I}$ and has the following properties:

-   $\\hat{Z}\_k$ is an Ising spin operator which corresponds to the Pauli matrix $\\hat{\\sigma}\_{z}$. It couples to a general one-body term.

-   The dynamics of the Ising spins is given by $\\hat{\\mathcal{H}}\_{0,I}$. This term is not specified here; it has to be specified by the user and is only relevant when the Monte Carlo update probability is computed in the code (see Sec. \[sec:walk2\] for an example).

Note that the matrices **T**<sup>(*k**s*)</sup>, **V**<sup>(*k**s*)</sup> and **I**<sup>(*k**s*)</sup> explicitly depend on the flavor index *s* but not on the color index *σ*. The color index *σ* only appears in the second quantized operators such that the Hamiltonian is manifestly SU(*N*<sub>*c**o**l*</sub>) symmetric. We also require the matrices **T**<sup>(*k**s*)</sup>, **V**<sup>(*k**s*)</sup> and **I**<sup>(*k**s*)</sup> to be Hermitian.

Formulation of the QMC
----------------------

The formulation of the Monte Carlo simulation is based on the following.

-   We will discretize the imaginary time propagation: *β* = *Δ**τ**L*<sub>Trotter</sub>

-   We will use the discrete Hubbard-Stratonovich transformation:
    $$\\label{HS\_squares}
            e^{\\Delta \\tau  \\lambda  \\hat{A}^2 } =
            \\sum\_{ l = \\pm 1, \\pm 2}  \\gamma(l)
    e^{ \\sqrt{\\Delta \\tau \\lambda }
           \\eta(l)  \\hat{A} }
                    + {\\cal O} (\\Delta \\tau ^4)\\;,$$
     where the fields *η* and *γ* take the values:
    $$\\begin{aligned}
     \\gamma(\\pm 1)  = 1 + \\sqrt{6}/3, \\quad  \\eta(\\pm 1 ) = \\pm \\sqrt{2 \\left(3 - \\sqrt{6} \\right)}\\;,\\\\\\nonumber
      \\gamma(\\pm 2) = 1 - \\sqrt{6}/3, \\quad  \\eta(\\pm 2 ) = \\pm \\sqrt{2 \\left(3 + \\sqrt{6} \\right)}\\;.
    \\nonumber\\end{aligned}$$

-   We will work in a basis for the Ising spins where $\\hat{Z}\_k$ is diagonal: $\\hat{Z}\_{k}|s\_{k}\\rangle = s\_{k}|s\_{k}\\rangle$, where *s*<sub>*k*</sub> = ±1.

-   From the above it follows that the Monte Carlo configuration space *C* is given by the combined spaces of Ising spin configurations and of Hubbard-Stratonovich discrete field configurations:
    *C* = {*s*<sub>*i*, *τ*</sub>,*l*<sub>*j*, *τ*</sub> with *i*=1⋯*M*<sub>*I*</sub>, *j*=1⋯*M*<sub>*V*</sub>, *τ*=1⋯*L*<sub>*T**r**o**t**t**e**r*</sub>}
     Here, the Ising spins take the values *s*<sub>*i*, *τ*</sub> = ±1 and the Hubbard-Stratonovich fields take the values *l*<sub>*j*, *τ*</sub> = ±2, ±1.

### The partition function

With the above, the partition function of the model (\[eqn:general\_ham\]) can be written as follows.
$$\\begin{aligned}
\\label{eqn:partition\_1}
Z &=& \\Tr{\\left(e^{-\\beta \\hat{\\mathcal{H}} }\\right)}\\nonumber\\\\
  &=&   \\Tr{  \\left\[ e^{-\\Delta \\tau \\hat{\\mathcal{H}}\_{0,I}}   \\prod\_{k=1}^{M\_T}   e^{-\\Delta \\tau \\hat{T}^{(k)}}  
    \\prod\_{k=1}^{M\_V}   e^{ - \\Delta \\tau  U\_k \\left(  \\hat{V}^{(k)} \\right)^2}   \\prod\_{k=1}^{M\_I}   e^{  -\\Delta \\tau  \\hat{\\sigma}\_{k}  \\hat{I}^{(k)}} 
   \\right\]^{L\_{\\text{Trotter}}}}  + \\mathcal{O}(\\Delta\\tau^{2})\\nonumber \\\\
   &=&
   \\sum\_{C} \\left( \\prod\_{k=1}^{M\_V} \\prod\_{\\tau=1}^{L\_{\\mathrm{Trotter}}} \\gamma\_{k,\\tau} \\right) e^{-S\_{0,I} \\left( \\left\\{ s\_{i,\\tau} \\right\\}  \\right) }\\times \\nonumber\\\\
   &\\quad&
    \\Trf{ \\left\\{  \\prod\_{\\tau=1}^{L\_{\\mathrm{Trotter}}} \\left\[   \\prod\_{k=1}^{M\_T}   e^{-\\Delta \\tau \\hat{T}^{(k)}}  
    \\prod\_{k=1}^{M\_V}   e^{  \\sqrt{ -\\Delta \\tau  U\_k} \\eta\_{k,\\tau} \\hat{V}^{(k)} }   \\prod\_{k=1}^{M\_I}   e^{  -\\Delta \\tau s\_{k,\\tau}  \\hat{I}^{(k)}}  \\right\]\\right\\} }+ \\mathcal{O}(\\Delta\\tau^{2})\\;.\\end{aligned}$$
 In the above, the trace *T**r* runs over the Ising spins as well as over the fermionic degrees of freedom, and *T**r*<sub>*F*</sub> only over the fermionic Fock space. *S*<sub>0, *I*</sub>({*s*<sub>*i*, *τ*</sub>}) is the action corresponding to the Ising Hamiltonian, and is only dependent on the Ising spins so that it can be pulled out of the fermionic trace. At this point, and since for a given configuration *C* we are dealing with a free propagation, we can integrate out the fermions to obtain a determinant:
$$\\begin{aligned}
 &\\quad&\\Trf{ \\left\\{  \\prod\_{\\tau=1}^{L\_{\\mathrm{Trotter}}} \\left\[   \\prod\_{k=1}^{M\_T}   e^{-\\Delta \\tau \\hat{T}^{(k)}}  
    \\prod\_{k=1}^{M\_V}   e^{  \\sqrt{ - \\Delta \\tau  U\_k} \\eta\_{k,\\tau} \\hat{V}^{(k)} }   \\prod\_{k=1}^{M\_I}   e^{  -\\Delta \\tau s\_{k,\\tau}  \\hat{I}^{(k)}}  \\right\] \\right\\}} = \\nonumber\\\\
&\\quad& \\quad\\prod\\limits\_{s=1}^{N\_{\\mathrm{fl}}} \\left\[  e^{\\sum\\limits\_{k=1}^{M\_V} \\sum\\limits\_{\\tau = 1}^{L\_{\\mathrm{Trotter}}}\\sqrt{-\\Delta \\tau U\_k}  \\alpha\_{k,s} \\eta\_{k,\\tau} }
   \\right\]^{N\_{\\mathrm{col}}}\\times
\\nonumber\\\\
&\\quad&\\quad   \\prod\\limits\_{s=1}^{N\_{\\mathrm{fl}}} 
   \\left\[
    \\det\\left(  1 + 
     \\prod\_{\\tau=1}^{L\_{\\mathrm{Trotter}}}   \\prod\_{k=1}^{M\_T}   e^{-\\Delta \\tau {\\bf T}^{(ks)}}  
    \\prod\_{k=1}^{M\_V}   e^{  \\sqrt{ -\\Delta \\tau  U\_k} \\eta\_{k,\\tau} {\\bm V}^{(ks)} }   \\prod\_{k=1}^{M\_I}   e^{  -\\Delta \\tau s\_{k,\\tau}  {\\bm I}^{(ks)}}  
     \\right) \\right\]^{N\_{\\mathrm{col}}}\\;,\\end{aligned}$$
 where the matrices ${\\bf T}^{(ks)}$, ${\\bf V}^{(ks)}$, and ${\\bf I}^{(ks)}$ define the Hamiltonian \[Eq. (\[eqn:general\_ham\]) - (\[eqn:general\_ham\_i\])\]. All in all, the partition function is given by:
$$\\begin{aligned}
\\label{eqn:partition\_2}
    Z  &=&   \\sum\_{C}   e^{-S\_{0,I} \\left( \\left\\{ s\_{i,\\tau} \\right\\}  \\right) }     \\left( \\prod\_{k=1}^{M\_V} \\prod\_{\\tau=1}^{L\_{\\mathrm{Trotter}}} \\gamma\_{k,\\tau} \\right)
    e^{ N\_{\\mathrm{col}}\\sum\\limits\_{s=1}^{N\_{\\mathrm{fl}}} \\sum\\limits\_{k=1}^{M\_V} \\sum\\limits\_{\\tau = 1}^{L\_{\\mathrm{Trotter}}}\\sqrt{-\\Delta \\tau U\_k}  \\alpha\_{k,s} \\eta\_{k,\\tau} } 
  \\times   \\nonumber \\\\
  &\\quad&
      \\prod\_{s=1}^{N\_{\\mathrm{fl}}}\\left\[\\det\\left(  1 + 
     \\prod\_{\\tau=1}^{L\_{\\mathrm{Trotter}}}   \\prod\_{k=1}^{M\_T}   e^{-\\Delta \\tau {\\bm T}^{(ks)}}  
    \\prod\_{k=1}^{M\_V}   e^{  \\sqrt{ -\\Delta \\tau  U\_k} \\eta\_{k,\\tau} {\\bm V}^{(ks)} }   \\prod\_{k=1}^{M\_I}   e^{  -\\Delta \\tau s\_{k,\\tau}  {\\bm I}^{(ks)}}  
     \\right) \\right\]^{N\_{\\mathrm{col}}}  \\nonumber \\\\ 
     & \\equiv&  \\sum\_{C} e^{-S(C) }\\;.\\end{aligned}$$
 In the above, one notices that the weight factorizes in the flavor index. The color index raises the determinant to the power *N*<sub>*c**o**l*</sub>. This corresponds to an explicit *S**U*(*N*<sub>*c**o**l*</sub>) symmetry for each configuration. This symmetry is manifest in the fact that the single particle Green functions are color independent, again for each given configuration *C*.

### Observables

In the auxiliary field QMC approach, the single particle Green function plays a crucial role. It determines the Monte Carlo dynamics and is used to compute observables:
$$\\label{eqn:obs}
\\langle \\hat{O}  \\rangle  = \\frac{ \\text{Tr}   \\left\[ e^{- \\beta \\hat{H}}  \\hat{O}   \\right\] }{ \\text{Tr}   \\left\[ e^{- \\beta \\hat{H}}  \\right\] } =   \\sum\_{C}   P(C) 
   \\langle \\langle \\hat{O}  \\rangle \\rangle\_{(C)} , \\text{   with   } 
  P(C)   = \\frac{ e^{-S(C)}}{\\sum\_C e^{-S(C)}}\\;,$$
 and $\\langle \\langle \\hat{O}  \\rangle \\rangle\_{(C)} $ denotes the observed value of $\\hat{O}$ for a given configuration *C*. For a given configuration *C* one can use Wicks theorem to compute *O*(*C*) from the knowledge of the single particle Green function:
$$G( x,\\sigma,s, \\tau |    x',\\sigma',s', \\tau')   =       \\langle \\langle {\\cal T} \\hat{c}^{\\phantom\\dagger}\_{x \\sigma s} (\\tau)  \\hat{c}^{\\dagger}\_{x' \\sigma' s'} (\\tau') \\rangle \\rangle\_{C}$$
 where $ {\\cal T} $ corresponds to the imaginary time ordering operator. The corresponding equal time quantity reads,
$$G( x,\\sigma,s, \\tau |    x',\\sigma',s', \\tau)   =       \\langle \\langle {\\cal T} \\hat{c}^{\\phantom\\dagger}\_{x \\sigma s} (\\tau)  \\hat{c}^{\\dagger}\_{x' \\sigma' s'} (\\tau) \\rangle \\rangle\_{C}$$
 Since for a given HS field translation invariance in imaginary time is broken, the Green function has an explicit *τ* and *τ*′ dependence. On the other hand it is diagonal in the flavor index, and independent on the color index. The later reflects the explicit SU(N) color symmetry present at the level of individual HS configurations.

To compute equal time as well as time-displaced observables, one can make use of Wicks theorem. A convenient formulation of this theorem for QMC simulations reads:
$$\\begin{aligned}
& & \\langle \\langle 	{\\cal T}   c^{\\dagger}\_{\\underline x\_{1}}(\\tau\_{1}) c^{\\phantom\\dagger}\_{{\\underline x}'\_{1}}(\\tau'\_{1})  
\\cdots c^{\\dagger}\_{\\underline x\_{n}}(\\tau\_{n}) c^{\\phantom\\dagger}\_{{\\underline x}'\_{n}}(\\tau'\_{n}) 
\\rangle \\rangle\_{C} =
\\nonumber \\\\ & & \\det  
\\begin{bmatrix}
   \\langle \\langle   {\\cal T}   c^{\\dagger}\_{\\underline x\_{1}}(\\tau\_{1}) c^{\\phantom\\dagger}\_{{\\underline x}'\_{1}}(\\tau'\_{1})  \\rangle \\rangle\_{C} & 
    \\langle \\langle  {\\cal T}   c^{\\dagger}\_{\\underline x\_{1}}(\\tau\_{1}) c^{\\phantom\\dagger}\_{{\\underline x}'\_{2}}(\\tau'\_{2})  \\rangle \\rangle\_{C}  & \\dots   &   
    \\langle \\langle   {\\cal T}   c^{\\dagger}\_{\\underline x\_{1}}(\\tau\_{1}) c^{\\phantom\\dagger}\_{{\\underline x}'\_{n}}(\\tau'\_{n})  \\rangle \\rangle\_{C}  \\\\
    \\langle \\langle   {\\cal T}   c^{\\dagger}\_{\\underline x\_{2}}(\\tau\_{2}) c^{\\phantom\\dagger}\_{{\\underline x}'\_{1}}(\\tau'\_{1})  \\rangle \\rangle\_{C}  & 
      \\langle \\langle   {\\cal T}   c^{\\dagger}\_{\\underline x\_{2}}(\\tau\_{2}) c^{\\phantom\\dagger}\_{{\\underline x}'\_{2}}(\\tau'\_{2})  \\rangle \\rangle\_{C}  & \\dots  &
       \\langle \\langle   {\\cal T}   c^{\\dagger}\_{\\underline x\_{2}}(\\tau\_{2}) c^{\\phantom\\dagger}\_{{\\underline x}'\_{n}}(\\tau'\_{n})  \\rangle \\rangle\_{C}   \\\\
    \\vdots & \\vdots &  \\ddots & \\vdots \\\\
    \\langle \\langle   {\\cal T}   c^{\\dagger}\_{\\underline x\_{n}}(\\tau\_{n}) c^{\\phantom\\dagger}\_{{\\underline x}'\_{1}}(\\tau'\_{1})  \\rangle \\rangle\_{C}   & 
     \\langle \\langle   {\\cal T}   c^{\\dagger}\_{\\underline x\_{n}}(\\tau\_{n}) c^{\\phantom\\dagger}\_{{\\underline x}'\_{2}}(\\tau'\_{2})  \\rangle \\rangle\_{C}   & \\dots  & 
     \\langle \\langle   {\\cal T}   c^{\\dagger}\_{\\underline x\_{n}}(\\tau\_{n}) c^{\\phantom\\dagger}\_{{\\underline x}'\_{n}}(\\tau'\_{n})  \\rangle \\rangle\_{C}
 \\end{bmatrix}\\end{aligned}$$
 In the subroutines and of the module (see Sec. \[sec:obs\]) the user is provided with the equal time and time displaced correlation function. Using the above formulation of Wicks theorem, arbitrary correlation functions can be computed. We note however, that the program is limited to the calculation of observables that contain only two different imaginary times.

### Reweighting and the sign problem

In general, the action *S*(*C*) will be complex, thereby inhibiting a direct Monte Carlo sampling of *P*(*C*). This leads to the infamous sign problem. When the average sign is not too small, we can nevertheless compute observables within a reweighting scheme. Here we adopt the following scheme. First note that the partition function is real such that:
$$Z =   \\sum\_{C}  e^{-S(C)}    =  \\sum\_{C}  \\overline{e^{-S(C)}} = \\sum\_{C}  \\Re \\left\[e^{-S(C)} \\right\].$$
 Thereby[1] and with the definition
$$\\label{Sign.eq}
	 \\text{ sign }(C)   =  \\frac{   \\Re \\left\[e^{-S(C)} \\right\]  } {\\left| \\Re \\left\[e^{-S(C)} \\right\]  \\right|  }\\;,$$
 the computation of the observable \[Eq. (\[eqn:obs\])\] is re-expressed as follows:
$$\\begin{aligned}
\\label{eqn:obs\_rw}
\\langle \\hat{O}  \\rangle  &=&  \\frac{\\sum\_{C}  e^{-S(C)} \\langle \\langle \\hat{O}  \\rangle \\rangle\_{(C)} }{\\sum\_{C}  e^{-S(C)}}       \\nonumber \\\\ 
                          &=&  \\frac{\\sum\_{C}   \\Re \\left\[e^{-S(C)} \\right\]    \\frac{e^{-S(C)}} {\\Re \\left\[e^{-S(C)} \\right\]}  \\langle \\langle \\hat{O}  \\rangle \\rangle\_{(C)} }{\\sum\_{C}   \\Re \\left\[e^{-S(C)} \\right\]}    \\nonumber \\\\ 
          &=&
   \\frac{
     \\left\\{
      \\sum\_{C}  \\left| \\Re \\left\[e^{-S(C)} \\right\]  \\right|   \\text{ sign }(C)   \\frac{e^{-S(C)}} {\\Re \\left\[e^{-S(C)} \\right\]}  \\langle \\langle \\hat{O}  \\rangle \\rangle\_{(C)}  \\right\\}/
            \\sum\_{C}  \\left| \\Re \\left\[ e^{-S(C)} \\right\] \\right|  
          }  
          { 
          \\left\\{ \\sum\_{C}  \\left|  \\Re \\left\[ e^{-S(C)} \\right\]   \\right|   \\text{ sign }(C) \\right\\}/
            \\sum\_{C}   \\left| \\Re \\left\[ e^{-S(C)} \\right\] \\right|  
          } \\nonumber\\\\
          &=&
  	 \\frac{  \\left\\langle  \\text{ sign }   \\frac{e^{-S}} {\\Re \\left\[e^{-S} \\right\]}  \\langle \\langle \\hat{O}  \\rangle \\rangle  \\right\\rangle\_{\\overline{P}} } { \\langle \\text{sign}   \\rangle\_{\\overline{P}}}  \\;.      \\end{aligned}$$
 The average sign is
$$\\label{eqn:sign\_rw}
	 \\langle \\text{sign} \\rangle\_{\\overline{P}} =    \\frac { \\sum\_{C}  \\left|  \\Re \\left\[ e^{-S(C)} \\right\]   \\right|   \\text{ sign }(C) }  {  \\sum\_{C}   \\left| \\Re \\left\[ e^{-S(C)} \\right\] \\right|  } \\;,$$
 and we have $\\langle \\text{sign} \\rangle\_{\\overline{P}} \\in \\mathbb{R}$ per definition. According to Eq. (\[eqn:obs\_rw\]) and Eq. (\[eqn:sign\_rw\]), the Monte Carlo simulation samples the probability distribution
$$\\overline{P}(C) = \\frac{ \\left|  \\Re \\left\[ e^{-S(C)} \\right\] \\right| }{\\sum\_C \\left|  \\Re \\left\[ e^{-S(C)} \\right\]  \\right| }\\;.$$

Implementation of the model
===========================

In general, the module defines the model Hamiltonian, the lattice under consideration and the desired observables (Table \[table:hamiltonian\]). We have collected a number of example Hamiltonians, lattices and observables in the file . They are described in the Sec. \[sec:walk1\] - \[sec:walk2\]. To implement a user-defined model, only the module has to be set up. Accordingly, this documentation focusses almost entirely on this module and the subprograms it includes. The remaining parts of the code may be treated as as a black box.

To specify the Hamiltonian, one needs an and type as well as a type for the observables. These three data structures will be described in the following.

| Subprogram | Description                                                                                                             | Section                      |
|:-----------|:------------------------------------------------------------------------------------------------------------------------|:-----------------------------|
|            | Reads in model and lattice parameters from the file `parameters`.                                                       |                              |
|            | And it sets the Hamiltonian by calling `Ham_latt`, `Ham_hop`, and `Ham_V`.                                              |                              |
|            | Sets the hopping term $\\hat{\\mathcal{H}}\_{T}$ by calling `Op_make` and `Op_set`.                                     | \[sec:op\], \[sec:specific\] |
|            | Sets the interaction terms $\\hat{\\mathcal{H}}\_{V}$ and $\\hat{\\mathcal{H}}\_{I}$ by calling `Op_make` and `Op_set`. | \[sec:op\], \[sec:specific\] |
|            | Sets the lattice by calling `Make_Lattice`.                                                                             | \[sec:latt\]                 |
|            | A function which returns an update ratio for the Ising term $\\hat{\\mathcal{H}}\_{I,0}$.                               | \[sec:s0\]                   |
|            | Asigns memory storage to the observables                                                                                |                              |
|            | Computes the scalar observables and equal-time correlation functions.                                                   | \[sec:obs\]                  |
|            | Computes time-displaced correlation functions.                                                                          | \[sec:obs\]                  |
| `Init_obs` | Initializes the observables to zero.                                                                                    |                              |
| `Pr_obs`   | Writes the observables to the disk by calling `Print_bin`.                                                              |                              |

The `Operator` type
-------------------

The fundamental data structure in the code is the derived data type . This type is used to define the Hamiltonian (\[eqn:general\_ham\]). In general, the matrices ${\\bf T}^{(ks)}$, ${\\bf V}^{(ks)}$ and ${\\bf I}^{(ks)}$ are sparse Hermitian matrices. Consider the matrix ${\\bm X}$ of dimension *N*<sub>*d**i**m*</sub> × *N*<sub>*d**i**m*</sub>, as an representative of each of the above three matrices. Let us denote with {*z*<sub>1</sub>,⋯,*z*<sub>*N*</sub>} a subset of *N* indices, for which
$$X\_{x,y}  =
\\left\\{\\begin{matrix}  X\_{x,y}  &  \\text{ if }   x,  y  \\in \\left\\{ z\_1, \\cdots z\_N \\right\\}\\\\ 
                                  0         &  \\text{ otherwise } 
      \\end{matrix}\\right.$$
 We define the *N* × *N*<sub>*d**i**m*</sub> matrices **P** as
*P*<sub>*i*, *x*</sub> = *δ*<sub>*z*<sub>*i*</sub>, *x*</sub> ,
 where *i* ∈ \[1, ⋯, *N*\] and *x* ∈ \[1, ⋯, *N*<sub>*d**i**m*</sub>\]. The matrix **P** picks out the non-vanishing entries of **X**, which are contained in the rank-*N* matrix **O**. Thereby:
**X** = **P**<sup>*T*</sup>**O****P** ,
 such that:
$$X\_{x,y} = \\sum\\limits\_{i,j}^{N}  P\_{i,x}  O\_{i,j} P\_{j,y}=\\sum\\limits\_{i,j}^{N} \\delta\_{z\_{i},x}  O\_{ij} \\delta\_{z\_{j},y} \\;.$$
 Since the **P** matrices have only one non-vanishing entry per column, they can be stored as a vector $\\vec{P}$:
*P*<sub>*i*</sub> = *z*<sub>*i*</sub>.
 There are many useful identities which emerge from this structure. For example:
$$e^{\\bm{X}} =  e^{\\bm{P}^{T} \\bm{O} \\bm{P}}   = \\sum\_{n=0}^{\\infty}  \\frac{\\left( \\bm{P}^{T} \\bm{O} \\bm{P} \\right)^n}{n!} =  \\bm{P}^{T} e^{ \\bm{O} } \\bm{P}$$
 since
**P****P**<sup>*T*</sup> = 1<sub>*N* × *N*</sub>.

In the code, we define a structure called to capture the above. This type bundles several components that are needed to define and use an operator matrix in the program.

### Specification of the model

| Variable          | Type    | Description                                                   |
|:------------------|:--------|:--------------------------------------------------------------|
|                   | Integer | effective dimension *N*                                       |
|                   | Complex | matrix **O** of dimension *N* × *N*                           |
|                   | Integer | projection matrix **P** encoded as a vector of dimension *N*. |
|                   | Complex | coupling strength *g*                                         |
|                   | Complex | constant *α*                                                  |
|                   | Integer | parameter to set the type of HS transformation                |
|                   |         | (1 = Ising, 2 = Discrete HS, for perfect square)              |
| `Op_X%U`          | Complex | matrix containing the eigenvectors of **O**                   |
| `Op_X%E`          | Real    | eigenvalues of **O**                                          |
| `Op_X%N_non_zero` | Integer | number of non-vanishing eigenvalues of **O**                  |

In this section we show how to specify the Hamiltonian (\[eqn:general\_ham\]) in the code. More precisely, we specify the Hamiltonian by setting the matrices $ e^{-\\Delta \\tau {\\bm T}^{(ks)}}$, $e^{  \\sqrt{ \\Delta \\tau  U\_k} \\eta\_{k\\tau} {\\bm V}^{(ks)} }$, and $e^{  -\\Delta \\tau s\_{k\\tau}  {\\bm I}^{(ks)}}$ that appear in the partition function (\[eqn:partition\_2\]). To do so, we consider the general expression
*e*<sup>*g* *ϕ*<sub>*k**τ*</sub>(`type`) (**X**+*α*)</sup> ,
 and store the following quantities in a variable of type (see Table \[table:operator\]): the matrix **X** \[see Eq. (\[eqn:xeqpdop\])\], the constants *g* and *α*, and, optionally, also the type of the fields *ϕ*<sub>*k**τ*</sub>(type). Either the fields stem from the Ising term (: *ϕ*<sub>*k**τ*</sub> = *s*<sub>*k**τ*</sub>), or they result from the discrete Hubbard-Stratonovich transformation (`type=2`: *ϕ*<sub>*k**τ*</sub> = *η*<sub>*k**τ*</sub>)). In general, we need several arrays of variables of type . Since the implementation exploits the *S**U*(*N*<sub>*c**o**l*</sub>) invariance of the Hamiltonian, we have dropped the color index *σ* in the following.

-   For the hopping Hamiltonian (\[eqn:general\_ham\_t\]), we have to set the exponentiated hopping matrices $ e^{-\\Delta \\tau {\\bm T}^{(ks)}}$:

    In this case **X**<sup>(*k**s*)</sup> = **T**<sup>(*k**s*)</sup>. Precisely, a single variable `Op_T` describes the operator matrix
    $$\\left( \\sum\_{x,y}^{N\_{\\mathrm{dim}}} \\hat{c}^{\\dagger}\_x T\_{xy}^{(ks)} \\hat{c}^{\\phantom{\\dagger}}\_{y}  \\right)  \\;,$$
     where *k* = \[1, *M*<sub>*T*</sub>\] and *s* = \[1, *N*<sub>*f**l*</sub>\]. To make contact with the general expression (\[eqn:exponent\_mat\]) we set *g* = −*Δ**τ*, *α* = 0. In case of the hopping matrix, the type variable `Op_T%type` is neglected by the code. All in all, the corresponding array of structure variables is `Op_T(M_T,N_{fl})`.

-   For the interaction Hamiltonian (\[eqn:general\_ham\_v\]), which is of perfect-square type, we have to set the exponentiated matrices $e^{  \\sqrt{ -  \\Delta \\tau  U\_k} \\eta\_{k\\tau} {\\bm V}^{(ks)} }$:

    In this case, ${\\mathbf X}  = \\mathbf{V}^{(ks)}$. A single variable `Op_V` describes the operator matrix:
    $$\\left\[ \\left( \\sum\_{x,y}^{N\_{\\mathrm{dim}}} \\hat{c}^{\\dagger}\_x V\_{x,y}^{(ks)} \\hat{c}^{\\phantom{\\dagger}}\_{y}  \\right)  + \\alpha\_{ks} \\right\]  \\;,$$
     where *k* = \[1, *M*<sub>*V*</sub>\] and *s* = \[1, *N*<sub>*f**l*</sub>\]. To make contact with the general expression (\[eqn:exponent\_mat\]), we set *α* = *α*<sub>*k**s*</sub> and $g = \\sqrt{-\\Delta \\tau  U\_k}$. The discrete Hubbard-Stratonovich decomposition which is used for the perfect-square interaction, is selected by setting the type variable to `Op_V%type` = 2. All in all, the required structure variables `Op_V` are defined using the array `Op_V(M_V,N_{fl})`.

-   For the Ising interaction Hamiltonian (\[eqn:general\_ham\_i\]), we have to set the exponentiated matrices $e^{  -\\Delta \\tau s\_{k\\tau}  {\\bm I}^{(ks)}}$:

    In this case, **X** = **I**<sup>(*k*, *s*)</sup>. A single variable `Op_V` then describes the operator matrix:
    $$\\left( \\sum\_{x,y}^{N\_{\\mathrm{dim}}} \\hat{c}^{\\dagger}\_x I\_{xy}^{(ks)} \\hat{c}^{\\phantom{\\dagger}}\_{y}  \\right)  \\;,$$
     where *k* = \[1, *M*<sub>*I*</sub>\] and *s* = \[1, *N*<sub>*f**l*</sub>\]. To make contact with the general expression (\[eqn:exponent\_mat\]), we set *α* = 0 and *g* = −*Δ**τ*. The Ising interaction is specified by setting the type variable `Op_V%type=1`. All in all, the required structure variables are contained in the array `Op_V(M_{I},N_{fl})`.

-   In case of a full interaction \[perfect-square term (\[eqn:general\_ham\_v\]) and Ising term (\[eqn:general\_ham\_i\])\], we define the corresponding doubled array `Op_V(M_V+M_I,N_{fl}) ` and set the variables separately for both ranges of the array according to the above.

The `Lattice` type
------------------

We have a lattice module which can generate one and two dimensional Bravais lattices. Note that the orbital structure of each unit cell, has to be specified by the user in the Hamiltonian module. The user has to specify unit vectors $\\vec{a}\_1$ and $\\vec{a}\_2$ as well as the size of the lattice. The size is characterized by two vectors $\\vec{L}\_1$ and $\\vec{L}\_2$ and the lattice is placed on a torus:
$$\\hat{c}\_{\\vec{i} + \\vec{L}\_1 }  = \\hat{c}\_{\\vec{i} + \\vec{L}\_2 }  = \\hat{c}\_{\\vec{i}}$$
 The function call

    Call Make_Lattice( L1, L2, a1,  a2, Latt )

will generate the lattice `Latt` of type `Lattice`. Note that the structure of the unit cell has to be provided by the user. The reciprocal lattice vectors are defined by:
$$\\label{Latt.G.eq}
	\\vec{a}\_i  \\cdot \\vec{g}\_i = 2 \\pi \\delta\_{i,j},$$
 and the Brillouin zone corresponds to the Wigner Seitz cell of the lattice. With $\\vec{k} = \\sum\_{i} \\alpha\_i  \\vec{g}\_i $, the k-space quantization follows from:
$$\\begin{bmatrix}
	\\vec{L}\_1 \\cdot \\vec{g}\_1  &  \\vec{L}\_1 \\cdot \\vec{g}\_2  \\\\
	\\vec{L}\_2  \\cdot \\vec{g\_1} & \\vec{L}\_2 \\cdot  \\vec{g}\_2  
\\end{bmatrix}
\\begin{bmatrix}
   \\alpha\_1 \\\\
   \\alpha\_2
\\end{bmatrix}
=
2 \\pi 
\\begin{bmatrix}
   n \\\\
   m
\\end{bmatrix}$$
 such that
$$\\begin{aligned}
\\label{k.quant.eq}
     \\vec{k} =  n \\vec{b}\_1  + m \\vec{b}\_2 \\text{  with  }   & &   \\vec{b}\_1 = \\frac{2 \\pi}{ (\\vec{L}\_1 \\cdot \\vec{g}\_1)  (\\vec{L}\_2 \\cdot  \\vec{g}\_2 )  - (\\vec{L}\_1 \\cdot \\vec{g}\_2) (\\vec{L}\_2  \\cdot \\vec{g\_1} ) }   \\left\[  (\\vec{L}\_2 \\cdot  \\vec{g}\_2) \\vec{g}\_1 -   (\\vec{L}\_2  \\cdot \\vec{g\_1} ) \\vec{g}\_2 \\right\] \\text{   and  } \\nonumber \\\\ 
        & & \\vec{b}\_2 = \\frac{2 \\pi}{ (\\vec{L}\_1 \\cdot \\vec{g}\_1)  (\\vec{L}\_2 \\cdot  \\vec{g}\_2 )  - (\\vec{L}\_1 \\cdot \\vec{g}\_2) (\\vec{L}\_2  \\cdot \\vec{g\_1} ) }   
           \\left\[  (\\vec{L}\_1 \\cdot  \\vec{g}\_1) \\vec{g}\_2 -   (\\vec{L}\_1  \\cdot \\vec{g\_2} ) \\vec{g}\_1 \\right\] \\end{aligned}$$

| Variable                   | Type    | Description                                                                                                                           |
|:---------------------------|:--------|:--------------------------------------------------------------------------------------------------------------------------------------|
|                            | Real    | Unit vectors $\\vec{a}\_1$, $\\vec{a}\_2$                                                                                             |
|                            | Real    | Vectors $\\vec{L}\_1$, $\\vec{L}\_2$ that define the topology of the lattice.                                                         |
|                            |         | Tilted lattices are thereby possible to implement.                                                                                    |
| `Latt%N`                   | Integer | Number of lattice points, *N*<sub>unit cell</sub>                                                                                     |
| `Latt%list`                | Integer | maps each lattice point *i* = 1, ⋯, *N*<sub>unit cell</sub> to a real space vector                                                    |
|                            |         | denoting the position of the unit cell:                                                                                               |
|                            |         | $\\vec{R}\_i$ = `list(i,1)` $\\vec{a}\_1$ + `list(i,2)` $\\vec{a}\_2$ $  \\equiv i\_1  \\vec{a}\_1 + i\_2  \\vec{a}\_2 $              |
| `Latt%invlist`             | Integer | `Invlist`(*i*<sub>1</sub>, *i*<sub>2</sub>)=*i*                                                                                       |
| `Latt%nnlist`              | Integer | *j* = `nnlist`(*i*, *n*<sub>1</sub>, *n*<sub>2</sub>), *n*<sub>1</sub>, *n*<sub>2</sub> ∈ \[ − 1, 1\]                                 |
|                            |         | $\\vec{R}\_j = \\vec{R}\_i + n\_1 \\vec{a}\_1  + n\_2 \\vec{a}\_2 $                                                                   |
| `Latt%imj`                 | Integer | $ \\vec{R}\_{imj(i,j)}  =  \\vec{R}\_i -  \\vec{R}\_j$. *i**m**j*, *i*, *j* ∈ 1, ⋯, *N*<sub>unit cell</sub>                           |
| `Latt%BZ1_p`, `Latt%BZ2_p` | Real    | Reciprocal space vectors $\\vec{g}\_i$ (See Eq. \[Latt.G.eq\])                                                                        |
| `Latt%b1_p`, `Latt%b1_p`   | Real    | k-quantization (See Eq. \[k.quant.eq\])                                                                                               |
| `Latt%listk`               | Integer | maps each reciprocal lattice point *k* = 1, ⋯, *N*<sub>unit cell</sub>                                                                |
|                            |         | to a reciprocal space vector                                                                                                          |
|                            |         | $\\vec{k}\_k= \\texttt{listk(k,1)} \\vec{b}\_1 +  \\texttt{listk(k,2)} \\vec{b}\_2  \\equiv k\_1  \\vec{b}\_1 +   k\_2  \\vec{b}\_2 $ |
| `Latt%invlistk`            | Integer | `Invlistk`(*k*<sub>1</sub>, *k*<sub>2</sub>)=*k*                                                                                      |
| `Latt%b1_perp_p`,          |         |                                                                                                                                       |
| `Latt%b2_perp_p`           | Real    | Orthonormal vectors to $\\vec{b}\_i$. For internal use.                                                                               |

The module equally handles the Fourier transformation. For example the subroutine carries out the transformation:
$$S(\\vec{k}, :,:,:) =  \\frac{1}{N\_{unit \\,cell}}  \\sum\_{\\vec{i},\\vec{j}}   e^{-i \\vec{k} \\cdot \\left( \\vec{i}-\\vec{j} \\right)} S(\\vec{i}  - \\vec{j}, :,:,:)$$
 and the inverse Fourier transform
$$S(\\vec{r}, :,:,:) =  \\frac{1}{N\_{unit \\,cell}}  \\sum\_{\\vec{k} \\in BZ }   e^{ i \\vec{k} \\cdot \\vec{r}} S(\\vec{k}, :,:,:).$$
 In the above, the unspecified dimensions of structure factor can refer to imaginary time, and orbital indices.

The observable types `Obser_Vec` and `Obser_Latt`
-------------------------------------------------

Our definition of the model includes observables \[Eq. (\[eqn:obs\_rw\])\]. We have defined two observable types: `Obser_vec` for a array of scalar observables such as the energy and `Obser_Latt` for correlation functions that have the lattice symmetry. In the latter case, translation symmetry can be used to provide improved estimators and to reduce the size of the I/O. We also obtain improved estimators by taking measurements in the imaginary-time interval `[LOBS_ST,LOBS_EN]` (see the parameter file in Sec. \[sec:input\]) thereby exploiting the invariance under translation in imaginary time. Note that the translation symmetries in space and in time are broken for a given configuration *C* but restored by the Monte Carlo sampling. In general, the user will define bins, each bins having a given amount of sweeps. Within a sweep we run sequentially trough the HS and Ising fields from time slice 1 to *L*<sub>Trotter</sub> and back. The results of each bin is written in a file and analyzed at the end of the run.

To accomplish the reweighting of observables (see Sec. \[sec:reweight\]), for each configuration the measurement of an observable has to be multiplied by the factors `ZS` and `ZP`:
$$\\begin{aligned}
\\texttt{ZS} &=& \\text{sign}(C)\\;\\\\
\\texttt{ZP} &=& \\frac{e^{-S(C)}} {\\Re \\left\[e^{-S(C)} \\right\]}\\;,\\end{aligned}$$
 They are computed from the Monte Carlo phase of a configuration,
$$\\label{eqn:phase}
	\\texttt{phase}   =   \\frac{e^{-S(C)}}{ \\left| e^{-S(C) }\\right| }\\;,$$
 which is provided by the main program.

Note that each observable structure also includes the average sign \[Eq. (\[eqn:sign\_rw\])\].

### Scalar observables

This data type is described in Table \[table:Obser\_vec\] and is useful to compute an array of scalar observables. Consider a variable `Obs` of type `Obser_vec`. At the beginning of each bin, a call to `Obser_Vec_Init` in the module `observables_mod.f90` will set `Obs%N=0`, `Obs%Ave_sign =0` and `Obs%Obs_vec(:)=0`. Each time the main program calls the routine `Obser` in the `Hamiltonian` module, the counter `Obs%N` is incremented by unity, the sign (see Eq. \[Sign.eq\]) is cumulated in the variable `Obs%Ave_sign`, and the desired the observables (multiplied by the sign and $\\frac{e^{-S(C)}} {\\Re \\left\[e^{-S(C)} \\right\]}$, see Sec. \[Observables.General\]) are cumulated in the vector `Obs%Obs_vec`.

|                  |        |                                                           |                                                                                                                                |
|:-----------------|:-------|:----------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------|
| Variable         | Type   | Description                                               | Contribution of                                                                                                                |
|                  |        |                                                           | configuration *C*                                                                                                              |
| `Obs%N`          | Int.   | Number of measurements                                    |                                                                                                                                |
| `Obs%Ave_sign`   | Real   | Cumulated sign \[Eq. (\[eqn:sign\_rw\])\]                 | sign(*C*)                                                                                                                      |
| `Obs%Obs_vec(:)` | Compl. | Cumulated vector of observables \[Eq. (\[eqn:obs\_rw\])\] | $ \\langle \\langle \\hat{O}(:) \\rangle \\rangle\_{C}\\frac{e^{-S(C)}} {\\Re \\left\[e^{-S(C)} \\right\]} \\text{ sign }(C) $ |
| `Obs%File_Vec`   | Char.  | Name of output file                                       |                                                                                                                                |

At the end of the bin, a call to `Print_bin_Vec` in module `observables_mod.f90` will append the result of the bin in the file `File_Vec`\_scal. Note that this subroutine will automatically append the suffix \_scal to the the filename `File_Vec`. This suffix is important to allow automatic analysis of the data at the end of the run.

###  Equal time and time-displaced correlation functions

|                                            |        |                                                |                                                                                                                                                                                         |
|:-------------------------------------------|:-------|:-----------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Variable                                   | Type   | Description                                    | Contribution of                                                                                                                                                                         |
|                                            |        |                                                | configuration *C*                                                                                                                                                                       |
| `Obs%N`                                    | Int.   | Number of measurements                         |                                                                                                                                                                                         |
| `Obs%Ave_sign`                             | Real   | Cumulated sign \[Eq. (\[eqn:sign\_rw\])\]      | sign(*C*)                                                                                                                                                                               |
| `Obs%Obs_latt`                             | Compl. | Cumul. correl. fct. \[Eq. (\[eqn:obs\_rw\])\]  | $ \\langle \\langle \\hat{O}\_{\\vec{i},\\alpha} (\\tau) \\hat{O}\_{\\vec{j},\\beta} \\rangle \\rangle\_{C} \\; \\frac{e^{-S(C)}} {\\Re \\left\[e^{-S(C)} \\right\]}  \\text{sign}(C) $ |
| $(\\vec{i}-\\vec{j},\\tau,\\alpha,\\beta)$ |        |                                                |                                                                                                                                                                                         |
| `Obs%Obs_latt0(\alpha)`                    | Compl. | Cumul. expect. value \[Eq. (\[eqn:obs\_rw\])\] | $ \\langle \\langle \\hat{O}\_{\\vec{i},\\alpha} \\rangle \\rangle\_{C}\\frac{e^{-S(C)}} {\\Re \\left\[e^{-S(C)} \\right\]}  \\text{ sign }(C) $                                        |
| `Obs%File_Latt`                            | Char.  | Name of output file                            |                                                                                                                                                                                         |

This data type (see Table \[table:Obser\_latt\]) is useful so as to deal with imaginary time displaced as well as equal time correlation functions of the form:
$$\\label{eqn:s}
	S\_{\\alpha,\\beta}(\\vec{k},\\tau) =   \\frac{1}{N\_{\\text{unit cell}}} \\sum\_{\\vec{i},\\vec{j}}  e^{- \\vec{k} \\cdot \\left( \\vec{i}-\\vec{j}\\right) } \\left( \\langle \\hat{O}\_{\\vec{i},\\alpha} (\\tau) \\hat{O}\_{\\vec{j},\\beta} \\rangle  - 
	  \\langle \\hat{O}\_{\\vec{i},\\alpha} \\rangle \\langle   \\hat{O}\_{\\vec{i},\\beta}  \\rangle \\right).$$
 Here, translation symmetry of the Bravais lattice is explicitly taken into account. The correlation function splits in a correlated part $S\_{\\alpha,\\beta}^{\\mathrm{(corr)}}(\\vec{k},\\tau)$ and a background part $S\_{\\alpha,\\beta}^{\\mathrm{(back)}}(\\vec{k})$:
$$\\begin{aligned}
  S\_{\\alpha,\\beta}^{\\mathrm{(corr)}}(\\vec{k},\\tau)
  &=&
   \\frac{1}{N\_{\\text{unit cell}}} \\sum\_{\\vec{i},\\vec{j}}  e^{- i\\vec{k} \\cdot \\left( \\vec{i}-\\vec{j}\\right) }  \\langle \\hat{O}\_{\\vec{i},\\alpha} (\\tau) \\hat{O}\_{\\vec{j},\\beta} \\rangle\\label{eqn:s\_corr}\\;,\\\\
         S\_{\\alpha,\\beta}^{\\mathrm{(back)}}(\\vec{k})
  &=&
   \\frac{1}{N\_{\\text{unit cell}}} \\sum\_{\\vec{i},\\vec{j}}  e^{- i\\vec{k} \\cdot \\left( \\vec{i}-\\vec{j}\\right) }  \\langle \\hat{O}\_{\\vec{i},\\alpha} (\\tau)\\rangle \\langle \\hat{O}\_{\\vec{j},\\beta} \\rangle\\nonumber\\\\
  &=& 
  N\_{\\text{unit cell}}\\, \\langle \\hat{O}\_{\\alpha} \\rangle \\langle \\hat{O}\_{\\beta} \\rangle \\, \\delta(\\vec{k})\\label{eqn:s\_back}\\;,\\end{aligned}$$
 where translation invariance in space and time has been exploited to obtain the last line. The background part depends only on the expectation value $\\langle \\hat{O}\_{\\alpha} \\rangle$, for which we use the following estimator
$$\\label{eqn:o}
\\langle \\hat{O}\_{\\alpha} \\rangle \\equiv \\frac{1}{N\_{\\text{unit\\,cell}}} \\sum\\limits\_{\\vec{i}} \\langle \\hat{O}\_{\\vec{i},\\alpha} \\rangle\\;.$$

Consider a variable `Obs` of type `Obser_latt`. At the beginning of each bin a call to `Obser_Latt_Init` in the module `observables_mod.f90` will initialize the elements of `Obs` to zero. Each time the main program calls the `Obser` or `ObserT` routines one cumulates $ \\langle \\langle \\hat{O}\_{\\vec{i},\\alpha} (\\tau) \\hat{O}\_{\\vec{j},\\beta} \\rangle \\rangle\_{C} \\; \\frac{e^{-S(C)}} {\\Re \\left\[e^{-S(C)} \\right\]}  \\text{sign}(C) $ in `Obs%Obs_latt(\vec{i}-\vec{j},\tau,\alpha,\beta)` and $ \\langle \\langle \\hat{O}\_{\\vec{i},\\alpha} \\rangle \\rangle\_{C}\\frac{e^{-S(C)}} {\\Re \\left\[e^{-S(C)} \\right\]}  \\text{ sign }(C) $ in `Obs%Obs_latt0(\alpha)`. At the end of each bin, a call to `Print_bin_Latt` in the module `observables_mod.f90` will append the result of the bin in the specified file `Obs%File_Latt`. Note that the routine `Print_bin_Latt` carries out the Fourier transformation and prints the results in k-space. We have adopted the following name convention. For equal time observables , that is the second dimension of the array `Obs%Obs_latt(\vec{i}-\vec{j},\tau,\alpha,\beta)` is equal to unity, the routine `Print_bin_Latt` attaches the suffix \_eq to `Obs%File_Latt`. For time displaced correlation functions we use the suffix \_tau.

Updating schemes
================

The program allows for different types of updating schemes. Given a configuration *C* we propose a new one, *C*′, with probability *T*<sub>0</sub>(*C* → *C*′) and accept it according to
$$P(C \\rightarrow C') =  \\text{min}  \\left( 1, \\frac{T\_0(C' \\rightarrow C) W(C')}{T\_0(C \\rightarrow C') W(C)} \\right)$$
 so as to guarantee the stationarity condition. Here, *W*(*C*)=|ℜ\[*e*<sup>−*S*(*C*)</sup>\]|.

<span>@ l l l @</span> Variable & Type & Description
`Propose_S0` & Logical & If true, proposes local moves according to the probability *e*<sup>−*S*<sub>0</sub></sup>
`Global_moves` & Logical & If true, allows for global moves.
`N_Global ` & Integer & Number of global moves per sweep of single spin flips.

The Default: sequential single spin flips
-----------------------------------------

The default updating scheme is a sequential single spin flip algorithm. Consider the Ising spin *s*<sub>*i*, *τ*</sub>, we will flip it with probability one such that for this local move the proposal matrix is symmetric. If we are considering the Hubbard-Stratonovich field *l*<sub>*i*, *τ*</sub> we will propose with probability 1/3 one of the other three possible fields. Again, for this local move, the proposal matrix is symmetric. Hence in both cases we will accept or reject the move according to
$$P(C \\rightarrow C') =  \\text{min}  \\left( 1, \\frac{ W(C')}{W(C)} \\right)$$
 In is worth noting that this type of sequential spin flip updating does not satisfy detailed balance but the more fundamental stationarity condition.

Sampling of *e*<sup>−*S*<sub>0</sub></sup>
------------------------------------------

Consider an Ising spin at space time *i*, *τ* and the configuration *C*. Flipping this spin will generate the configuration *C*′ and we will propose the move according to
$$T\_0(C \\rightarrow C')  =  \\frac{e^{-S\_0(C')}}{ e^{-S\_0(C')} + e^{-S\_0(C)} }   = 1 - \\frac{1}{1 +  e^{-S\_0(C')} /e^{-S\_0(C)}}$$
 Note that the function `S0` in the `Hamitonian_example.f90` module computes precisely the ratio *e*<sup>−*S*<sub>0</sub>(*C*′)</sup>/*e*<sup>−*S*<sub>0</sub>(*C*)</sup> so that *T*<sub>0</sub>(*C* → *C*′) does not require any further programming. Thereby one will accept the proposed move with the probability:
$$P(C \\rightarrow C') =  \\text{min}  \\left( 1,  \\frac{e^{-S\_0(C)}   W(C')}{ e^{-S\_0(C)'} W(C)} \\right).$$
 With Eq. \[eqn:partition\_2\] one sees that the bare action *S*<sub>0</sub>(*C*) determining the dynamics of the Ising spin in the absence of coupling to the fermions does not enter the Metropolis acceptance rejection step. This sampling scheme is used if the logical variable `Propose_S0` is switched on.

Global updates
--------------

The code equally allows for global updates. The user will have to provide two other functions in the module `Hamiltonian_Examples.f90`.

The subroutine `Global_move(T0_Proposal_ratio,nsigma_old)` proposes a global move.
`nsigma_old(M_V+ M_I, Ltrot)` is a two dimensional array containing the full configuration *C*. On output, the new configuration, C’,– determined by the user – is to be stored in the array `nsigma(M_V+ M_I, Ltrot)`. `nsigma(M_V+ M_I, Ltrot)` is a global variable declared in the module, `Hamiltonian`. Equally, on output, the variable `T0_Proposal_ratio` contains the proposal ratio
$$\\frac{T\_0(C' \\rightarrow C)}{T\_0(C \\rightarrow C') }$$
 Since we allow for a stochastic generation of the global move, it may very well be that no change is proposed. In this case, `T0_Proposal_ratio` takes the value 0 upon exit, and `nsigma=nsigma_old`.

To compute the acceptance rejection ratio, the user will equally have to provide the function
`Delta_S0_global(Nsigma_old)` that computes the ratio *e*<sup>−*S*<sub>0</sub>(*C*′)</sup>/*e*<sup>−*S*<sub>0</sub>(*C*)</sup>. Again the configuration *C*′ is given by the array `nsigma(M_V+ M_I, Ltrot)` which is a global variable declared in the module, `Hamiltonian`.

Note that global updates are expensive since they require a complete recalculation of the weight. We thereby allow the user to set a variable `N_Global` that allows to determine how many global updates per sweeps will be carried out.

File structure
==============

| Directory | Description                                 |
|:----------|:--------------------------------------------|
|           | Main program and subroutines                |
|           | Collection of mathematical routines         |
|           | Routines for error analysis                 |
|           | Example simulations for Hubbard-type models |
|           | Parameter files and scripts                 |
|           | Documentation of the QMC code.              |

The code package consists of the program directories , and . The sample simulations corresponding to the walkthroughs of Sec. \[sec:walk1\] - \[sec:walk2\] are included in . The package content is summarized in Table \[table:files\].

Input files
-----------

| File | Description                                                                    |
|:-----|:-------------------------------------------------------------------------------|
|      | Sets the parameters for lattice, model, QMC process, and the error analysis.   |
|      | List of integer numbers to initialize the random number generator and          |
|      | to start a simulation from scratch.                                            |
|      | Input files for the HS and Ising configuration, used to continue a simulation. |

The input files are listed in Table \[table:input\]. The parameter file has the following form, using as an example the *S**U*(2)-symmetric Hubbard model on a square lattice (see Sec. \[sec:walk1\] for a detailed walkthrough):

    ===============================================================================
    !  Variables for the Hubb program
    !-------------------------------------------------------------------------------
    &VAR_lattice
    L1 = 4                    ! Length in direction a_1
    L2 = 4                    ! Length in direction a_2
    Lattice_type = "Square"	  ! a_1 = (1,0),  a_2=(0,1),  Norb=1, N_coord=2
    !Lattice_type ="Honeycomb"! a_1 = (1,0),  a_2 =(1/2,sqrt(3)/2), Norb=2, N_coord=3
    Model = "Hubbard_SU2"     ! Sets Nf=1, N_sun=2. HS field couples to the density
    !Model = "Hubbard_Mz"     ! Sets Nf=2, N_sun=1. HS field couples to the 
                              ! z-component of magnetization.  
    !Model="Hubbard_SU2_Ising"! Sets Nf_1, N_sun=2 and runs only for the square lattice
                              ! Hubbard model  coupled to transverse Ising field
    /
    &VAR_Hubbard              ! Variables for the Hubbard model
    ham_T   = 1.D0            ! Hopping parameter
    ham_chem= 0.D0            ! chemical potential
    ham_U   = 4.D0            ! Hubbard interaction
    Beta    = 5.D0            ! inverse temperature
    dtau    = 0.1D0           ! Thereby Ltrot=Beta/dtau
    /

    &VAR_Ising                ! Model parameters for the Ising code
    Ham_xi = 1.d0             ! Only needed if Model="Hubbard_SU2_Ising"
    Ham_J  = 0.2d0
    Ham_h  = 2.d0
    /

    &VAR_QMC                  ! Variables for the QMC run
    Nwrap   = 10              ! Stabilization. Green functions will be computed from scratch 
                              ! after each time interval  Nwrap*Dtau
    NSweep  = 500             ! Number of sweeps
    NBin    = 2               ! Number of bins
    Ltau    = 1               ! 1 for calculation of time displaced Green functions. 0 otherwise
    LOBS_ST = 1               ! Start measurements at time slice LOBS_ST
    LOBS_EN = 50              ! End   measurements at time slice LOBS_EN
    CPU_MAX = 0.1             ! Code will stop after CPU_MAX hours. 
                              ! If not specified, code will stop after Nbin bins.
    /                          
    &VAR_errors               ! Variables for analysis programs
    n_skip  = 1               ! Number of bins that will be skipped. 
    N_rebin = 1               ! Rebinning  
    N_Cov   = 0               ! If set to 1 covariance will be computed
                              ! for unequal time correlation functions.                   
    /            

Output files
------------

| File | Description                                                                      |
|:-----|:---------------------------------------------------------------------------------|
|      | After completion of the simulation, this file documents parameters of the model, |
|      | the QMC run and simulation metrics (precision, acceptance rate, CPU time).       |
|      | Results of equal-time measurements of scalar observables.                        |
|      | The placeholder stands for the observables , and .                               |
|      | Results of equal-time and time-displaced measurements of correlation functions.  |
|      | The placeholder stands for , and .                                               |
|      | Output files for the HS and Ising configuration.                                 |

The output of the measured data is organized in bins. One bin corresponds to the geometric average over a fixed number of individual measurements which depends on the chosen measurement interval on the imaginary time axis and on the number of Monte Carlo sweeps. If the user run a parallelized version of the code, the average also extends over the number of MPI threads. The standard output files are listed in Table \[table:output\].

The formatting of the output for a single bin depends on the observable type: or .

-   Observables of type : For each additional bin, a single new line is added to the output file. In case of an observable with components, the formatting is

        N_size + 1    <measured value, 1> ... <measured value, N_size>    <measured sign>

    The counter variable refers to the number of measurements per line, including the phase measurement. This format is required by the error analysis routine (see Sec. \[sec:analysis\]). Scalar observables like kinetic energy, potential energy, total energy and particle number are treated as a vector of size .

-   Observables of type : For each additional bin, a new data block is added to the output file. The block consists of the expectation values \[Eq. (\[eqn:o\])\] contributing to the background part \[Eq. (\[eqn:s\_back\])\] of the correlation function, and the correlated part \[Eq. (\[eqn:s\_corr\])\] of the correlation function. For imaginary-time displaced correlation functions, the formatting of the block follows this scheme:

    `<measured sign>  <N_orbital>  <N_unit_cell> <N_time_slices> <dtau>`
    `do alpha = 1, N_orbital`
    `    `$\\langle\\hat{O}\\sb{\\alpha}\\rangle\\ $
    `enddo`
    `do i = 1, N_unit_cell`
    `   <reciprocal lattice vector k(i)>`
    `   do tau = 1, N_time_slices`
    `      do alpha = 1, N_orbital`
    `         do beta = 1, N_orbital`
    `            `$\\langle\\,S\\sb{\\alpha,\\beta}\\sp{(corr)}(k(i),\\tau)\\rangle$
    `         enddo`
    `      enddo`
    `   enddo`
    `enddo`

    The same block structure is used for equal-time correlation functions, except for the entries and which are not present in the latter. Using this structure for the bins as input, the full correlation function $S\_{\\alpha,\\beta}(\\vec{k},\\tau)$ \[Eq. (\[eqn:s\])\] is then calculated by calling the error analysis routine (see Sec. \[sec:analysis\])

### The `info` file and stabilization

The finite temperature auxiliary field QMC algorithm is known to be numerically unstable. The origin the numerical instabilities arises from the imaginary time propagation which invariably leads to exponentially small and exponentially large scales. Numerical stabilization of the code is delicate and has been pioneered in Ref.  for the finite temperature algorithm and in Refs.  for the zero temperature projective algorithm. As shown in Ref.  scales can be omitted in the ground state algorithm – thus rendering it very stable – but have to be taken into account in the finite temperature code. Apart from runtime information, the file `info` contains important information concerning the stability of the code. For example, in the directory an example simulation of the 4 × 4 Hubbard model at *U*/*t* = 4 and *β**t* = 10, the `info` file contains the lines

`Precision Green  Mean, Max :    1.2918865817224671E-014   4.0983018995027644E-011`
`Precision Phase, Max       :    5.0272908791449966E-012`
`Precision tau    Mean, Max :    8.4596701790588625E-015   3.5033530012121281E-011`

showing the mean and maximum difference between the <span>*wrapped* </span> and from scratched computed equal and time displaced Green functions . A stable code should produce results where the mean difference is smaller than the stochastic error. The above example shows a very stable simulation since the Green function is of order 1. Numerical stabilization is delicate and there is no guarantee that it will work for all models. For example switching to a HS field which couples to the z-component of the magnetization will yield (see directory ):

`Precision Green  Mean, Max :    5.0823874429126405E-011   5.8621144596315844E-006`
`Precision Phase, Max       :    0.0000000000000000     `
`Precision tau    Mean, Max :    1.5929357848647394E-011   1.0985132530727526E-005 `

This is still an excellent precision but nevertheless a couple of order of magnitudes less precise than a HS decomposition coupling to the charge. If the numerical stabilization turns out to be bad, one option is to reduce the value of the parameter `Nwrap` in the parameter file. For performing the stabilization of the involved matrix multiplications we rely on routines from lapack. Hence it is very likely that your results may change significantly if you switch libraries. In order to offer a simple baseline to which people can quickly switch if they want to see whether their results depend on the library used for linear algebra routines we have included parts of the lapack-3.6.1 reference implementation from <http://www.netlib.org/lapack/>. You can switch to the QR decomposition related routines from the lapack reference implementation by including the switch `-DQRREF` into their PROGRAMMCONFIGURATION string. To use these routines you need to link against a lapack library that implements at least the lapack-3.4.0 interface. [2]

To provide further flexibility, we have incorporated different stabilization schemes. Our default strategy is quick and generically works well but we have encountered some models where it fails. If this applies to your model, you can use the switch `-DSTAB1` in the `set_env.sh` file and recompile the code.

Scripts
-------

| Script | Description                                            | Section          |
|:-------|:-------------------------------------------------------|:-----------------|
|        | Sets the environment variables for the compiler        |                  |
|        | and the libraries.                                     | \[sec:running\]  |
|        | Copies the output configurations of HS and Ising spins |                  |
|        | to the respective input files.                         | \[sec:running\]  |
|        | Starts the error analysis.                             | \[sec:analysis\] |

Walkthrough: the *S**U*(2)-Hubbard model on a square lattice
============================================================

To implement a Hamiltonian, the user has to provide a module which specifies the lattice, the model, as well as the observables he/she wishes to compute. In this section, we describe the module
`Hamiltonian_Examples.f90` which contains an implementation of the Hubbard model on the square lattice. A sample run for this model can be found in .

The Hamiltonian reads
$$\\label{eqn\_hubbard\_sun}
\\mathcal{H}= 
\\sum\\limits\_{\\sigma=1}^{2} 
\\sum\\limits\_{x,y =1 }^{N\_{\\text{unit cell}}} 
  c^{\\dagger}\_{x \\sigma} T\_{x,y}c^{\\phantom\\dagger}\_{y \\sigma} 
+ \\frac{U}{2}\\sum\\limits\_{x}\\left\[
\\sum\\limits\_{\\sigma=1}^{2}
\\left(  c^{\\dagger}\_{x \\sigma} c^{\\phantom\\dagger}\_{x \\sigma}  -1/2 \\right) \\right\]^{2}\\;.$$
 We can make contact with the general form of the Hamiltonian by setting: *N*<sub>*f**l*</sub> = 1, *N*<sub>*c**o**l*</sub> ≡ `N_SUN` = 2, *M*<sub>*T*</sub> = 1, *T*<sub>*x**y*</sub><sup>(*k**s*)</sup> = *T*<sub>*x*, *y*</sub>, *M*<sub>*V*</sub> = *N*<sub>unit cell</sub>, $U\_{k}       =   -\\frac{U}{2}$, *V*<sub>*x**y*</sub><sup>(*k**s*)</sup> = *δ*<sub>*x*, *y*</sub>*δ*<sub>*x*, *k*</sub>, $\\alpha\_{ks}   = - \\frac{1}{2}  $ and *M*<sub>*I*</sub> = 0.

Setting the Hamiltonian: `Ham_set` 
-----------------------------------

The main program will call the subroutine `Ham_set` in the module `Hamiltonian_Hub.f90`. The latter subroutine defines the public variables


      Type (Operator), dimension(:,:), allocatable  :: Op_V 
      Type (Operator), dimension(:,:), allocatable  :: Op_T
      Integer, allocatable :: nsigma(:,:)
      Integer              :: Ndim,  N_FL,  N_SUN,  Ltrot

which specify the model. The array `nsigma` contains the HS field. The routine `Ham_set` will first read the parameter file, then set the lattice, `Call Ham_latt`, set the hopping `Call Ham_hop` and set the interaction `call Ham_V`. The parameters are read in from the file `parameters`, see Sec. \[sec:input\].

### The lattice: `Call Ham_latt` 

The choice `Lattice_type = Square` sets $\\vec{a}\_1 =  (1,0) $ and $\\vec{a}\_2 =  (0,1) $ and for an *L*<sub>1</sub> × *L*<sub>2</sub> lattice $\\vec{L}\_1 = L\_1 \\vec{a}\_1$ and $\\vec{L}\_2 = L\_2 \\vec{a}\_2$. The call to ` Call Make_Lattice( L1, L2, a1, a2, Latt)` will generate the lattice `Latt` of type `Lattice` . For the Hubbard model on the square lattice, the number of orbitals per unit cell is given by `NORB=1` such that *N*<sub>*d**i**m*</sub> ≡ *N*<sub>unit cell</sub> ⋅ `NORB` = `Latt%N`.

### The hopping term: `Call Ham_hop`

The hopping matrix is implemented as follows. We allocate an array of dimension 1 × 1 of type operator called `Op_T` and set the dimension for the hopping matrix to *N* = *N*<sub>*d**i**m*</sub>. One allocates and initializes this type by a single call to the subroutine `Op_make`:


    call Op_make(Op_T(1,1),Ndim)

Since the hopping does not break down into small blocks, we have ${\\bm P}=\\mathds{1}$ and


    Do i= 1,Ndim
      Op_T(1,1)%P(i) = i
    Enddo

We set the hopping matrix with


    DO I = 1, Latt%N
       Ix = Latt%nnlist(I,1,0)
       Iy = Latt%nnlist(I,0,1)
       Op_T(1,1)%O(I  ,Ix) = cmplx(-Ham_T,   0.d0,kind(0.D0))
       Op_T(1,1)%O(Ix,I  ) = cmplx(-Ham_T,   0.d0,kind(0.D0))
       Op_T(1,1)%O(I  ,Iy) = cmplx(-Ham_T,   0.d0,kind(0.D0))
       Op_T(1,1)%O(Iy, I ) = cmplx(-Ham_T,   0.d0,kind(0.D0))
       Op_T(1,1)%O(I  ,I ) = cmplx(-Ham_chem,0.d0,kind(0.D0))
    ENDDO

Here, the integer function ` j= Latt%nnlist(I,n,m)` is defined in the lattice module and returns the index of the lattice site $ \\vec{I} +  n \\vec{a}\_1 +  m \\vec{a}\_2$. Note that periodic boundary conditions are already taken into account. The hopping parameter, `Ham_T` as well as the chemical potential `Ham_chem` are read from the parameter file. Note that although a checkerboard decomposition is not used here, it can be implemented by considering a larger number of sparse hopping matrices

### The interaction term: `Call Ham_V`

To implement this interaction, we allocate an array of `Operator` type. The array is called `Op_V` and has dimensions *N*<sub>*d**i**m*</sub> × *N*<sub>*f**l*</sub> = *N*<sub>*d**i**m*</sub> × 1. We set the dimension for the interaction term to *N* = 1, and allocate and initialize this array of type `Operator` by repeatedly calling the subroutine `Op_make`:


    do i  = 1,Ndim
       call Op_make(Op_V(i,1),1)
    enddo

For each lattice site *i*, the matrices ${\\bm P}$ are of dimension 1 × *N*<sub>*d**i**m*</sub> and have only one non-vanishing entry. Thereby we can set:


    Do i = 1,Ndim
       Op_V(i,1)%P(1)   = i
       Op_V(i,1)%O(1,1) = cmplx(1.d0,0.d0, kind(0.D0))
       Op_V(i,1)%g      = sqrt(cmplx(-dtau*ham_U/(dble(N_SUN)),0.D0,kind(0.D0)))
       Op_V(i,1)%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
       Op_V(i,1)%type   = 2
    Enddo

so as to completely define the interaction term.

Observables
-----------

At this point, all the information for the simulation to start has been provided. The code will sequentially go through the operator list `Op_V` and update the fields. Between time slices `LOBS_ST` and `LOBS_EN` the main program will call the routine `Obser(GR,Phase,Ntau)` which is provided by the user and handles equal time correlation functions. If `Ltau=1` the the main program will call the routine `ObserT(NT, GT0,G0T,G00,GTT, PHASE) ` which is again

The user will have to implement the observables he/she wants to compute. Here we will describe how to proceed.

### Allocating space for the observables: `Call Alloc_obs(Ltau) `

For four scalar or vector observables, the user will have to declare the following:


    Allocate ( Obs_scal(4) )
    Do I = 1,Size(Obs_scal,1)
       select case (I)
       case (1)
          N = 2;  Filename ="Kin"
       case (2)
          N = 1;  Filename ="Pot"
       case (3)
          N = 1;  Filename ="Part"
       case (4)
          N = 1,  Filename ="Ener"
       case default
          Write(6,*) ' Error in Alloc_obs '  
       end select
       Call Obser_Vec_make(Obs_scal(I),N,Filename)
    enddo

Here, `Obs_scal(1)` contains a vector of two observables so as to account for the x -and -y components of the kinetic energy for example.

For equal time correlation functions we allocate `Obs_eq` of type `Obser_Latt`. Here we include the calculation of spin-spin and density-density correlation functions alongside equal time Green functions.


    Allocate ( Obs_eq(4) )
    Do I = 1,Size(Obs_eq,1)
       select case (I)
       case (1)
          Ns = Latt%N; No = Norb;  Filename ="Green"
       case (2)
          Ns = Latt%N; No = Norb;  Filename ="SpinZ"
       case (3)
          Ns = Latt%N; No = Norb;  Filename ="SpinXY"
       case (4)
          Ns = Latt%N; No = Norb;  Filename ="Den"
       case default
          Write(6,*) ' Error in Alloc_obs '  
       end select
       Nt = 1
       Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
    enddo
     

For the Hubbard model `Norb = 1` and for equal time correlation functions `Nt = 1`. If `Ltau = 1` then the code will allocate space for time displaced quantities. The same structure as for equal time correlation functions will be used albeit with `Nt = Ltrot + 1`. At the beginning of each bin, the main program will set the bin observables to zero by calling the routine `Init_obs(Ltau)`. The user does not have to edit this routine.

### Measuring equal time observables: `Obser(GR,Phase,Ntau)`

The equal time green function,
$$\\texttt{GR(x,y},\\sigma{\\texttt)}  = \\langle c^{\\phantom{\\dagger}}\_{x,\\sigma} c^{\\dagger}\_{y,\\sigma}  \\rangle,$$
 the phase factor `phase` \[Eq. (\[eqn:phase\])\] and time slice `Ntau` is provided by the main program.

Here, *x* and *y* label unit-cell as well as the orbital within the unit cell. For the Hubbard model described here, *x* corresponds to the unit cell. The Green function does not depend on the color index, and is diagonal in flavor. For the SU(2)-symmetric implementation there is only one flavor, *σ* = 1 and the Green function is independent on the spin index. This renders the calculation of the observables particularly easy.

An explicit calculation of the potential energy $ \\langle U \\sum\_{\\vec{i}}  \\hat{n}\_{\\vec{i},\\uparrow}   \\hat{n}\_{\\vec{i},\\downarrow}  \\rangle $ reads


    Obs_scal(2)%N     = Obs_scal(2)%N + 1
    Obs_scal(2)%Ave_sign = Obs_scal(2)%Ave_sign + Real(ZS,kind(0.d0))
    Do i = 1,Ndim
       Obs_scal(2)%Obs_vec(1) = Obs_scal(2)%Obs_vec(1) + (1-GR(i,i,1))**2 * Ham_U * ZS * ZP
    Enddo

Here `ZS` =  sign(*C*) \[see Eq. (\[Sign.eq\])\], $ \\texttt{ZP} =   \\frac{e^{-S(C)}} {\\Re \\left\[e^{-S(C)} \\right\]}   $ \[see Eq. (\[eqn:phase\])\] and `Ham_U` corresponds to the Hubbard *U* term.

Equal time correlations are also computed in this routine. As an explicit example, we consider the equal time density-density fluctuations:
$$\\langle n\_{\\vec{i},\\alpha}   n\_{\\vec{j},\\beta} \\rangle   -  \\langle n\_{\\vec{i},\\alpha} \\rangle  \\langle    n\_{\\vec{j},\\beta}  \\rangle$$
 For the calculation of such quantities, it is convenient to define:
`GRC(x,y,s)` = *δ*<sub>*x*, *y*</sub> − `GR(y,x,s) `
 such that `GRC(x,y,s)` corresponds to $ \\langle \\langle  \\hat{c}\_{x,s}^{\\dagger}\\hat{c}\_{y,s}^{\\phantom\\dagger} \\rangle \\rangle $.


    Obs_eq(4)%N     = Obs_eq(4)%N + 1       ! Even if it is redundant, each observable carries 
    Obs_eq(4)%Ave_sign = Obs_eq(4)%Ave_sign + Real(ZS,kind(0.d0))  ! its own counter and sign.
    Do I1 = 1,Ndim
       I    = List(I1,1)                    ! = I  For the Hubbard model  on the square
       no_I = List(I1,2)                    ! = 1  lattice there is one orbital per unit-cell. 
       Do J1 = 1,Ndim
          J    = List(J1,1)
          no_J = List(J1,2)
          imj = latt%imj(I,J)
          Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) + &
                         &     (    GRC(I1,I1,1) * GRC(J1,J1,1) * N_SUN * N_SUN      + &
                         &          GRC(I1,J1,1) * GR(I1,J1,1) * N_SUN   ) * ZP * ZS 
       Enddo
       Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +   GRC(I1,I1,1) * N_SUN * ZP * ZS
    Enddo

At the end of each bin the main program will call the routine ` Pr_obs(LTAU)`. This routine will append the result of the bins in the specified file, with appropriate suffix.

### Measuring time-displaced observables: `ObserT(NT, GT0,G0T,G00,GTT, PHASE) `

This subroutine is called by the main program at the beginning of each sweep, provided that `LTAU` is set to unity. `NT` runs from `0` to `Ltrot` and denotes the imaginary time difference. For a given time displacement, the main program provides:
$$\\begin{aligned}
\\texttt{GT0(x,y,s) }  &=&   \\phantom{+} \\langle \\langle \\hat{c}^{\\phantom\\dagger}\_{x,s} (Nt \\Delta \\tau)   \\hat{c}^{\\dagger}\_{y,s} (0)   \\rangle \\rangle = \\langle \\langle {\\cal T} \\hat{c}^{\\phantom\\dagger}\_{x,s} (Nt \\Delta \\tau)   \\hat{c}^{\\dagger}\_{y,s} (0)   \\rangle \\rangle  \\nonumber \\\\
\\texttt{G0T(x,y,s) }   &=&  -   \\langle \\langle   \\hat{c}^{\\dagger}\_{y,s} (Nt \\Delta \\tau)    \\hat{c}^{\\phantom\\dagger}\_{x,s} (0)    \\rangle \\rangle =
    \\langle \\langle {\\cal T} \\hat{c}^{\\phantom\\dagger}\_{x,s} (0)    \\hat{c}^{\\dagger}\_{y,s} (Nt \\Delta \\tau)   \\rangle \\rangle  \\nonumber  \\\\
  \\texttt{G00(x,y,s) }  &=&    \\phantom{+} \\langle \\langle \\hat{c}^{\\phantom\\dagger}\_{x,s} (0)   \\hat{c}^{\\dagger}\_{y,s} (0)   \\rangle \\rangle    \\nonumber \\\\
    \\texttt{GTT(x,y,s) }  &=&   \\phantom{+} \\langle \\langle \\hat{c}^{\\phantom\\dagger}\_{x,s} (Nt \\Delta \\tau)   \\hat{c}^{\\dagger}\_{y,s} (Nt \\Delta \\tau)   \\rangle \\rangle    \\end{aligned}$$
 In the above we have omitted the color index since the Green functions are color independent. The time displaced spin-spin correlations: $ 4 \\langle \\langle \\hat{S}^{z}\_{\\vec{i}} (\\tau)  \\hat{S}^{z}\_{\\vec{j}} (0)\\rangle \\rangle   $ are thereby given by:
$$4 \\langle \\langle \\hat{S}^{z}\_{\\vec{i}} (\\tau)  \\hat{S}^{z}\_{\\vec{j}} (0)\\rangle \\rangle   = - 2 \\; \\texttt{G0T(J,I,1) } \\texttt{GT0(I,J1) }$$
 Note that the above holds for the SU(2) HS transformation discussed in this chapter. The handling of time-displaced correlation functions is identical to that of equal time correlations.

Walkthrough: the *M*<sub>*z*</sub>-Hubbard model on a square lattice
====================================================================

The Hubbard Hamiltonian can equally be written as:
$$\\label{eqn\_hubbard\_Mz}
\\mathcal{H}=
\\sum\\limits\_{\\sigma=1}^{2} 
\\sum\\limits\_{x,y =1 }^{N\_{unit\\; cells }} 
  c^{\\dagger}\_{x \\sigma} T\_{x,y}c^{\\phantom\\dagger}\_{y \\sigma} 
- \\frac{U}{2}\\sum\\limits\_{x}\\left\[
c^{\\dagger}\_{x, \\uparrow} c^{\\phantom\\dagger}\_{x \\uparrow}  -   c^{\\dagger}\_{x, \\downarrow} c^{\\phantom\\dagger}\_{x \\downarrow}  \\right\]^{2}\\;.$$
 We can make contact with the general form of the Hamiltonian (see Eq. \[eqn:general\_ham\]) by setting: *N*<sub>*f**l*</sub> = 2, *N*<sub>*c**o**l*</sub> ≡ `N_SUN` = 1, *M*<sub>*T*</sub> = 1, *T*<sub>*x**y*</sub><sup>(*k**s*)</sup> = *T*<sub>*x*, *y*</sub>, *M*<sub>*V*</sub> = *N*<sub>unit cell</sub>, $U\_{k}       =   \\frac{U}{2}$, *V*<sub>*x**y*</sub><sup>(*k*, *s* = 1)</sup> = *δ*<sub>*x*, *y*</sub>*δ*<sub>*x*, *k*</sub>, *V*<sub>*x**y*</sub><sup>(*k*, *s* = 2)</sup> = −*δ*<sub>*x*, *y*</sub>*δ*<sub>*x*, *k*</sub>, *α*<sub>*k**s*</sub> = 0 and *M*<sub>*I*</sub> = 0. The coupling of the HS to the z-component of the magnetization breaks the SU(2) spin symmetry. Nevertheless the z-component of spin remains a good quantum number such that the imaginary time propagator – for a given HS field – is block diagonal in this quantum number. This corresponds to the flavor index which runs from one to two labelling spin up and spin down degrees of freedom. In the parameter file listed in Sec. \[sec:input\] setting the model variable to `Hubbard_Mz` will carry the simulation in the above representation. With respect to the SU(2) case, the changes required in the `Hamiltonian_Examples.f90` module are minimal and essentially effect only the interaction term, and calculation of observables. We note that in this formulation the hopping matrix can be favor dependent such that a Zeeman magnetic field can be introduced. If the chemical potential is set to zero, this will not generate a negative sign problem . A sample run for this model can be found in `Examples/Hubbard_Mz_Square/`.

The interaction term: `Call Ham_V` 
-----------------------------------

The interaction term is now given by:


    Allocate(Op_V(Ndim,N_FL))
    do nf = 1,N_FL
       do i  = 1, Ndim
          Call Op_make(Op_V(i,nf),1)
       enddo
    enddo
    Do nf = 1,N_FL
       nc = 0
       X = 1.d0
       if (nf == 2) X = -1.d0
       Do i = 1,Ndim
          nc = nc + 1
          Op_V(nc,nf)%P(1) = I
          Op_V(nc,nf)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.D0))
          Op_V(nc,nf)%g      = X*SQRT(CMPLX(DTAU*ham_U/2.d0, 0.D0, kind(0.D0))) 
          Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
          Op_V(nc,nf)%type   = 2
          Call Op_set( Op_V(nc,nf) )
       Enddo
    Enddo

In the above, one will see explicitly that there is a sign difference between the coupling of the HS field in the two flavor sectors.

The measurements: `Call Obser, Call ObserT` 
--------------------------------------------

Since the spin up and spin down Green functions differ for a given HS configuration, the Wick decomposition will take a different form. In particular, the equal time spin-spin correlation functions, $ 4 \\langle \\langle \\hat{S}^{z}\_{\\vec{i}}   \\hat{S}^{z}\_{\\vec{j}} \\rangle \\rangle   $, calculated in the subroutine `Obser`, will take the form:
$$\\begin{aligned}
   4 \\langle \\langle \\hat{S}^{z}\_{x}   \\hat{S}^{z}\_{y} \\rangle \\rangle   =  & &  \\texttt{  GRC(x,y,1) \* GR(x,y,1) + GRC(x,y,2) \* GR(x,y,2) + }  \\nonumber \\\\ 
& & \\texttt{   (GRC(x,x,2) - GRC(x,x,1))\*(GRC(y,y,2) - GRC(y,y,1))}  \\nonumber
  \\end{aligned}$$
 Here, `GRC` is defined in Eq. \[GRC.eq\]. Equivalent changes will have to be carried out for other equal time and time displaced observables.

Apart from these modifications, the program will run in exactly the same manner as for the SU(2) case.

Walkthrough: the *S**U*(2)-Hubbard model on the honeycomb lattice
=================================================================

The Hamilton module `Hamiltonian_Examples.f90` can also carry out simulations for the the Hubbard model on the Honeycomb lattice by setting in the parameter file (see Sec. \[sec:input\]) `Lattice_type = Honeycomb`. A sample run for this model can be found in `Examples/Hubbard_SU2_Honeycomb/`.

Working with multi-orbital unit cells: `Call Ham_Latt` 
-------------------------------------------------------

This model is an example of a multi-orbital unit cell, and the aim of this section is to document how implement this in the code. The Honeycomb lattice is a triangular Bravais lattice the two orbitals per unit cell. The routine `Ham_Latt` will set:


    Norb    = 2
    N_coord = 3
    a1_p(1) =  1.D0   ; a1_p(2) =  0.d0
    a2_p(1) =  0.5D0  ; a2_p(2) =  sqrt(3.D0)/2.D0             
    L1_p    =  dble(L1) * a1_p
    L2_p    =  dble(L2) * a2_p
                

and then call ` Make_Lattice( L1_p, L2_p, a1_p, a2_p, Latt ) ` so as to generate the triangular lattice. The coordination number of this lattice is ` N_coord=3 ` and the number of orbitals per unit cell corresponds to `NORB=2`. The total number of orbitals is thereby: `N_{\mathrm{dim}}=Latt%N*NORB`. To easily keep track of the orbital and unit cell, we define a super-index as shown below:


    Allocate (List(Ndim,2), Invlist(Latt%N,Norb))
    nc = 0
    Do I = 1,Latt%N           ! Unit-cell index 
       Do no = 1,Norb         ! Orbital index
          nc = nc + 1         ! Super index labeling unit cell and orbital
          List(nc,1) = I      ! Unit-cell  of super index  nc
          List(nc,2) = no     ! Orbital of super inde nc
          Invlist(I,no) = nc  ! Super index for given  unit cell and orbital
      Enddo
    Enddo

With the above lists one can run trough all the orbitals and at each time keep track of the unit-cell and orbital index. We note that when translation symmetry is completely absent one can work with on unit cell, and the number of orbitals will then correspond to the number of lattice sites.

The hopping term: `Call Ham_Hop` 
---------------------------------

Some care has to be taken when setting the hopping matrix. In the Hamilton module we do this in the following way.


    DO I = 1, Latt%N                              ! Loop over unit cell 
       do no = 1,Norb                             ! Runs over orbitals and sets chemical potential
          I1 = invlist(I,no)  
          Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
       enddo
       I1 = Invlist(I,1)                          ! Orbital A of unit cell I
       Do nc1 = 1,N_coord                         ! Loop over coordination  number
          select case (nc1)
          case (1)
             J1 = invlist(I,2)                    ! Orbital B of unit cell i
          case (2)
             J1 = invlist(Latt%nnlist(I,1,-1),2)  ! Orbital B of unit cell i + a_1 - a_2
          case (3)
             J1 = invlist(Latt%nnlist(I,0,-1),2)  ! Orbital B of unit cell i - a_1 
          case default
             Write(6,*) ' Error in  Ham_Hop '  
          end select
          Op_T(nc,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0))   
          Op_T(nc,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
       Enddo
    Enddo

As apparent from the above, hopping matrix elements are non-zero only between the A and B sublattices.

Observables: `Call Obser`, `Call ObserT`
----------------------------------------

In the multi-orbital case, the correlation functions have additional orbital indices. This is automatically taken care of in the routines `Call Obser` and `Call ObserT` since we considered the Hubbard model on the square lattice to correspond to a multi-orbital unit cell albeit with the special choice of one orbital per unit cell.

Walkthrough: the *S**U*(2)-Hubbard model on a square lattice coupled to a transverse Ising field
================================================================================================

The model we consider here is very similar to the above, but has an additional coupling to a transverse field.
$$\\begin{aligned}
\\label{eqn\_hubbard\_sun\_Ising}
\\mathcal{H}= & & 
\\sum\\limits\_{\\sigma=1}^{2} 
\\sum\\limits\_{x,y } 
  c^{\\dagger}\_{x \\sigma} T\_{x,y}c^{\\phantom\\dagger}\_{y \\sigma} 
+ \\frac{U}{2}\\sum\\limits\_{x}\\left\[
\\sum\\limits\_{\\sigma=1}^{2}
\\left(  c^{\\dagger}\_{x \\sigma} c^{\\phantom\\dagger}\_{x \\sigma}  -1/2 \\right) \\right\]^{2}   
+  \\xi \\sum\_{\\sigma,\\langle x,y \\rangle} \\hat{Z}\_{\\langle x,y \\rangle}  \\left( c^{\\dagger}\_{x \\sigma} c^{\\phantom\\dagger}\_{y \\sigma}  + h.c. \\right)  \\nonumber \\\\ 
 & & - h \\sum\_{\\langle x,y \\rangle} \\hat{X}\_{\\langle x,y \\rangle}   - J \\sum\_{\\langle \\langle x,y \\rangle \\langle x',y' \\rangle \\rangle} 
  \\hat{Z}\_{\\langle x,y \\rangle}   \\hat{Z}\_{\\langle x',y' \\rangle} \\end{aligned}$$
 We can make contact with the general form of the Hamiltonian by setting: *N*<sub>*f**l*</sub> = 1, *N*<sub>*c**o**l*</sub> ≡ `N_SUN` = 2, *M*<sub>*T*</sub> = 1, *T*<sub>*x**y*</sub><sup>(*k**s*)</sup> = *T*<sub>*x*, *y*</sub>, *M*<sub>*V*</sub> = *N*<sub>unit cell</sub> ≡ *N*<sub>*d**i**m*</sub>, $U\_{k}       =   -\\frac{U}{2}$, *V*<sub>*x**y*</sub><sup>(*k**s*)</sup> = *δ*<sub>*x*, *y*</sub>*δ*<sub>*x*, *k*</sub>, $\\alpha\_{ks}   = - \\frac{1}{2}  $ and *M*<sub>*I*</sub> = 2*N*<sub>unit cell</sub>. The last two terms of the above Hamiltonian describes a transverse Ising field model on the bonds of the square lattice. This type of Hamiltonian has recently been extensively discussed . Here we adopt the notation of Ref. . Note that ⟨⟨*x*, *y*⟩⟨*x*′,*y*′⟩⟩ denotes nearest neighbor bonds. The modifications required to generalize the Hubbard model code to the above model are two-fold.

Firstly, one has to specify the function `Real (Kind=8) function S0(n,nt)` and secondly modify the interaction `Call Ham_V`.

A sample run for this model can be found in `Examples/Hubbard_SU2_Ising_Square/`.

The interaction term: `Call Ham_V`
----------------------------------

The dimension of `Op_V` is now (*M*<sub>*I*</sub> + *M*<sub>*V*</sub>)×*N*<sub>*f**l*</sub> = (3 \* *N*<sub>*d**i**m*</sub>)×1. We set the effective dimension for the Hubbard term to *N* = 1 and to *N* = 2 for the Ising term. The allocation of this array of operators reads:


    do i  = 1,N_coord*Ndim    !  Runs  over bonds for Ising variable
      call Op_make(Op_V(i,1),2)
    enddo
    do i  =  N_coord*Ndim+1, (N_coord+1)*Ndim   ! Runs over sites for Hubbard
      call Op_make(Op_V(i,1),1)
    enddo

The first `N_coord*Ndim` operators run through the 2N bonds of the square lattice and are given by:


    Do nc = 1,Ndim*N_coord   ! Runs over bonds.  Coordination number = 2.
                             ! For the square lattice Ndim = Latt%N
      
       I1 = L_bond_inv(nc,1) ! Site one of the bond.  
    	                 ! L_bond_inv is setup in Setup_Ising_action
       if ( L_bond_inv(nc,2)  == 1 ) I2 = Latt%nnlist(I1,1,0)   ! Site two of the bond 
       if ( L_bond_inv(nc,2)  == 2 ) I2 = Latt%nnlist(I1,0,1) 
       Op_V(nc,1)%P(1) = I1
       Op_V(nc,1)%P(2) = I2
       Op_V(nc,1)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0))
       Op_V(nc,1)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
       Op_V(nc,1)%g = cmplx(-dtau*Ham_xi,0.D0,kind(0.D0))
       Op_V(nc,1)%alpha = cmplx(0d0,0.d0, kind(0.D0))
       Op_V(nc,1)%type =1
    Enddo

Here, `ham_xi` defines the coupling strength between the Ising and fermion degree of freedom. As for the Hubbard case, the first `Ndim` operators read:


    nc = N_coord*Ndim 
    Do nc  = i = 1, Ndim
        nc = nc + 1
        Op_V(nc,1)%P(1)   = i 
        Op_V(nc,1)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
        Op_V(nc,1)%g      = sqrt(cmplx(-dtau*ham_U/(DBLE(N_SUN)), 0.D0, kind(0.D0)))
        Op_V(nc,1)%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
        Op_V(nc,1)%type   = 2
    Enddo

The function `Real (Kind=8) function S0(n,nt)` 
-----------------------------------------------

As mentioned above, a configuration is given by
*C* = {*s*<sub>*i*, *τ*</sub>,*l*<sub>*j*, *τ*</sub> with *i*=1⋯*M*<sub>*I*</sub>,*j*=1⋯*M*<sub>*V*</sub>,*τ*=1,*L*<sub>*T**r**o**t**t**e**r*</sub>}
 and is stored in the integer array `nsigma(M_V + M_I, Ltrot)`. With the above ordering of Hubbard and Ising interaction terms, and a for a given imaginary time, the first fields corresponds to the Hubbard interaction and the next `2*Ndim` ones to the Ising interaction. The first argument of the function `S0`, `n`, corresponds to the index of the operator string `Op_V(n,1)`. If `Op_V(n,1)%type = 2`, `S0(n,nt)` returns 1. If `Op_V(n,1)%type = 1` then function `S0` returns
$$\\frac{e^{-S\_{0,I} \\left(  s\_{1,\\tau},  \\cdots,  - s\_{m,\\tau},  \\cdots s\_{M\_I,\\tau}   \\right) } }{e^{-S\_{0,I}  \\left(  s\_{1,\\tau},  \\cdots,   s\_{m,\\tau},  \\cdots s\_{M\_I,\\tau}   \\right)   } }$$
 That is, `S0(n,nt)` returns the ratio of the new to old weight of the Ising Hamiltonian upon flipping a single Ising spin *s*<sub>*m*, *τ*</sub>. Note that in this specific case ` m = n - Ndim `

Other models
============

The aim of this section is to briefly mention a small selection of other models that can be simulated within the ALF-project.

The Kondo lattice
-----------------

Simulating the Kondo lattice within the ALF-project requires rewriting of the model along the lines of Refs. . Adopting the notation of these articles, the Hamiltonian that one will simulate reads:
$$\\hat{\\mathcal{H}}  = 
	\\underbrace{-t \\sum\_{\\langle  \\vec{i},\\vec{j} \\rangle,\\sigma} \\left( \\hat{c}\_{\\vec{i},\\sigma}^{\\dagger}  \\hat{c}\_{\\vec{j},\\sigma}^{\\phantom\\dagger}   + \\text{H.c.} \\right) }\_{\\equiv \\hat{\\mathcal{H}}\_t} - \\frac{J}{4} 
	\\sum\_{\\vec{i}} \\left( \\sum\_{\\sigma} \\hat{c}\_{\\vec{i},\\sigma}^{\\dagger}  \\hat{f}\_{\\vec{i},\\sigma}^{\\phantom\\dagger}  + 
	                                                        \\hat{f}\_{\\vec{i},\\sigma}^{\\dagger}  \\hat{c}\_{\\vec{i},\\sigma}^{\\phantom\\dagger}   \\right)^{2}   +
        \\underbrace{\\frac{U}{2}   \\sum\_{\\vec{i}}   \\left( \\hat{n}^{f}\_{\\vec{i}} -1 \\right)^2}\_{\\equiv \\hat{\\mathcal{H}}\_U}.$$
 This from is included in Eq. \[eqn:general\_ham\_i\] such the above Hamiltonian can be implemented in our program package. The relation to the Kondo lattice model follows from expanding the square of the hybridization to obtain:
$$\\hat{\\mathcal{H}}  =\\hat{\\mathcal{H}}\_t   
	+ J \\sum\_{\\vec{i}}  \\left(  \\hat{\\vec{S}}^{c}\_{\\vec{i}} \\cdot  \\hat{\\vec{S}}^{f}\_{\\vec{i}}    +   \\hat{\\eta}^{z,c}\_{\\vec{i}} \\cdot  \\hat{\\eta}^{z,f}\_{\\vec{i}}  
		-  \\hat{\\eta}^{x,c}\_{\\vec{i}} \\cdot  \\hat{\\eta}^{x,f}\_{\\vec{i}}  -  \\hat{\\eta}^{y,c}\_{\\vec{i}} \\cdot  \\hat{\\eta}^{y,f}\_{\\vec{i}} \\right) 
	 + \\hat{\\mathcal{H}}\_U.$$
 where the *η*-operators relate to the spin-operators via a particle-hole transformation in one spin sector:
$$\\hat{\\eta}^{\\alpha}\_{\\vec{i}}  = \\hat{P}^{-1}  \\hat{S}^{\\alpha}\_{\\vec{i}} \\hat{P}  	\\; \\text{ with }  \\;   
	\\hat{P}^{-1}  \\hat{c}^{\\phantom\\dagger}\_{\\vec{i},\\uparrow} \\hat{P}  =   (-1)^{i\_x+i\_y} \\hat{c}^{\\dagger}\_{\\vec{i},\\uparrow}  \\; \\text{ and }  \\;   
	\\hat{P}^{-1}  \\hat{c}^{\\phantom\\dagger}\_{\\vec{i},\\downarrow} \\hat{P}  = \\hat{c}^{\\phantom\\dagger}\_{\\vec{i},\\downarrow}$$
 Since the $\\hat{\\eta}^{f} $ and $ \\hat{S}^{f} $ operators do not alter the parity \[$(-1)^{\\hat{n}^{f}\_{\\vec{i}}}$ \] of the *f*-sites,
$$\\left\[  \\hat{\\mathcal{H}}, \\hat{\\mathcal{H}}\_U \\right\] = 0.$$
 Thereby, and for positive values of *U* , doubly occupied or empty *f*-sites corresponding to even parity will be suppressed by a Boltzmann factor *e*<sup>−*β**U*/2</sup> is comparison to odd parity ones. Choosing *β**U* adequately will essentially allow to restrict the Hilbert space to odd parity *f*-sites. In this Hilbert space $\\hat{\\eta}^{x,f} = \\hat{\\eta}^{y,f} =  \\hat{\\eta}^{z,f} =0$ such that the Hamiltonian reduces to the Kondo lattice model.

SU(N) Hubbard-Heisenberg models
-------------------------------

SU(2N) Hubbard-Heisenberg models can be written as:
$$\\hat{\\mathcal{H}}  =  
 \\underbrace{ - t \\sum\_{ \\langle \\vec{i},\\vec{j} \\rangle }    \\left(  \\vec{\\hat{c}}^{\\dagger}\_{\\vec{i}}  \\vec{\\hat{c}}^{\\phantom{\\dagger}}\_{\\vec{j}} + \\text{H.c.} \\right) }\_{\\equiv \\hat{\\mathcal{H}}\_t} \\; \\; 
\\underbrace{ -\\frac{J}{2 N}  \\sum\_{ \\langle \\vec{i},\\vec{j} \\rangle  } \\left(
           \\hat{D}^{\\dagger}\_{ \\vec{i},\\vec{j} }\\hat{D}^{\\phantom\\dagger}\_{ \\vec{i},\\vec{j}}  +
            \\hat{D}^{\\phantom\\dagger}\_{ \\vec{i},\\vec{j} } \\hat{D}^{\\dagger}\_{ \\vec{i},\\vec{j} }  \\right) }\_{\\equiv\\hat{\\mathcal{H}}\_J}
            + 
 \\underbrace{\\frac{U}{N}  \\sum\_{\\vec{i}} \\left(
             \\vec{\\hat{c}}^{\\dagger}\_{\\vec{i}}  \\vec{\\hat{c}}^{\\phantom\\dagger}\_{\\vec{i}} -  {\\frac{N}{2} } \\right)^2}\_{\\equiv \\hat{\\mathcal{H}}\_U}$$
 Here, $ \\vec{\\hat{c}}^{\\dagger}\_{\\vec{i}} =
(\\hat{c}^{\\dagger}\_{\\vec{i},1},  \\hat{c}^{\\dagger}\_{\\vec{i},2}, \\cdots, \\hat{c}^{\\dagger}\_{\\vec{i}, N } ) $ is an *N*-flavored spinor, and $ \\hat{D}\_{ \\vec{i},\\vec{j}} = \\vec{\\hat{c}}^{\\dagger}\_{\\vec{i}}
\\vec{\\hat{c}}\_{\\vec{j}}  $. To use the present package to simulate this model, one will rewrite the *J*-term as a sum of perfect squares,
$$\\hat{\\mathcal{H}}\_J =  -\\frac{J}{4 N}  \\sum\_{  \\langle \\vec{i}, \\vec{j} \\rangle }
        \\left(\\hat{D}^{\\dagger}\_{  \\langle \\vec{i}, \\vec{j} \\rangle  } +  \\hat{D}\_{  \\langle \\vec{i}, \\vec{j} \\rangle }  \\right)^2  -
        \\left(\\hat{D}^{\\dagger}\_{   \\langle \\vec{i}, \\vec{j} \\rangle } -  \\hat{D}\_{  \\langle \\vec{i}, \\vec{j} \\rangle}  \\right)^2,$$
 so to manifestly bring it into the form of Eq. \[eqn:general\_ham\_i\]. It is amusing to note that setting the hopping *t* = 0, charge fluctuations will be suppressed by the Boltzmann factor $e^{\\beta U /N \\left(  \\vec{\\hat{c}}^{\\dagger}\_{\\vec{i}}  \\vec{\\hat{c}}^{\\phantom\\dagger}\_{\\vec{i}} -  {\\frac{N}{2} } \\right)^2 } $ since in this case $ \\left\[   \\hat{\\mathcal{H}}, \\hat{\\mathcal{H}}\_U \\right\]  = 0 $. This provides a route to use the auxiliary field QMC algorithm to simulate – free of the sign problem – SU(2N) Heisenberg models in the self-adjoint antisymmetric representation [3] For odd values of *N* recent progress in our understanding of the origins of the sign problem will allows us to simulate – without encountering the sign problem – a set of non-trivial Hamiltonians .

Hubbard model in the canonical ensemble
---------------------------------------

To simulate the Hubbard model in the canonical ensemble one can add the constraint:
$$\\hat{\\mathcal{H}}   = \\hat{\\mathcal{H}\_{tU}}     + \\underbrace{\\lambda \\left( \\hat{N} -  N \\right)^{2}}\_{\\equiv \\hat{H}\_\\lambda }$$
 In the limit *λ* → ∞ the uniform charge fluctuations,
$$S ( \\vec{q} = 0)   =  \\sum\_{\\vec{r}}   \\left\[ \\langle \\hat{n}\_{\\vec{r}}  \\hat{n}\_{\\vec{0}} \\rangle  - \\langle \\hat{n}\_{\\vec{r}}\\rangle \\langle  \\hat{n}\_{\\vec{0}} \\rangle  \\right\]$$
 are suppressed and the grand-canonical simulation maps onto a canonical one. To implement this in the code we have adopted the following strategy. Since $ \\left( \\hat{N} -  N \\right)^{2}  $ effectively corresponds to a long range interaction one may face the issue that the acceptance rate of a single HS flip becomes very small on large lattices. To circumvent this problem we have used the following decomposition:
$$e^{-\\beta \\hat{H}}  =   \\prod\_{\\tau = 1}^{L\_{\\text{Trotter}}} \\left\[  e^{-\\Delta \\tau \\hat{H}\_t} e^{-\\Delta \\tau \\hat{H}\_U}  
	\\underbrace{e^{-\\frac{\\Delta \\tau}{n\_{\\lambda}} \\hat{H}\_{\\lambda} } \\cdots e^{-\\frac{\\Delta \\tau}{n\_{\\lambda}} \\hat{H}\_{\\lambda} } }\_{n\_\\lambda \\text{-times } }\\right\].$$
 Thereby per time slice we need *n*<sub>*λ*</sub> fields to impose the constraint. Since for each field the coupling constant is suppressed by a factor *n*<sub>*λ*</sub>, we can monitor the acceptance. An implementation of this program can be found in `Prog/Hamiltonian_Hub_Canonical.f90` and a test run in the directory `Examples/Hubbard_Mz_Square_Can`

 Analysis programs 
===================

| Program        | Description                                                                                      |
|:---------------|:-------------------------------------------------------------------------------------------------|
| `cov_scal.f90` | In combination with the script `analysis.sh`, the bin files with suffix `_scal` are read in,     |
|                | and corresponding file with suffix `_scalJ` are produced. They contain the result of the         |
|                | Jackknife resampling.                                                                            |
| `cov_eq.f90`   | In combination with the script `analysis.sh`, the bin files with suffix `_eq` are read in,       |
|                | and corresponding files will suffix `_eqJR` and `_eqJK` are produced. They correspond to         |
|                | correlation functions in real and Fourier space, respectively.                                   |
| `cov_tau.f90`  | In combination with the script `analysis.sh`, the bin files `X_tau` are read in,                 |
|                | and the directories `X_kx_ky` are produced for all `kx` and `ky` greater or equal to zero.       |
|                | Here `X` is a place holder from `Green`, `SpinXY`, etc as specified in ` Alloc_obs(Ltau)`        |
|                | (See section \[Alloc\_obs\_sec\]). Each directory contains a file `g_kx_ky` containing the       |
|                | time displaced correlation function traced over the orbitals. It also contains the               |
|                | covariance matrix if `N_cov` is set to unity in the parameter file listed in Sec. \[sec:input\]. |
|                | Equally, a directory `X_R0` for the local time displaced correlation function is generated.      |

Here we briefly discuss the analysis programs which read in bins and carry out the error analysis. Error analysis is based on the central limit theorem, which required bins to be statistically independent, and also the existence of a well-defined variance of the distribution. The former will be the case if bins are longer than the auto-correlation time. The latter has to be checked by the user, since in general the distribution variance depends on the model and on the observable. In the parameter file listed in Sec. \[sec:input\], the user can specify the how many initial bins should be omitted (variable `n_skip`). This number should be comparable to the auto-correlation time. The re-binning variable `N_rebin` will merge `N_rebin` bins into a single one. If the autocorrelation time is smaller than the effective bin size, then the error should be independent on the bin size and thereby on the variable `N_rebin`. Our analysis is based on the Jackknife resampling. As listed in Table, \[table:analysis\_programs\] we provide three programs to account for the three observable types. The programs can be found in the directory `Analysis` and are executed by running the bash shell script `analysis.sh`

| File                      | Description                                              |
|:--------------------------|:---------------------------------------------------------|
| `parameters`              | Contains also variables for the error analysis:          |
|                           | `n_skip`, `N_rebin` and `N_Cov` (see Sec. \[sec:input\]) |
| `X_scal`, `Y_eq`, `Y_tau` | Monte Carlo bins (see Table \[table:output\])            |

| File                  | Description                                                                              |     |
|:----------------------|:-----------------------------------------------------------------------------------------|:----|
| `X_scalJ`             | Jackknife mean and error of `X`, where `X` stands for `Kin, Pot, Part`, and `Ener`.      |     |
| `Y_eqJR` and `Y_eqJK` | Jackknife mean and error of `Y`, where `Y` stands for `Green, SpinZ, SpinXY`, and `Den`. |     |
|                       | The suffixes `R` and `K` refers to real and reciprocal space, respectively.              |     |
| `Y_R0/g_R0`           | Time-resolved and spatially local Jackknife mean and error of `Y`,                       |     |
|                       | where `Y` stands for `Green, SpinZ, SpinXY`, and `Den`.                                  |     |
| `Y_kx_ky/g_kx_ky`     | Time-resolved and $\\vec{k}$-dependent Jackknife mean and error of `Y`,                  |     |
|                       | where `Y` stands for `Green, SpinZ, SpinXY`, and `Den`.                                  |     |

In the following, we describe the formatting of the output files mentioned in Table \[table:analysis\_output\].

-   For the scalar quantities `X`, the output files `X_scalJ` have the following formatting:

    `Effective number of bins, and bins:           <N_bin - n_skip>          <N_bin>`
    `OBS :    1      <mean(X)>      <error(X)>`
    `OBS :    2      <mean(sign)>   <error(sign)>`

-   For the equal-time correlation functions `Y`, the formatting of the output files `Y_eqJR` and `Y_eqJK` follows this structure:

    `do i = 1, N_unit_cell`
    `   <k_x(i)>   <k_y(i)>`
    `   do alpha = 1, N_orbital`
    `   do beta  = 1, N_orbital`
    `      alpha   beta   Re<mean(Y)>   Re<error(Y)>   Im<mean(Y)>   Im<error(Y)>`
    `   enddo`
    `   enddo`
    `enddo`

    where `Re` and `Im` refer to the real and imaginary part, respectively.

-   The time-displaced correlation functions `Y` are written to the output files `Y_R0/g_R0`, when measured locally in space, and to the output files `Y_kx_ky/g_kx_ky` when they are measured $\\vec{k}$-resolved. Both output files have the following formatting:

    `do i = 0, Ltau`
    `   tau(i)   <mean( Tr[Y] )>   <error( Tr[Y])>`
    `enddo`

    where `Tr` corresponds to the trace over the orbital degrees of freedom.

Running the code
================

In this section we describe the steps to compile and run the code and to perform the error analysis of the data.

Compilation
-----------

The environment variables are defined in the bash script `set_env.sh` as follows:


    # Description of PROGRAMMCONFIGURATION:
    # -DMPI selects MPI.
    # Setting nothing  compiles without mpi.
    # -DQRREF selects  a reference implementation of the QR decomposition. 
    # Setting nothing selects system lapack for the QR decomposition.
    # -DSTAB1 selects an alternative stabilization scheme.
    # Setting nothing selects the default stabilizatiion
    PROGRAMMCONFIGURATION=""
    f90="gfortran"
    export f90
    F90OPTFLAGS="-O3"
    export F90OPTFLAGS
    FL="-c ${F90OPTFLAGS} ${PROGRAMMCONFIGURATION}"
    export FL
    DIR=`pwd`
    export DIR
    Libs=${DIR}"/Libraries/"
    export Libs
    LIB_BLAS_LAPACK="-llapack -lblas"
    export LIB_BLAS_LAPACK

In the above, the GNU Fortan compiler `gfortran` is set.[4] The program can be compiled and ran either in single-thread mode (default) or in multi-threading mode (define `-DMPI`) using the MPI standard for parallelization. To compile the libraries, the analysis programs and the quantum Monte Carlo program, the following steps should be executed:

1.  Export the environment variables:

        source set_env.sh

2.  Compile the libraries and the error analysis routines

        cd Libraries
        make
        cd ..
        cd Analysis
        make
        cd ..

3.  Compile the quantum Monte Carlo code

        cd Prog
        make
        cd ..

Starting a simulation
---------------------

To start a simulation from scratch, the following files have to be present: `parameters` and `seeds`. To run a single-thread simulation for one of the Hubbard model described in Sec. \[sec:walk1\] - \[sec:walk2\], issue the command

    ./Prog/Examples.out

To restart the code using an existing simulation as a starting point, first run the script `out_to_in.sh` to set the input configuration files.

Error analysis
--------------

To perform an error analysis, based on the jackknife scheme, of the Monte Carlo bins for all observables run the script `analysis.sh` (see Sec. \[sec:analysis\]).

Performance
===========

Next to the entire computational time is spent in BLAS routines such that the performance of the code will depend on the implementation of this library. We have found that the code performs well, and that an efficient OpenMP version can be obtained merely by loading the corresponding BLAS and LAPACK routines.

Conclusions and future directions
=================================

In its present form, the ALF-project allows to simulate a very large class of non-trivial models efficiently and at a minimal programming cost. There are many possible extensions which deserve to be considered in future releases. The Hamiltonians we presently defining are imaginary time independent. This however, can be easily generalized to time dependent Hamiltonians thus allowing, for example, to access entanglement properties of interacting fermionic systems . Generalizations to include global moves are equally desirable. This is a prerequisite to play with recent ideas of self-learning algorithms so as to possibly avoid critical slowing down. At present we are restricted to discrete fields such that implementations of the long range Coulomb repulsion as introduced in is not included in the package. Extensions to continuous fields are certainly possible, but require an efficient upgrading scheme. Finally, a ground state projective formulation is equally desirable.

Acknowledgments
===============

We are very grateful to S. Beyl, M. Hohenadler, F. Parisen Toldin, M. Raczkowski, J. Schwab, T. Sato, Z. Wang and M. Weber, for constant support during the development of this project. FFA would also like to thank T. Lang and Z.Y. Meng for developments of the auxiliary field code as well as T. Grover. MB thanks the the Bavarian Competence Network for Technical and Scientific High Performance Computing (KONWIHR) for financial support. FG and JH thank the SFB-1170 for financial support under projects Z03 and C01. FFA thanks the DFG-funded FOR1807 and FOR1346 for partial financial support. Calculations to extensively test this package were carried out on SuperMUC at the Leibniz Supercomputing Centre and on JURECA at the Jülich Supercomputing Centre (JSC). We thank those institutions for generous computer allocations.

License
=======

The ALF code is provided as an open source software such that it is available to all and we hope that it will be useful. If you benefit from this code we ask that you acknowledge the ALF collaboration as mentioned on our homepage [alf.physik.uni-wuerzburg.de](alf.physik.uni-wuerzburg.de). The git repository at [alf.physik.uni-wuerzburg.de](alf.physik.uni-wuerzburg.de) gives us the tools to create a small but vibrant community around the code and provides a suitable entrypoint for future contributors and future developments. The homepage is also the place where the original source files can be found. With the coming public release it was necessary to add copyright headers to our source files. The Creative Commons licenses are a good way to share our documentation and it is also well accepted by publishers. Therefore this documentation is licensed to you under a CC-BY-SA license. This means you can share it and redistribute it as long as you cite the original source and license your changes under the same license. The details are in the file license.CCBYSA that you should have received with this documentation. The source code itself is licensed under a GPL license to keep the source as well as any future work in the community. To express our desire for a proper attribution we decided to make this a visible part of the license. To that end we have exercised the rights of section 7 of GPL version 3 and have amended the license terms with an additional paragraph that expresses our wish that if an author has benfitted from this code that he/she should consider giving back a citation as specified on [alf.physik.uni-wuerzburg.de](alf.physik.uni-wuerzburg.de). This is not something that is meant to restrict your freedom of use, but something that we strongly expect to be good scientific conduct. The original GPL license can be found in the file license.GPL and the additional terms can be found in license.additional. In favour to our users, *ALF* contains part of the lapack implementation version 3.6.1 from <http://www.netlib.org/lapack>. Lapack is licensed under the modified BSD license whose full text can be found in license.BSD.
With that being said, we hope that ALF will prove to you to be a suitable and highly performant tool that enables you to perform Monte Carlo studies of solid state models of unprecedented complexity.
The ALF project’s contributors.

COPYRIGHT
---------

Copyright  2016, The *ALF* Project.
The ALF Project Documentation is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License. You are free to share and benefit from this documentation as long as this license is preserved and proper attribution to the authors is given. For details see the ALF project homepage [alf.physik.uni-wuerzburg.de](alf.physik.uni-wuerzburg.de) and the file `license.CCBYSA`.

[1] The attentive reader will have noticed that for arbitrary Trotter decompositions, the imaginary time propagator is not necessarily Hermitian. Thereby, the above equation is correct only up to corrections stemming from the controlled Trotter systematic error.

[2] We have encountered some compiling issues with this flag. In particular the older intel ifort compiler version 10.1 fails for all optimization levels.

[3] This corresponds to a Young tableau with single column and *N*/2 rows.

[4] A known issue with the alternative Intel Fortran compiler `ifort` is the handling of automatic, temporary arrays which `ifort` allocates on the stack. For large system sizes and/or low temperatures this may lead to a runtime error. One solution is to demand allocation of arrays above a certain size on the heap instead of the stack. This is accomplished by the `ifort` compiler flag `-heap-arrays [n]` where `[n]` is the minimal size (in kilobytes, for example `n=1024`) of arrays that are allocated on the heap.
