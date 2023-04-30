# FluidSim_PT
This project is a based on the works of Mathias Muelle from TenMinutePhysics ( https://www.youtube.com/@TenMinutePhysics) and is an extension of their already written code. The idea behind this project is to improve the already existing simulator by adding functions that simulate the behavior of a fluid more accurately.
The aim was not to build a simulator that performs on the same level as commercial products, but to build an easily editable and expandable framework. This project is meant as a tool to teach the application of numerical approaches to difficult topics, such as the Navier-Stokes equation.
Therefore it is written in plain JavaScript up to this point and requires no additional dependencies.
For a more extensive library of examples, e.g. a narrowing pipe, I will probably switch the project over to p5.js.

The project is still not entirely finished but incorporates methods to allow for the modelling the effects of a pressure gradient, viscosity and diffusion on the flow of the fluid.

It is still assumed, that the fluid is incompressible, viscosity and density are constant and that the process is isothermal.

This project is a based on the works of [10min physics] and is an extension of their already written code. The idea behind this project is to improve the already existing simulator by adding functions that simulate the behavior of a fluid more accurately.
The aim was not to build a simulator that performs on the same level as commercial products, but to build an easily editable and expandable framework. This project is meant as a tool to teach the application of numerical approaches to difficult topics, such as the Navier-Stokes equation.
Therefore it is written in plain JavaScript up to this point and requires no additional dependencies.
For a more extensive library of examples, e.g. a narrowing pipe, I will probably switch the project over to p5.js.

The project is still not entirely finished but incorporates methods to allow for the modelling the effects of a pressure gradient, viscosity and diffusion on the flow of the fluid.

It is still assumed, that the fluid is incompressible, viscosity and density are constant and that the process is isothermal.

Used resources: 

http://www.thevisualroom.com/index.html

https://matthias-research.github.io/pages/tenMinutePhysics/index.html

https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/17-fluidSim.html

https://john-s-butler-dit.github.io/NumericalAnalysisBook/Chapter%2009%20-%20Elliptic%20Equations/901_Poisson%20Equation-Laplacian.html

------------------To-Do-List--------------------------------------------

    //TODO Sliding obstacle over source "blocks" it, moving obstacle into source does not.
    //TODO Change all finite difference methods to central difference method
    //TODO Maybe make difference methods modular
    //TODO Decouple render speed from canvas size
    //TODO Add method to control and compare different numerical methods
    //TODO Let canvas size vary but sim size stays same, only visually stretched
    //TODO 3D with WebGl or p5.js
    //TODO Update and expand GUI
    //TODO decouple simulation speed from rendering speed
    //TODO Add dimensions to all physical units
    //TODO Make number of cells independent of sim cell size
    //TODO Check if dynamic viscosity and kinematic viscosity have not been falsly used.
    //! dynamic visocsity = kinematic viscosity * density. my = ny * rho
    //TODO solving incompressibillity causes overflow of float numbers. Like 3.5454...e300



