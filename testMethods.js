//! viscosity = const , means we're dealing only with newtonian fluids because I'm not a mad man

function computeLaplacian(velocity, dx) { //with second-order central difference method
    //Laplacion of a function f(x,) with a step size/grid spacing of h.
    //Lap_h_f(x,y)=( f(x-h, y) + f(x+h,y) + f(x,y-h) + f(x,y+h) - 4*f(x,y) )/(h*h)
}

// Compute the viscosity term
const du = computeLaplacian(u);
const dv = computeLaplacian(v);
const viscousTermU = viscosity * du;
const viscousTermV = viscosity * dv;

function computeAdvection(u, v, q) {
    // Compute advection of the Scalar-Field q by the velocity field u and v
}

function computeDivergence() {
    //Compute divergence of a field at a point i,j
    //Alt compute divergence
}

function solvePressureField(pressureField, u, v, numIterations) { //isU is boolean, pressure is determined differently if the horizontal component is concerned
    var n = this.numY;
    var alpha = -dt / (this.density*this.cellSize * this.cellSize);  // rate of change of pressure over time due to the diffusion of the fluid
    var beta = 4.0;                                                  // Coefficient to increase more heavily weigh the coeffSum, for stability and convergence reasons.

    for (var iter = 0; iter < numIterations; iter++) {
        for (var i = 1; i < this.numX - 1; i++) {
            for (var j = 1; j < this.numY - 1; j++) {
                var index = i * n + j;
                var divergence = computeDivergence(u,v);

                var sx0 = this.s[(i - 1) * n + j];
                var sx1 = this.s[(i + 1) * n + j];
                var sy0 = this.s[i * n + j - 1];
                var sy1 = this.s[i * n + j + 1];
                var numberFreeEdges = sx0 + sx1 + sy0 + sy1;		// sum up all values of s. This equals the Number of edges that are NOT walls.
                if (numberFreeEdges == 0.0 || this.s[index] == 0.0)						// If the cell itself or all surrounding cells ar walls
                {
                    continue;
                }

                var scaledDivergence = alpha * divergence[index]; //Determines influence of the divergence on the pressure field (right hand sight of the poisson equation)
                var coeffSum = pressureField[(i - 1) * n + j] + pressureField[(i + 1) * n + j] + pressureField[i * n + (j - 1)] + pressureField[i * n + (j + 1)]; //Sum up all surrounding pressure values

                var newPressure = (scaledDivergence + beta * coeffSum) / (beta * numberFreeEdges);
                pressureField[index] = newPressure;        //Set the newPressure at the index.
            }
        }
    }
}

function computePressureTerm(p, divergence, true);


function updateVelocity() {             //TODO change to 1D Arrays, update starting point for boundaries
    for (let i = 0; i < N; i++) {
        for (let j = 0; j < N; j++) {
            u[i][j] += dt * (-advectiveTermU[i][j] + viscousTermU[i][j] + pressureTermU[i][j]);
            v[i][j] += dt * (-advectiveTermV[i][j] + viscousTermV[i][j] + pressureTermV[i][j]);
        }
    }
}

function cellInObstacle(x,y){
    //if in obstacle return true
    //if not return false
}
