//! viscosity = const , means we're dealing only with newtonian fluids because I'm not a mad man

function computeLaplacian() { //with second-order central difference method
    //Laplacion of a function f(x,) with a step size/grid spacing of h.
    //Lap_h_f(x,y)=( f(x-h, y) + f(x+h,y) + f(x,y-h) + f(x,y+h) - 4*f(x,y) )/(h*h)

    for (var i = 1; i < this.numXCells - 1; i++) {
        for (var j = 1; j < this.numYCells - 1; j++) {
            var index = i * n + j;

            laplacianHorizontal[index] = (this.u[(i - 1) * n + j] + this.u[(i + 1) * n + j] + this.u[i * n + j - 1] + this.u[i * n + j + 1] - 4 * this.u[i * n + j]) / (this.h * this.h);
            laplacianVertical[index] = (this.v[(i - 1) * n + j] + this.v[(i + 1) * n + j] + this.v[i * n + j - 1] + this.v[i * n + j + 1] - 4 * this.v[i * n + j]) / (this.h * this.h);
        }
        return laplacianHorizontal;
        return laplacianVertical;
    }
}



function computeAdvectiveTerm(u, v, q) {
    // Compute advection of the Scalar-Field q by the velocity field u and v

    const advectiveTermU = computeAdvection(u, v, u);
    const advectiveTermV = computeAdvection(u, v, v);
}

function computePressureTerm(p, divergence,isHorizontal) {
    // pressureDiffrence = p(x+h)-p(x) ) / h
   // pressureTermV(i,j) = pressureDifference(i,j) * dt / (density * cellSize);
}

// Compute the viscosity term
const du = computeLaplacian(u);
const dv = computeLaplacian(v);
const viscousTermU = viscosity * du;
const viscousTermV = viscosity * dv;

//Compute current divergence. Total amount of Out-/Inflow of fluid
function computeDivergence() {
    for (var i = 1; i < this.numXCells - 1; i++) {
        for (var j = 1; j < this.numYCells - 1; j++) {
            var index = i * n + j;
            var divergenceValue = (this.horizontalVelocity[(i + 1) * n + j] - this.horizontalVelocity[i * n + j] + this.verticalVelocity[i * n + j + 1] - this.verticalVelocity[i * n + j]) / this.simCellSize;
            divergence[index] = divergenceValue                             //TODO innitialize field...somwhere
        }
    }
    return divergence;
}

function solvePressureField(pressureField, u, v, numIterations) {
    var n = this.numYCells;
    var alpha = -dt / (this.density * this.simCellSizecellSize * this.simCellSizecellSize);  // rate of change of pressure over time due to the diffusion of the fluid
    var beta = 4.0;                                                  // Coefficient to increase more heavily weigh the coeffSum, for stability and convergence reasons.

    for (var iter = 0; iter < numIterations; iter++) {
        for (var i = 1; i < this.numXCells - 1; i++) {
            for (var j = 1; j < this.numYCells - 1; j++) {
                var index = i * n + j;
                var div = divergence[index];                         //Calls the value of the divergence array at the index from the computeDivergence function

                var sx0 = this.cellType[(i - 1) * n + j];
                var sx1 = this.cellType[(i + 1) * n + j];
                var sy0 = this.cellType[i * n + j - 1];
                var sy1 = this.cellType[i * n + j + 1];
                var numberFreeEdges = sx0 + sx1 + sy0 + sy1;		// sum up all values of s. This equals the Number of edges that are NOT walls.
                if (numberFreeEdges == 0.0 || this.cellType[index] == 0.0)						// If the cell itself or all surrounding cells ar walls
                {
                    continue;
                }

                var scaledDivergence = alpha * div[index]; //Determines influence of the divergence on the pressure field (right hand sight of the poisson equation)
                var coeffSum = pressureField[(i - 1) * n + j] + pressureField[(i + 1) * n + j] + pressureField[i * n + (j - 1)] + pressureField[i * n + (j + 1)]; //Sum up all surrounding pressure values

                var newPressure = (scaledDivergence + beta * coeffSum) / (beta * numberFreeEdges);
                pressureField[index] = newPressure;        //Set the newPressure at the index.
            }
        }
    }
}



newforceIncompressibility(numIters, dt, viscosity)
{
    //Multiple iterations through all cells
    var n = this.numY;
    var currentPressure = densityValue * this.h / dt;
    var diffusionRate = viscosity * dt / (this.h * this.h);

    for (var iter = 0; iter < numIters; iter++) {    //Multiple iterations
        for (var i = 1; i < this.numX - 1; i++) {		// Go through all Cells
            for (var j = 1; j < this.numY - 1; j++) {
                var index = i * n + j;
                if (this.s[i * n + j] == 0.0)			// If there is a wall at s[i*n+j] go to the next cell
                {
                    continue;
                }

                var sx0 = this.s[(i - 1) * n + j];
                var sx1 = this.s[(i + 1) * n + j];
                var sy0 = this.s[i * n + j - 1];
                var sy1 = this.s[i * n + j + 1];
                var numberFreeEdges = sx0 + sx1 + sy0 + sy1;		// sum up all values of s.  This s then equals the Number of edges that are NOT walls.
                if (numberFreeEdges == 0.0)						// If alle cells are walls, go to the next cell
                {
                    continue;
                }

                var div = computeDivergence[index];                         //Calls the value of the divergence array at the index from the computeDivergence function

                var relativeDivergence = -div / numberFreeEdges;					// Equal part through each adjacent cells
                relativeDivergence *= scene.overRelaxation;
                this.pressureField[i * n + j] += currentPressure * relativeDivergence;

                // Diffusion term added for viscosity
                var laplacianU = (this.u[(i - 1) * n + j] + this.u[(i + 1) * n + j] + this.u[i * n + j - 1] + this.u[i * n + j + 1] - 4 * this.u[i * n + j]) / (this.h * this.h);
                var laplacianV = (this.v[(i - 1) * n + j] + this.v[(i + 1) * n + j] + this.v[i * n + j - 1] + this.v[i * n + j + 1] - 4 * this.v[i * n + j]) / (this.h * this.h);

                this.u[i * n + j] += diffusionRate * laplacianU - sx0 * relativeDivergence;
                this.u[(i + 1) * n + j] += diffusionRate * laplacianU + sx1 * relativeDivergence;
                this.v[i * n + j] += diffusionRate * laplacianV - sy0 * relativeDivergence;
                this.v[i * n + j + 1] += diffusionRate * laplacianV + sy1 * relativeDivergence;
            }
        }
    }
}



function updateVelocity() {             //TODO change to 1D Arrays, update starting point for boundaries



    for (let i = 0; i < N; i++) {
        for (let j = 0; j < N; j++) {
            index = i * n + j;
            u[index] += dt * (-advectiveTermU[index] + viscousTermU[index] + pressureTermU[index]);
            v[index] += dt * (-advectiveTermV[index] + viscousTermV[index] + pressureTermV[index]);
        }
    }
}

function cellInObstacle(x, y) {
    //if in obstacle return true
    //if not return false
}
