//! viscosity = const , means we're dealing only with newtonian fluids because I'm not a mad man

function computeLaplacian(field) { //Compute the lapalcian of a 2D VECTOR-FIELD! with second-order central difference method. Atm. vertical and horizontal component must be stored in separate 1D arrays
    //Laplacion of a function f(x,) with a step size/grid spacing of h. Note that a scalar field, like the pressure Field has no laplacian associated with it.
    //Lap_h_f(x,y)=( f(x-h, y) + f(x+h,y) + f(x,y-h) + f(x,y+h) - 4*f(x,y) )/(h*h)
    var n = this.numYCells;
    var laplacian = new Float32Array(this.numCells);
    for (var i = 1; i < this.numXCells - 1; i++) {
        for (var j = 1; j < this.numYCells - 1; j++) {
            var index = i * n + j;

            laplacian[index] = (field[(i - 1) * n + j] + field[(i + 1) * n + j] + field[i * n + j - 1] + field[i * n + j + 1] - 4 * field[i * n + j]) / (this.simCellSize * this.simCellSize)
        }

    }
    return laplacian;
}



/*function computeAdvectiveTerm(u, v, field) {
//Optional: For a diffrenct approach where first every change to the velocity field is calculated and the velocity field than updated, but not advected in the sense of a flowing fluid.

    // Compute advection of the Scalar-Field 'Field' by the velocity field u and v
    const n = this.numYCells;
    const result = new Array(n * n).fill(0);
    for (let i = 1; i < n - 1; i++) {
        for (let j = 1; j < n - 1; j++) {
            const index = i * n + j;

            //Determine Gradient of the velocity field u,v
            //Aproximate partial derivatives by central difference.
            let du_dx = (u[j*n + (i+1)] - u[j*n + (i-1)]) / 2*this.simCellSize;
            let dv_dy = (v[(j+1)*n + i] - v[(j-1)*n + i]) / 2*this.simCellSize;
            let gradAtPoint = [du_dx , dv_dy];

            //calculated dot-product of the gradient with the field 'field'
            for (let k = 0; k < 3; k++)
            {

            }

        }
    }
    return result;

}
*/

function computePressureTerm(p, isHorizontal) { //compute the influence of the pressure gradient on the velocity field by a forward finite difference method.
    // pressureDifference = p(x+h)-p(x) ) / h
    // pressureTermV(i,j) = pressureDifference(i,j) * dt / (density * cellSize);
    for (var i = 1; i < this.numXCells - 1; i++) {
        for (var j = 1; j < this.numYCells - 1; j++) {
            var index = i * n + j;

            if(isHorizontal = true){
                //Compute horizontal pressure term
                var pressureDifferenceU = (this.pressureField[i * n + j]-this.pressureField[(i + 1) * n + j]) / this.simCellSize;
                var pressureTermU = pressureDifferenceU * dt / (density * this.simCellSize);
                pressureTermU[index] = pressureTermU;
                //Note: [(i + 1) * n + j] is right of [i * n + j]
            }
            else {
                //Compute vertical pressure term
                var pressureDifferenceV = (this.pressureField[i * n + j]-this.pressureField[i * n + (j + 1)]) / this.simCellSize;
                var pressureTermV = pressureDifferenceV * dt / (density * this.simCellSize);
                pressureTermV[index] = pressureTermV;
                //Note: [i * n + (j + 1)] is visually "below" [i * n + j]. Note that the y-direction of the grid points downwards
            }
        }
    }
    return (isHorizontal ? pressureTermU : pressureTermV);
}

// Compute the viscosity term
const laplacianU = computeLaplacian(u);
const laplacianV = computeLaplacian(v);
const viscousTermU = viscosity * laplacianU;
const viscousTermV = viscosity * laplacianV;

//Compute current divergence. Total amount of Out-/Inflow of fluid
function computeDivergence(horizontal, vertical) {
    for (var i = 1; i < this.numXCells - 1; i++) {
        for (var j = 1; j < this.numYCells - 1; j++) {
            var index = i * n + j;
            var divergenceValue = (horizontal[(i + 1) * n + j] - horizontal[i * n + j] + vertical[i * n + j + 1] - vertical[i * n + j]) / this.simCellSize;
            divergence[index] = divergenceValue                             //TODO innitialize field...somwhere
        }
    }
    return divergence;
}

function solvePressureField(pressureField, u, v, numIterations) {
    //TODO add static pressure
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
                var numberFreeEdges = sx0 + sx1 + sy0 + sy1;		                            // sum up all values of s. This equals the Number of edges that are NOT walls.
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

                var laplacianUAtIndex = laplacianU[index];
                var laplacianVAtIndex = laplacianV[index];

                // Diffusion term added for viscosity
                this.xVel[i * n + j] += diffusionRate * laplacianUAtIndex - sx0 * relativeDivergence;
                this.xVel[(i + 1) * n + j] += diffusionRate * laplacianUAtIndex + sx1 * relativeDivergence;
                this.yVel[i * n + j] += diffusionRate * laplacianVAtIndex - sy0 * relativeDivergence;
                this.yVel[i * n + j + 1] += diffusionRate * laplacianVAtIndex + sy1 * relativeDivergence;
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
    //Check Cell Type
    //if in obstacle return true
    //if not return false
}

function changeCellType(obstacleType){ //Addition to setobstacle function
    //Change cell type to walls for all coordinates within the outline of an obstacle
    //Get coordinates of obstacles, determine outline of obstacle
    //go through all cells and change cell-type to wall if within outline
    //if statement dependent on obstacleType
    //Convert obstacle coordinates to sim coordinates
}

function drawObstacle(obstacleID){ //TODO Obstacles should not override each other

    //Draws a circle
    c.beginPath();
    c.arc(convertSimToCanvasX(scene.obstacleX), convertSimToCanvasY(scene.obstacleY), canvasScale * radius, 0.0, 2.0 * Math.PI);
    c.closePath();
    c.fill()
    c.lineWidth = 3.0; //!Unnecessary
    c.strokeStyle = "#000000";          //! Make a black edge. cool
    c.beginPath();
    c.arc(convertSimToCanvasX(scene.obstacleX), convertSimToCanvasY(scene.obstacleY), canvasScale * radius, 0.0, 2.0 * Math.PI);
    c.closePath();
    c.stroke();
    c.lineWidth = 1.0;  //!Unnecessary

    //Draw a rectangle
}
