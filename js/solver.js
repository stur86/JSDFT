solver = {};

solver.grid_n = 40;
solver.mass = 1.0;

solver.x_scale = 1.0;      // Bohr
solver.V_scale = 10.0;     // Hartree

solver.Etol = 0.02;        // Hartree

solver.max_it = 20;

solver.init = function() {
    this.sol = null;    // A good way to keep track of the fact we haven't solved anything yet
    this.set_particlen(1);
    this.density = numeric.rep([this.grid_n], 0.0); // Starting density: empty
    this.x_step = this.x_scale*2.0/this.grid_n;
    this.V = numeric.rep([this.grid_n], 0.0);
    this.dftLOOP = null;
}

solver.build_potential = function() {

    // Rebuild total potential from current existing path

    this.V = [];

    var pnode = potential_plot.Vline.node();
    var path_l = pnode.getTotalLength();

    // Now travel it to find the potential
    grid_i = 0.0;
    for (var px = 0; px < path_l; ++px) {
        p_atL = pnode.getPointAtLength(px);
        var x = p_atL.y/potential_plot.w;
        var y = p_atL.x/potential_plot.h;
        while (x > grid_i/solver.grid_n) {
            this.V.push((y-0.5)*2.0*this.V_scale);
            ++grid_i;
        }
    }


}

solver.solve = function() {
    // Well, duh.

    // If it's ALREADY solving, get outta here

    if (this.dftLOOP != null) {
        return;
    }

    // Clear log

    logBox.clear();

    // Number of particles

    if (this.part_n == 1) {
        // Pretty simple 

        logBox.log("Solving exactly for 1 particle");

        // First, build potential
        this.build_potential();
        // Then the Hamiltonian
        var d = (this.x_scale/this.grid_n)*2.0;
        this.K = kin_op_gen(this.grid_n, this.mass, d, false);
        // Diagonalize
        this.sol = mnumerov_solve(solver.grid_n, this.K, this.V);

        // And build the density
        this.density = numeric.pow(this.sol[0].evec, 2.0);

        // Finally, replot
        wavef_plot.redraw();        

        logBox.log("Solution complete");
        logBox.log("Energy: " + this.sol[0].eval + " Ha");

    }
    else {
        // DFT is required

        logBox.log("Initializing DFT");

        // First, build potential
        this.build_potential();
        // Then the Hamiltonian
        this.K = kin_op_gen(this.grid_n, this.mass, this.x_step, false);
        this.sol = mnumerov_solve(this.grid_n, this.K, this.V);

        // Build density

        this.E_0 = 0.0;
        this.density = numeric.rep([this.grid_n], 0.0);
        for (var i = 0; i < this.part_n; ++i) {
            this.E_0 += this.sol[i].eval;
            numeric.addeq(this.density, numeric.pow(this.sol[i].evec, 2.0));
        }

        logBox.log("Wavefunction initialized, commencing DFT convergence...");

        this.E_1 = NaN;
        this.dE = 0.0;

        this.it = 0;
        
        dft_iteration = function() {

            if (solver.is_computing) {
                return;
            }

            solver.is_computing = true;

            // Redefine potential
            V_dft = numeric.clone(solver.V);
            for (var i = 0; i < solver.grid_n; ++i)
            {
                for (var j = 0; j < solver.grid_n; ++j)
                {
                    if (j == i)
                    {
                        continue;
                    }

                    V_dft[i] += solver.density[i]*solver.density[j]/(Math.abs(j-i)*solver.x_step);   // Coulombic term
                    // Self interaction correction
                    for (var k = 0; k < solver.part_n; ++k) {
                        V_dft[i] -= Math.pow(solver.sol[k].evec[i], 2.0)*Math.pow(solver.sol[k].evec[j], 2.0)/(Math.abs(j-i)*solver.x_step);
                    }
                }

                V_dft[i] -= solver.part_n*GLDA.ec(solver.part_n, 1.0/(2.0*solver.density[i]))/1000.0;        // LDA correlation functional term
            }

            solver.sol = mnumerov_solve(solver.grid_n, solver.K, V_dft);

            // Build density

            solver.E_1 = 0.0;

            solver.density = numeric.rep([solver.grid_n], 0.0);
            for (var i = 0; i < solver.part_n; ++i) {
                solver.E_1 += solver.sol[i].eval;
                numeric.addeq(solver.density, numeric.pow(solver.sol[i].evec, 2.0));
            }

            ++solver.it;

            solver.dE = Math.abs(solver.E_1-solver.E_0);

            logBox.log("Iteration n. " + solver.it + " completed - E = " + solver.E_1 + " - dE = " + solver.dE);

            solver.E_0 = solver.E_1;

            dft_completed = (solver.it >= solver.max_it || solver.dE <= solver.Etol);
            wavef_plot.redraw();

            if (dft_completed)
            {
                // Completing actions
                clearInterval(solver.dftLOOP);
                solver.dftLOOP = null;
                logBox.log("DFT convergence completed after " + solver.it + " iterations");
                if (solver.it >= solver.max_it && solver.dE > solver.Etol)
                {
                    logBox.log("Energy did not converge in the given limit - try increasing the number of iterations or the tolerance");
                }
                else
                {
                    logBox.log("Final energy: " + solver.E_1 + " Ha");
                }
            }

            solver.is_computing = false;
        }

        solver.is_computing = false;
        solver.dftLOOP = setInterval(dft_iteration, 100);

    }


}

solver.set_particlen = function(n) {

    $("#current_pn").html(n);
    solver.part_n = n;

}

solver.set_particlemass = function(m) {

   solver.mass = m;

}