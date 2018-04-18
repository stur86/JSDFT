// Solving the Schroedinger equation with a matrix Numerov algorithm 

function invB_gen(n, period)
{
	var theta = [];
	var phi = [];
	var B  = numeric.mul(10.0/12.0, numeric.identity(n));

	var inv12 = 1.0/12.0;

	for (var j = 0; j < 2; ++j)
	{
		var j_s = j*2-1;

		B[1][0] = inv12;
		B[n-2][n-1] = inv12;

		if (period)
		{
			B[n-1][0] = inv12;
			B[0][n-1] = inv12;
		}

		for (var i = 1; i < n-1; ++i)
		{
			B[i+j_s][i] = inv12;
		}
	}

	invB = numeric.inv(B);

	return invB;

}

function A_gen(n, d, period)
{
	var invd2 = 1.0/(d*d);
	var m2invd2 = -2.0/(d*d);

	var A = numeric.mul(m2invd2, numeric.identity(n));

	for (var j = 0; j < 2; ++j)
	{
		var j_s = j*2-1;

		A[1][0] = invd2;
		A[n-2][n-1] = invd2;

		if (period)
		{
			A[n-1][0] = invd2;
			A[0][n-1] = invd2;
		}

		for (var i = 1; i < n-1; ++i)
		{
			A[i+j_s][i] = invd2;
		}
	}

	return A;
}

function kin_op_gen(n, m, d, period)
{

	var invB = invB_gen(n, period);
	var A    = A_gen(n, d, period);

	K = numeric.mul(-1.0/(2.0*m), numeric.dot(invB, A));	// Using atomic units

	return K;

}

function hamiltonian_op(n, K, V)
{
	H = numeric.clone(K);

	for (var i = 0; i < n; ++i)
	{
		H[i][i] += V[i];
	}

	return H;
}

function mnumerov_solve(n, K, V)
{

	H = hamiltonian_op(n, K, V);
	res = numeric.eig(H);

	sol = [];

	evec_mat = numeric.transpose(res.E.x);

	for (var i = 0; i < n; ++i)
	{
		sol.push({'eval': res.lambda.x[i], 'evec': evec_mat[i]});
	}

	// Now sort them

	sol.sort(function(a,b) {return a.eval-b.eval;});

	return sol;

}