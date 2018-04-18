GLDA = {};
// Reduced version, good enough for now
GLDA.par_table = {"2": [-0.18425928344365855, -2.2877442354587085, 370.17982716632997, -6.6748200888549096, 49.688907055374557, -5.0147203939002187, 10.321653580246998], 
				  "3": [-0.2418829110238255, -3.1097706142339319, 506.37025876331887, -9.1092393300183492, 66.18709647090148, -6.8913910929527518, 13.428196857289967], 
				  "4": [-0.26875129182764079, -3.5166606043007915, 642.26358327816911, -10.302049066620755, 82.750946514863344, -7.811627272522081, 16.559268501341514]};

/*
Discarded, not working if not from http
$.getJSON("data/LDA_3exp_pars.json", function(data) {
	GLDA.par_table = data;
})
*/

GLDA.gamma = 0.5772156649;	// Euler-Mascheroni constant

GLDA.ex = function(n, rs) {

	e_m2 = (n*n-1.0)/(n*n)*Math.pow(Math.PI, 2.0)/24.0;
	e_m1 = (n-1.0)/n*Math.log()

}

GLDA.ec = function(n, rs) {

	if (this.par_table[n] == null) {
		// Invalid n
		return 0.0;
	}

	E = this.par_table[n][0];
	L = 2.0*n*rs
	for (var i = 1; i < 4; ++i) {
		E += this.par_table[n][2*i-1]*Math.exp(-L/this.par_table[n][2*i]);
	}

	return E;
}