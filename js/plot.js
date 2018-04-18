// Handling actual plotting

wavef_plot = {};

wavef_plot.scale = [-0.05, 0.2];			// y scale

wavef_plot.init = function() {
	wavef_plot.svg = d3.select('#wavef_svg');

	// Clean up everything
	this.svg.html("");

	wavef_plot.redraw();
}

wavef_plot.update_scale = function(s) {

	wavef_plot.scale = [-s/4.0, s];
	wavef_plot.redraw();

}

wavef_plot.update_shadows = function() {

	// Get the size
	var w = $('#wavef_svg').width();
	var h = $('#wavef_svg').height();

	var dV = [];
	// Maximum possible derivative
	var max_dV = 2.0*solver.V_scale;

	if (get_style() == "billiards") {
		solver.build_potential();
		dV = numeric.sub(solver.V.slice(1, solver.grid_n), solver.V.slice(0, solver.grid_n-1));
	}

	var shadow_sel = this.svg.selectAll('.shadow_rect').data(dV);

	shadow_sel.enter().append('rect').classed('shadow_rect', true);

	shadow_sel.attr('x', function(d, i) { return Math.floor(i*w/(solver.grid_n-1.0))})
			.attr('y', 0)
			.attr('width', function(d, i) { return Math.floor((i+1)*w/(solver.grid_n-1.0))-Math.floor(i*w/(solver.grid_n-1.0)); })
			.attr('height', h)
			.attr('fill', function(d) { return (d>0?'#ccffcc':'#003300')})
			.attr('opacity', function(d) { return Math.abs(d)/(1.2*max_dV)});

	shadow_sel.exit().remove();
}

wavef_plot.redraw = function() {

	// Get the size
	var w = $('#wavef_svg').width();
	var h = $('#wavef_svg').height();

	// Define scales

	this.fracx = d3.scale.linear().domain([0.0, 1.0]).range([0.0, w]);
	this.fracy = d3.scale.linear().domain([0.0, 1.0]).range([h, 0.0]);

	// Fetch current style, act upon it

	switch (get_style()) {
		case "particles":
			// "Particles" style - the simplest one, but uses a grid
			// Draw grid lines

			var grid = this.svg.selectAll("#wfp_grid").data([0]);
			grid.enter().append('g').attr('id', 'wfp_grid');

			// Draw lines, evenly spaced
			xl = numeric.linspace(0, 1, 20);
			yl = numeric.linspace(0, 1, 12);

			this.gridx = grid.selectAll('.gridline_x').data(xl)

			this.gridx.enter().append('line').classed('gridline gridline_x', true);

			this.gridx.attr('x1', function(d) {return wavef_plot.fracx(d);})
				.attr('x2', function(d) {return wavef_plot.fracx(d);})
				.attr({'y1': 0, 'y2': h});

			this.gridx.exit().remove();

			this.gridy = grid.selectAll('.gridline_y').data(xl)

			this.gridy.enter().append('line').classed('gridline gridline_y', true);

			this.gridy.attr('y1', function(d) {return wavef_plot.fracy(d);})
				.attr('y2', function(d) {return wavef_plot.fracy(d);})
				.attr({'x1': 0, 'x2': w});

			this.gridy.exit().remove();

			// Create the axis

			var x_axis_scale = d3.scale.linear().domain([-0.8*solver.x_scale, 0.8*solver.x_scale]).range([0.1*w, 0.9*w]);
			var x_axis = this.svg.selectAll('.x.axis').data([0]);
			x_axis.enter().append('g')
				.classed('x axis', true)
				.attr('transform', 'translate(0,' + h*0.9 + ')');
			x_axis.call(d3.svg.axis().scale(x_axis_scale).ticks(4));

			x_label = this.svg.selectAll('.x.label').data(['x (Ang)']);
			x_label.enter().append('text').classed('x label', true);
			x_label.attr('transform', 'translate(' + w*0.5 + ',' + h*0.98 + ')')
				.attr('text-anchor', 'middle')
				.text(function(d) {return d;});

			// Now for vertical ones

			var scale_range = this.scale[1] - this.scale[0]
			var rho_axis_scale = d3.scale.linear().domain([this.scale[0]+scale_range*0.1, this.scale[0]+scale_range*0.9]).range([0.9*h, 0.1*h]);
			var rho_axis = this.svg.selectAll('.rho.axis').data([0]);
			rho_axis.enter().append('g')
				.classed('rho axis', true)
				.attr('transform', 'translate(' + w*0.05 + ',0)');
			var tickvals = [];
			for (var tval = 0.05; tval < this.scale[1]; tval += 0.05) {
				tickvals.push(tval);
			}
			rho_axis.call(d3.svg.axis().scale(rho_axis_scale).ticks(4).orient('left').tickValues(tickvals));

			rho_label = this.svg.selectAll('.rho.label').data(['Density']);
			rho_label.enter().append('text').classed('rho label', true);
			rho_label.attr('transform', 'translate(' + w*0.015 + ',' + h*0.5 + ') rotate(-90)')
				.attr('text-anchor', 'middle')
				.text(function(d) {return d;});

			break;

		case "billiards":

			var billiard_model = d3.select('#billiard_model').html();
			// Get the starting coordinates to correct for
			var billiard_bbox = d3.select('#billiard_model').select('g').node().getBBox();
			// First set the scale factor to fit four of these in the screen (with some margin)
			var billiard_scale = 0.20*h/billiard_bbox.height;
			// How many billiards?
			var billiard_n = Math.ceil(w/(billiard_bbox.width*billiard_scale));
			// Find the intensities for each billiard
			var billiard_data = []

			for (var i = 0; i < solver.part_n; ++i) {
				var tot_rho = 0.0
				for (var j = 0; j < billiard_n; ++j) {
					// Find the cumulative density for this one
					var rho_ij = 0.0;
					if (solver.sol != null)
						rho_ij = numeric.sum(numeric.pow(solver.sol[i].evec.slice(j*solver.grid_n/billiard_n, (j+1)*solver.grid_n/billiard_n), 2.0));
					tot_rho += rho_ij;
					//console.log(rho_ij);
					billiard_data.push({
						'i': i,
						'x': j*billiard_bbox.width,
						'y': i*billiard_bbox.height+(i*h*0.05+h*0.025)+(4-solver.part_n)*h/6.0,
						'rho': rho_ij,
					});
				}
			}

			var billiard_sel = this.svg.selectAll('.quantum_billiards').data(billiard_data);

			billiard_sel.enter().append('g').classed('quantum_billiards', true).attr('opacity', 0.0);

			billiard_sel.html(billiard_model)
				.transition()
				.attr('transform', function(d) { 
					return 'scale(' + billiard_scale + ') translate (' + (-billiard_bbox.x+d.x) + ',' +  + (-billiard_bbox.y+d.y) + ')';
				})
				.attr('opacity', function(d) {
					return Math.sqrt(d.rho);
				});

			billiard_sel.select('#ball_numtext').html(function(d) { return d.i+1;});
			billiard_sel.select('#ball_colorbar').attr("class", function(d) { return 'color_' + d.i;});

			billiard_sel.exit().transition().attr('opacity', 0.0).remove();

			break;
	}

	// Now draw the density line, common to everyone

	rholine_func = d3.svg.line()
			.x(function(d) {return wavef_plot.fracx(d.x)})
			.y(function(d) {return wavef_plot.fracy(d.y)})
			.interpolate('monotone');

	// Then to build the actual data to feed

	var rhodata = [];
	var scale_range = this.scale[1] - this.scale[0];

	for (var i = 0; i < solver.grid_n; ++i) {
		rhodata.push({'x': i/solver.grid_n, 'y': (solver.density[i]-wavef_plot.scale[0])/scale_range});
	}
	// Add a final point
	rhodata.push({'x': 1.0, 'y': (solver.density[0]-wavef_plot.scale[0])/scale_range});

	this.rholine = this.svg.select(".rholine");

	if (this.rholine.empty()) {
		this.rholine = this.svg.append("path")
					.attr("d", rholine_func(rhodata))
					.classed("rholine", true);
	}
	else {
		this.rholine.transition().duration(500).attr("d", rholine_func(rhodata));
	}


}

// Potential plot 

potential_plot = {};

potential_plot.init = function() {

	potential_plot.svg = d3.select("#potential_svg");
	potential_plot.points = [{'x': 0.0, 'y': 0.9, 'id': "startp"}, {'x': 1.0, 'y': 0.9, 'id': "endp"}];			// A series of xy pairs defining the potential

	// Clean up everything
	this.svg.html("");

	potential_plot.id_count = 0;
	potential_plot.redraw();

	// events

	this.svg.on("mousedown", function(d, i) {

		// Are we hovering or not?
		is_hover = !potential_plot.svg.selectAll(".potdot").filter(".hovered").empty();
		if (!is_hover || d3.event.which != 1) {
			// Left or right click?
			switch(d3.event.which) {
				case 1:
					potential_plot.add_point(d3.mouse(this));
					break;
				case 3:
					potential_plot.remove_point(d3.mouse(this));
					break;
			}
		}

	});

	// We also need to prevent default behaviour for right click
	this.svg.on("contextmenu", function(data, index) {
	     //handle right click
		//potential_plot.remove_point(d3.mouse(this));	     
	     //stop from showing browser menu
	    d3.event.preventDefault();
	});

	// And generate a drag behaviour for all dots

	this.drag_behav = d3.behavior.drag()
		.on("dragstart", function(d) {
			// Find the current index of this in the points
			thisd3 = d3.select(this);
			thisid = thisd3.attr("id");
			for (var i = 1; i < potential_plot.points.length-1; ++i) {
				if (thisid == potential_plot.points[i].id) {
					thisd3.attr("data-p_i", i);
					break;
				}
			}
		})
		.on("drag", function(d) {
			thisd3 = d3.select(this);
			thisd3.attr({'cy': d3.event.y});
			potential_plot.points[thisd3.attr("data-p_i")].y = 1.0-d3.event.y/potential_plot.h;
			potential_plot.redraw(true);
		});

}

potential_plot.redraw = function (fast) {

	fast = fast || false;

	// Get the size
	this.w = $('#potential_svg').width();
	this.h = $('#potential_svg').height();

	// Define scales

	this.fracx = d3.scale.linear().domain([0.0, 1.0]).range([0.0, this.w]);
	this.fracy = d3.scale.linear().domain([0.0, 1.0]).range([this.h, 0.0]);


	switch (get_style()) {
		case "particles":
			// "Particles" style - the simplest one

			// Axis

			var V_axis_scale = d3.scale.linear().domain([-0.8*solver.V_scale, 0.8*solver.V_scale]).range([0.9*this.h, 0.1*this.h]);
			var V_axis = this.svg.selectAll('.V.axis').data([0]);
			V_axis.enter().append('g')
				.classed('V axis', true)
				.attr('transform', 'translate(' + this.w*0.05 + ',0)');
			V_axis.call(d3.svg.axis().scale(V_axis_scale).ticks(4).orient('left'));

			V_label = this.svg.selectAll('.V.label').data(['Potential (Ha)']);
			V_label.enter().append('text').classed('rho label', true);
			V_label.attr('transform', 'translate(' + this.w*0.015 + ',' + this.h*0.5 + ') rotate(-90)')
				.attr('text-anchor', 'middle')
				.text(function(d) {return d;});

			break;
	}

	Vline_func = d3.svg.line()
						.x(function(d) {return potential_plot.fracy(1.0-d.y)})
						.y(function(d) {return potential_plot.fracx(d.x)})
						.interpolate('cardinal')
						.tension(0.85);

	this.Vline = this.svg.select(".potline");

	if (this.Vline.empty()) {
		this.Vline = this.svg.insert("path", ":first-child")
					.attr("transform", "rotate(-90," + this.fracy(0.5) + "," + this.fracy(0.5) + ")")
					.attr("d", Vline_func(this.points))
					.classed("potline", true);
	}
	else {
		if (!fast) {
			this.Vline.transition().duration(500).attr("d", Vline_func(this.points)).each("end", function() {
				wavef_plot.update_shadows()
			});
		}
		else {
			this.Vline.attr("d", Vline_func(this.points));			
			wavef_plot.update_shadows();
		}
	}


}

potential_plot.add_point = function(xy) {

	xy_frac = [xy[0]/this.w, 1.0 - xy[1]/this.h];

	// Add a potential dot

	// First add it in the array
	for (var i = 0; i < this.points.length; ++i) {
		if (this.points[i].x > xy_frac[0]) {
			this.points.splice(i, 0, {'x': xy_frac[0], 'y': xy_frac[1], 'id': "potdot_" + this.id_count});
			break;
		}
	}

	var newdot = this.svg.append("circle")
		.classed("potdot", true)
		.attr({	"id": "potdot_" + this.id_count,
			   	"cx": xy[0], "cy": xy[1], "r": 10});

	newdot.on("mouseenter", function() {d3.select(this).classed("hovered", true)});
	newdot.on("mouseleave", function() {d3.select(this).classed("hovered", false)});
	newdot.call(this.drag_behav);

	this.id_count += 1;

	this.redraw();

}

potential_plot.remove_point = function(xy) {

	// Check if it's hovering on a point, and in case remove it
	var to_remove = this.svg.selectAll(".potdot").filter(".hovered");

	if (to_remove.empty()) {
		return;
	}
	// Remove the point
	to_remove_id = to_remove.attr("id");
	for (var i = 0; i < this.points.length; ++i) {
		if (this.points[i].id == to_remove_id) {
			this.points.splice(i, 1);
			break;
		}
	}
	to_remove.remove();

	this.redraw();

}

logBox = {}

logBox.clear = function() {
    $('.log_box').html('');
}

logBox.log = function(s) {

    $('.log_box').html($('.log_box').html() + s + '<br>');
	$('.log_box').scrollTop(1e10);
}