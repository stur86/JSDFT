var SliderControl = function(id, minVal, maxVal, updateMethod) {

	// First build the drag behaviour

    this.drag_behav = d3.behavior.drag()
	  .on("drag", function(d) {
	    var mywidth = $(this).width();
	    var parentwidth = $(this.parentNode).width();
	    var newleft = d3.event.x - mywidth/2.0;
	    newleft = newleft < 0? 0 : (newleft > parentwidth-mywidth? parentwidth-mywidth : newleft);
	    d3.select(this).style('left', newleft);
	  })
	  .on("dragend", function(d) {
	    var mywidth = $(this).width();
	    var parentwidth = $(this.parentNode).width();
	    var frac = parseFloat(d3.select(this).style('left'))/(parentwidth-mywidth);
	  	updateMethod(minVal+frac*(maxVal-minVal));
	  });

	this.bar_d3 = d3.select('#' + id + '.scrollbar_css');
	this.handle_d3 = this.bar_d3.select('.scrollhandle_css');
	this.bar_d3.select('.scrolllabel_left').html(minVal);
	this.bar_d3.select('.scrolllabel_right').html(maxVal);
	this.handle_d3.call(this.drag_behav);

}