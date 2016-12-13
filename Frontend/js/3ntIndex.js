			//Width and height
			var marginFigure1 = {top: 20, right: 20, bottom: 50, left: 100},
    			widthFigure1 = 700 - marginFigure1.left - marginFigure1.right,
    			paddingFigure1=10,
    			heightFigure1 = 250 - marginFigure1.top - marginFigure1.bottom;
    			
    			
    		var bisectposition = d3.bisector(function(d) { return d.position; }).right;

    		
    		//Create scale functions	
    		var xFigure1 = d3.scale.linear()
    						.range([paddingFigure1, widthFigure1 - paddingFigure1 * 2]);

			var yFigure1 = d3.scale.linear()
    						.range([heightFigure1-paddingFigure1*2, paddingFigure1]);
    			
			
			//Define X axis
			var xAxisFigure1 = d3.svg.axis()
    			.scale(xFigure1)
    			.orient("bottom")
    			.ticks(11);
			
			//Define Y axis
			var yAxisFigure1 = d3.svg.axis()
    			.scale(yFigure1)
    			.orient("left")
    			.ticks(6);
    		
			
			//Define line
			var valuelineFigure1 = d3.svg.line()
    			.x(function(d) { return xFigure1(d.position); })
    			.y(function(d) { return yFigure1(d.count); });
    			
    		
			
			//Create SVG element
			var svgFigure1 = d3.select("#nucleotidePeriodicity")
						.append("svg")
    					.attr("width", widthFigure1 + marginFigure1.left + marginFigure1.right)
    					.attr("height", heightFigure1 + marginFigure1.top + marginFigure1.bottom)
  						.append("g")
    					.attr("transform", "translate(" + marginFigure1.left + "," + marginFigure1.top + ")");
			
    
	d3.tsv("../Data/F1_2016_Weinberg_RPF.tsv", function(error, data) {
  		
  		if (error) throw error;
  		data.forEach(function(d) {
    		d.position = +d.Position;
    		d.count = +d.Counts;
    		d.end=+d.End;
  		});

  
				//x-axis	
  				svgFigure1.append("g")
      				.attr("class", "x axis")
      				.attr("transform", "translate(0," + (heightFigure1) + ")")
      				.call(xAxisFigure1);


				//x-axis line
    			svgFigure1.selectAll(".x.axis")	
  					.append("rect")
  	   				.attr("width", widthFigure1)
  	   				.attr("height",1)
  	   				.attr("fill","#000");
				
				//y-axis
  				svgFigure1.append("g")
      				.attr("class", "y axis")
      				.call(yAxisFigure1);
      				
      				// now add titles to the axes
        svgFigure1.append("text")
            .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
            .attr("transform", "translate("+(0-paddingFigure1*6)+","+(heightFigure1/2)+")rotate(-90)")  // text is drawn off the screen top left, move down and out and rotate
            .text("Mapped reads").style("font-size","16px");

       svgFigure1.append("text")
            .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
            .attr("transform", "translate("+ (widthFigure1/2) +","+(heightFigure1+paddingFigure1*3.6)+")")  // centre below axis
            .text("Distance from translation start/stop (nt)").style("font-size","16px");
      				  				
    				
function updateFigure1(value4) {
      				
var datanew = data.filter(function(d, key) { 
            return (d.end==value4);
    
    });
  
xFigure1.domain(d3.extent(datanew, function(d) { return d.position; }));
yFigure1.domain([0, d3.max(datanew, function(d) { return d.count; })]);


var freqFigure1 = svgFigure1.selectAll(".datanew")
      		.data(datanew);
      		
freqFigure1.select("g .datanew path") 
						.transition()
						 .ease("linear")
						.duration(600)
						.attr("class", "line")
						.attr("d", valuelineFigure1(datanew));

freqFigure1.select("g .datanew circle")
    				.transition()
             		.ease("linear")
					   .delay(function(d, i) {
						   return i / data.length * 1000;
					   })
					.duration(600)
        			.attr("r", 4.5)
        			.attr("cx", function(d) { return xFigure1(d.position); })
        			.attr("cy", function(d) { return yFigure1(d.count); })
        			.style("fill",function(d) { return d.color; })
        			.style("fill","rgb(239,138,98)");

    				
var freqgroup=freqFigure1.enter().append("g").attr("class", "datanew"); 
					
					freqgroup.append("path")
             		.attr("class", "line")
             		.style("stroke-opacity", 0.1)
             		.attr("d", valuelineFigure1(datanew));
             		
             		freqgroup.append("circle")
             		.attr("r", 4.5)
        			.attr("cx", function(d) { return xFigure1(d.position); })
        			.attr("cy", function(d) { return yFigure1(d.count); })
        			.style("fill",function(d) { return d.color; })
        			.style("fill","rgb(239,138,98)");
      		 

        			
var focus = freqFigure1.append("g") 
    				.style("display", "none");

				// append the x line on mouse
    			focus.append("line")
        			.attr("class", "x")
        			.style("stroke", "rgb(239,138,98)")
        			.style("stroke-dasharray", "3,3")
        			.style("opacity", 0.75)
        			.attr("y1", 0)
        			.attr("y2", heightFigure1);

    			// append the y line on mouse
    			focus.append("line")
        			.attr("class", "y")
        			.style("stroke", "rgb(239,138,98)")
        			.style("stroke-dasharray", "3,3")
        			.style("opacity", 0.75)
        			.attr("x1", widthFigure1)
        			.attr("x2", widthFigure1);

    			// append the circle at the intersection
    			focus.append("circle")
        			.attr("class", "y")
        			.style("fill", "rgb(239,138,98)")
        			.style("stroke", "rgb(239,138,98)")
        			.attr("r", 5.2);
            
    			// append the rectangle to capture mouse
    			freqFigure1.append("rect")
        			.attr("width", widthFigure1)
        			.attr("height", heightFigure1)
        			.style("fill", "none")
        			.style("pointer-events", "all")
        			.on("mouseover", function() { focus.style("display", null); })
        			.on("mouseout", function() { focus.style("display", "none"); })
        			.on("mousemove", mousemove);

 		 
freqgroup.transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 1000;
			})
      		.duration(600);
      		
freqFigure1.exit()
  .transition()
    .ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 1000;
			})
      .duration(600)
      .style("fill-opacity", 1)
      .remove();     		
      		      		      		

			
			
						//Upposition X axis
						svgFigure1.select(".x.axis")
							.transition()
							.ease("linear")
							.delay(function(d, i) {
								return i / datanew.length * 500;
							})
							.duration(600)
							.call(xAxisFigure1);
				
						//Upposition Y axis
						svgFigure1.select(".y.axis")
							.transition()
							.ease("linear")
							.delay(function(d, i) {
								return i / datanew.length * 500;
							})
							.duration(600)
							.call(yAxisFigure1);

function mousemove() {
					var x0 = xFigure1.invert(d3.mouse(this)[0]),
		    		i = bisectposition(datanew, x0, 1),
		    		d0 = datanew[i - 1],
		    		d1 = datanew[i],
		    		d = x0 - d0.position > d1.position - x0 ? d1 : d0;

					focus.select("circle.y")
		    		.attr("transform",
		          		"translate(" + xFigure1(d.position) + "," +
		                         yFigure1(d.count) + ")");

					focus.select("text.y1")
		    			.attr("transform",
		          "translate(" + xFigure1(d.position) + "," +
		                         yFigure1(d.count) + ")")
		    			.text(d.count);

					focus.select("text.y2")
		    			.attr("transform",
		          "translate(" + xFigure1(d.position) + "," +
		                         yFigure1(d.count) + ")")
		    			.text(d.count);

					focus.select("text.y3")
		    			.attr("transform",
		          "translate(" + xFigure1(d.position) + "," +
		                         yFigure1(d.count) + ")")
		    			.text(d.position);

					focus.select("text.y4")
		    			.attr("transform",
		          "translate(" + xFigure1(d.position) + "," +
		                         yFigure1(d.count) + ")")
		    			.text(d.position);

					focus.select(".x")
		    			.attr("transform",
		          "translate(" + xFigure1(d.position) + "," +
		                         yFigure1(d.count) + ")")
		               .attr("y2", heightFigure1 - yFigure1(d.count));

					focus.select(".y")
		    			.attr("transform",
		          "translate(" + widthFigure1 * -1 + "," +
		                         yFigure1(d.count) + ")")
		               .attr("x2", widthFigure1 + widthFigure1);
				}

};

	updateFigure1(5);
	
	d3.select(document.getElementById("3prime"))
	.on("click", function() {
		updateFigure1(3);
	});
	d3.select(document.getElementById("5prime"))
	.on("click", function() {
		updateFigure1(5);
	});


});




