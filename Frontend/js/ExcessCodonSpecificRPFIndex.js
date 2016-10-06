            //Width and height
			var marginFigure6 = {top: 20, right: 10, bottom: 20, left: 100},
    			widthFigure6 = 700 - marginFigure6.left - marginFigure6.right,
    			paddingFigure6=5,
    			paddingFigure6l=50
    			heightFigure6 = 300 - marginFigure6.top - marginFigure6.bottom;
    			
    			
			var colorFigure6 = d3.scale.category20()
					.domain(["Data", "CHX", "FF"])
  					.range(["rgb(239,138,98)", "#a6cee3", "#1f78b4"]);
    			
    		var xFigure6 = d3.scale.ordinal()
    					.rangeRoundBands([paddingFigure6, widthFigure6 - paddingFigure6 * 2], .1);

			var yFigure6 = d3.scale.linear()
    						.range([heightFigure6-paddingFigure6*2, paddingFigure6]);
    			
			
			//Define X axis
			var xAxisFigure6 = d3.svg.axis()
    			.scale(xFigure6)
    			.orient("bottom");
			
			//Define Y axis
			var yAxisFigure6 = d3.svg.axis()
    			.scale(yFigure6)
    			.orient("left");
    		
			
			//Define line
			var valuelineFigure6 = d3.svg.line()
    			.x(function(d) { return xFigure6(d.Codon); })
    			.y(function(d) { return yFigure6(d.excess_reads); });
    			
    		
			
			//Create SVG element
			var svgFigure6 = d3.select("#ExcessCodonSpecificRPF")
						.append("svg")
    					.attr("width", widthFigure6 + 2*marginFigure6.left + marginFigure6.right)
    					.attr("height", heightFigure6 + marginFigure6.top + marginFigure6.bottom)
  						.append("g")
    					.attr("transform", "translate(" + marginFigure6.left + "," + marginFigure6.top + ")");


			//Load data
			d3.tsv("../../Data/F6_Year_2016_Author_Weinberg_Dataset_RPF_data.tsv", function(error, data) {
  				if (error) throw error;

  				data.forEach(function(d) {
  					d.Year = d.Year;
  					d.Author = d.Author;
  					d.Dataset = d.Dataset;
    				d.Codon = d.Codon;
    				d.excess_reads = +d.excess_reads;
    				d.CHX=d.CHX
    				d.SD=+d.SD
  				});

  				//x-axis	
  				svgFigure6.append("g")
      				.attr("class", "x axis")
      				.attr("transform", "translate(0," + (heightFigure6-3*paddingFigure6) + ")")
      			.call(xAxisFigure6)
      			.selectAll("text")
    				.attr("y", 0)
    				.attr("x", 10)
    				.attr("dy", ".35em")
    				.attr("transform", "rotate(90)")
    				.style("text-anchor", "start");
  
    			//x-axis line
    			svgFigure6.selectAll(".x.axis")	
  					.append("rect")
  	   				.attr("width", widthFigure6)
  	   				.attr("height",1)
  	   				.attr("fill","#000")
				
				//y-axis
  				svgFigure6.append("g")
      				.attr("class", "y axis")
      				.attr("transform", "translate(" + 0 + "," + (0-paddingFigure6) + ")")
      				.call(yAxisFigure6);
      				
      				
      				
function updateFigure6(data) {
      				
var datanew = data;
 
 

    			var nucleotide = d3.nest()
        			.key(function(d) {return d.CHX;})
        			.entries(datanew);
        
      
				legendSpace6 =  widthFigure6/nucleotide.length;
  				colorFigure6.domain(["Data", "CHX", "FF"]);
  

  				xFigure6.domain(datanew.map(function(d) { return d.Codon; }));
  				yFigure6.domain([d3.min(datanew, function(d) { return d.excess_reads; }),d3.max(datanew, function(d) { return d.excess_reads; })]);


var freqFigure6 = svgFigure6.selectAll(".nucleotide")
      		.data(nucleotide);

 freqFigure6.select("g .nucleotide path")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(750)
      		.attr("class", "line")
             		.attr("id", function(d) { return "tag"+d.key; }) // assign ID
             		.attr("d", function(d) { return valuelineFigure6(d.values); })
      		        .style("stroke", function(d) { return colorFigure6(d.key); })
      		        .filter(function(d) { return !isNaN(d.values)});;
      		
   freqFigure6.select("g .nucleotide text")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(500);
      		   		
      		      		
var lineSvg6=freqFigure6.enter().append("g").attr("class", "nucleotide");
        
				//Create the line
    			lineSvg6.append("path")
             		.attr("class", "line")
             		.attr("id", function(d) { return "tag"+d.key; }) // assign ID
             		.attr("d", function(d) { return valuelineFigure6(d.values); })
      		        .style("stroke", function(d) { return colorFigure6(d.key); });

 

      lineSvg6.append("text")
    			.datum(function(d) { return {name: d.key, value: d.values[d.values.length - 1]}; })
				.attr("transform", function(d, i) { return "translate(" + ((widthFigure6-paddingFigure6l*11/10)) + "," + ((legendSpace6/7)+i*legendSpace6/7)+ ")"; })
      			.attr("x", paddingFigure6l/5)
      			.attr("y", 30)
      			.attr("class", "legend")
      			.text(function(d) {  return d.name; })
      			.style("stroke", function(d) { return colorFigure6(d.name); })
      			.style("fill", function(d) { return colorFigure6(d.name); })
      			.style("font-size","20px")
      			.on("click", function(name){
      				console.log(name.name);
      				
                 //Determine if current line is visible 
                 var active   = name.active ? false : true,
                 newOpacity = active ? 0 : 1;
                 //console.log(active);
                 //console.log(newOpacity); 
                 //Hide or show the elements based on the ID
                 d3.select("#tag"+name.name)
                     .transition().duration(100) 
                     .style("opacity", newOpacity);
                if (name.name=="FF"){d3.select("#tagFF")
                     .transition().duration(100) 
                     .style("opacity", newOpacity);} 
                if (name.name=="CHX"){d3.select("#tagCHX")
                     .transition().duration(100) 
                     .style("opacity", newOpacity);} 
                     
                if (name.name=="Data"){d3.select("#tagData")
                     .transition().duration(100) 
                     .style("opacity", newOpacity);}  
                 //Update whether or not the elements are active
                 name.active = active;
                 }) ;
                 
 
lineSvg6.transition()
      		.transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
      		.duration(750);
      		
     		
      	// EXIT
  // Remove old elements as needed.	
  freqFigure6.exit()
  .transition()
    .ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
      .duration(750)
      .style("fill-opacity", 1)
      .remove();
      

       //Update X axis
	svgFigure6.select(".x.axis")
		.transition()
					.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
		.duration(550)
		.attr("transform", "translate(0," + (heightFigure6-3*paddingFigure6) + ")")
		.call(xAxisFigure6)
		.selectAll("text")
    				.attr("y", 0)
    				.attr("x", 10)
    				.attr("dy", ".35em")
    				.attr("transform", "rotate(90)")
    				.style("text-anchor", "start");
		
		
    							
		svgFigure6.select(".xaxis_label")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(550);   
					
    							
	//Update Y axis
	svgFigure6.select(".y.axis")
			.transition()
						.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(550)
			.call(yAxisFigure6);
	

      				        	
 
 };
 
updateFigure6(data);




 				

});
        		