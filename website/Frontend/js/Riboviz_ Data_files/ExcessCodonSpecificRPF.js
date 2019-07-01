setTimeout(function() {

//Width and height
			var marginFigure5 = {top: 20, right: 10, bottom: 50, left: 100},
    			widthFigure5 = 700 - marginFigure5.left - marginFigure5.right,
    			paddingFigure5=5,
    			paddingFigure5l=50
    			heightFigure5 = 300 - marginFigure5.top - marginFigure5.bottom;
    			

					var colorFigure5 = d3.scale.category20()
					.domain(["Data", "FF", "CHX", ])
  					.range(["rgb(239,138,98)", "#a1d99b", "#9ecae1"]);
    		
    		//Create scale functions	
    		var xFigure5 = d3.scale.ordinal()
    					.rangeRoundBands([paddingFigure5, widthFigure5 - paddingFigure5 * 2], .1);
    						

			var yFigure5 = d3.scale.linear()
    						.range([heightFigure5-paddingFigure5*2, paddingFigure5]);
    			
			
			//Define X axis
			var xAxisFigure5 = d3.svg.axis()
    			.scale(xFigure5)
    			.orient("bottom");
			
			//Define Y axis
			var yAxisFigure5 = d3.svg.axis()
    			.scale(yFigure5)
    			.orient("left");
    		
			var areaF5 = d3.svg.area()
    			.x(function(d) { return xFigure5(d.Codon); })
    			.y0(function(d) { return yFigure5(d.Mean+d.SD); })
    			.y1(function(d) { return yFigure5(d.Mean-d.SD); });
    
			//Define line
			var valuelineFigure5 = d3.svg.line()
    			.x(function(d) { return xFigure5(d.Codon); })
    			.y(function(d) { return yFigure5(d.Mean); });
    			
    		
			
			//Create SVG element
			var svgFigure5 = d3.select("#ExcessCodonSpecificRPF")
						.append("svg")
    					.attr("width", widthFigure5 + 2*marginFigure5.left + marginFigure5.right)
    					.attr("height", heightFigure5 + marginFigure5.top + marginFigure5.bottom)
  						.append("g")
    					.attr("transform", "translate(" + marginFigure5.left + "," + marginFigure5.top + ")");

    
//x-axis	
  				svgFigure5.append("g")
      				.attr("class", "x axis")
      				.attr("transform", "translate(0," + (heightFigure5-3*paddingFigure5) + ")")
      			.call(xAxisFigure5)
      			.selectAll("text")
    				.attr("y", 0)
    				.attr("x", 10)
    				.attr("dy", ".35em")
    				.attr("transform", "rotate(90)")
    				.style("text-anchor", "start");
  
    			//x-axis line
    			svgFigure5.selectAll(".x.axis")	
  					.append("rect")
  	   				.attr("width", widthFigure5)
  	   				.attr("height",1)
  	   				.attr("fill","#000")
				
				//y-axis
  				svgFigure5.append("g")
      				.attr("class", "y axis")
      				.attr("transform", "translate(" + 0 + "," + (0-paddingFigure5) + ")")
      				.call(yAxisFigure5);
      			
      			// now add titles to the axes
        		svgFigure5.append("text")
            		.attr("text-anchor", "middle")  
            		.attr("transform", "translate("+(0-paddingFigure5*6.5)+","+(heightFigure5/2)+")rotate(-90)")  
            		.text("Log2 ratio of normalized reads").style("font-size","16px");

       			svgFigure5.append("text")
            		.attr("text-anchor", "middle")  
            		.attr("transform", "translate("+ (widthFigure5/2) +","+(heightFigure5+paddingFigure5*6.8)+")")  
            		.text("Codon").style("font-size","16px");	
      					
d3.selectAll(".form-control")
.on("change.5", change5);


function change5() {
	var year5=d3.select("#yearform").node().value;
	var author5=d3.select("#authorform").node().value;
	var thedataset5=d3.select("#dataform").node().value;
	
	var string="../Data/";
	
	var thefile=string.concat("F5_", year5, "_", author5,"_",thedataset5,".tsv");

	       
	d3.tsv(thefile, function(error, data) {
  				if (error) throw error;
  				data.forEach(function(d) {
    				d.Codon = d.Codon;
    				d.Mean = +d.Mean;
    				d.BG=d.BG;
    				d.SD=+d.SD;
  				});

				
      				
function updateFigure5(data) {
      				
var datanew = data; 

    			var nucleotide = d3.nest()
        			.key(function(d) {return d.BG;})
        			.entries(datanew);
        
      
				legendSpace5 =  widthFigure5/nucleotide.length;
  				colorFigure5.domain(["Data", "FF", "CHX"]);
  

  				xFigure5.domain(datanew.map(function(d) { return d.Codon; }));
  				yFigure5.domain([d3.min(datanew, function(d) { return d.Mean-d.SD; }),d3.max(datanew, function(d) { return d.Mean+d.SD; })]);


var freqFigure5 = svgFigure5.selectAll(".nucleotide")
      		.data(nucleotide);

freqFigure5.select("g .nucleotide .area")
      		.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(500)
      		.attr("class", "area")
      		.attr("id", function(d) { return "tag2"+d.key; }) // assign ID
      		.attr("d", function(d) { return areaF5(d.values); })
      		.style("stroke-fill", function(d) { return colorFigure5(d.key); })
      		.style("stroke-opacity", 0.15);
        

 freqFigure5.select("g .nucleotide .line")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(750)
      		.attr("class", "line")
             		.attr("id", function(d) { return "tag"+d.key; }) // assign ID
             		.attr("d", function(d) { return valuelineFigure5(d.values); })
      		        .style("stroke", function(d) { return colorFigure5(d.key); })
      		        .filter(function(d) { return !isNaN(d.values)});;
      		
   freqFigure5.select("g .nucleotide text")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			});
      		   		
      		      		
var lineSvg5=freqFigure5.enter().append("g").attr("class", "nucleotide");
        
        		lineSvg5.append("path")
      				.attr("class", "area")
      				.attr("id", function(d) { return "tag2"+d.key; }) // assign ID
      				.attr("d",  function(d) { return areaF5(d.values); })
      				.style("fill", function(d) { return colorFigure5(d.key); })
      				.style("fill-opacity", 0.15);
				//Create the line
    			lineSvg5.append("path")
             		.attr("class", "line")
             		.attr("id", function(d) { return "tag"+d.key; }) // assign ID
             		.attr("d", function(d) { return valuelineFigure5(d.values); })
      		        .style("stroke", function(d) { return colorFigure5(d.key); });

 

      lineSvg5.append("text")
    			.datum(function(d) { return {name: d.key, value: d.values[d.values.length - 1]}; })
				.attr("transform", function(d, i) { return "translate(" + ((widthFigure5-paddingFigure5l*11/10)) + "," + ((legendSpace5/7)+i*legendSpace5/7)+ ")"; })
      			.attr("x", paddingFigure5l/5)
      			.attr("y", 10)
      			.attr("class", "legend")
      			.text(function(d) {  return d.name; })
      			.style("stroke", function(d) { return colorFigure5(d.name); })
      			.style("fill", function(d) { return colorFigure5(d.name); })
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
                if (name.name=="non-CHX"){d3.select("#tagFF")
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
                 
 
lineSvg5.transition()
      		.transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
      		.duration(750);
      		
     		
      	// EXIT
  // Remove old elements as needed.	
  freqFigure5.exit()
  .transition()
    .ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
      .duration(750)
      .style("fill-opacity", 1)
      .remove();
      
    
  
       		// now should really update the axes
       //Update X axis
	svgFigure5.select(".x.axis")
		.transition()
					.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
		.duration(550)
		.attr("transform", "translate(0," + (heightFigure5-3*paddingFigure5) + ")")
		.call(xAxisFigure5)
		.selectAll("text")
    				.attr("y", 0)
    				.attr("x", 10)
    				.attr("dy", ".35em")
    				.attr("transform", "rotate(90)")
    				.style("text-anchor", "start");
		
		
    							
		svgFigure5.select(".xaxis_label")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(550);   
					
    							
	//Update Y axis
	svgFigure5.select(".y.axis")
			.transition()
						.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(550)
			.call(yAxisFigure5);
	

      				        	
 
 };
 


updateFigure5(data);

d3.select("#download4")
		.on("click", function (){
			window.open(thefile );
		});
		
}); //data

}; //change form

}, 500); //timeout       		