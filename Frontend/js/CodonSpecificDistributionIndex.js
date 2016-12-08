var marginFigure4 = {top: 20, right: 70, bottom: 70, left: 100},
		paddingFigure4=35,
    	widthFigure4 =700 - marginFigure4.left - marginFigure4.right,
    	heightFigure4 = 350 - marginFigure4.top - marginFigure4.bottom;


var xFigure4 = d3.scale.linear()
    	.range([0, widthFigure4]);

	var yFigure4 = d3.scale.linear()
    	.range([heightFigure4, paddingFigure4]).nice();

	var colorFigure4 = d3.scale.category20()
					.domain(["Data","FF", "CHX"])
  					.range(["rgb(239,138,98)", "#a1d99b", "#9ecae1"] ); 
	
	var xAxisFigure4 = d3.svg.axis()
    	.scale(xFigure4)
    	.orient("bottom").ticks(15);

	var yAxisFigure4 = d3.svg.axis()
    	.scale(yFigure4)
    	.orient("left").ticks(10);

    	
	var area = d3.svg.area()
    .x(function(d) { return xFigure4(d.Position); })
    .y0(function(d) { return yFigure4(d.Mean+d.SD); })
    .y1(function(d) { return yFigure4(d.Mean-d.SD); });
    
	var svglineFigure4 = d3.svg.line()
    	.x(function(d) { return xFigure4(d.Position); })
    	.y(function(d) { return yFigure4(d.Mean); });
    	

	var svgFigure4 = d3.select("#CodonBias").append("svg")
    	.attr("width", widthFigure4 + marginFigure4.left + marginFigure4.right)
    	.attr("height", heightFigure4 + marginFigure4.top + marginFigure4.bottom)
  		.append("g")
    		.attr("transform", "translate(" + marginFigure4.left + "," + marginFigure4.top + ")");
    		
	
	d3.tsv("../../Data/F4_2016_Weinberg_RPF.tsv", function(error, data) {
  		
  		if (error) throw error;
  		data.forEach(function(d) {
    		d.Position = +d.Position;
    		d.End = +d.End;
    		d.Mean = +d.Mean;
    		d.SD = +d.SD;
    		d.BG = d.BG;
    		d.DataType = d.DataType;
  		});

      	
		//attach x axis
  		svgFigure4.append("g")
      		.attr("class", "x axis")
      		.attr("transform", "translate(0," + heightFigure4 + ")")
      		.call(xAxisFigure4);
      		
      		//x-axis line
    			svgFigure4.selectAll(".x.axis")	
  					.append("rect")
  	   				.attr("width", widthFigure4)
  	   				.attr("height",1)
  	   				.attr("fill","#000")
  
    	svgFigure4.append("text")      // text label for the x axis
      				.attr("class", "xaxis_label")
      				.attr("transform", "translate(" + (widthFigure4 / 2) + " ," + (heightFigure4 + paddingFigure4*8) + ")")
        			.style("text-anchor", "middle")
        			.style("font-size","13px");
  	   
  		svgFigure4.append("g")
      		.attr("class", "y axis")
      		.call(yAxisFigure4);
      		
      		// now add titles to the axes
        svgFigure4.append("text")
            .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
            .attr("transform", "translate("+(0-paddingFigure4)+","+(heightFigure4/2)+")rotate(-90)")  // text is drawn off the screen top left, move down and out and rotate
            .text("Relative normalized reads");

       svgFigure4.append("text")
            .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
            .attr("transform", "translate("+ (widthFigure4/2) +","+(heightFigure4+paddingFigure4)+")")  // centre below axis
            .text("Codon position");
      


function updateFigure4(data, value1, value2, value3) {

var datanewfirst = data;

var datanew = datanewfirst.filter(function(d, key) { 
            return (d.End==value2 & d.DataType==value1 & (d.BG==value3[0] | d.BG==value3[1] | d.BG==value3[2])); 
    
    });
   


var nucleotide = d3.nest()
        .key(function(d) {return d.BG;})
        .entries(datanew);
  		
 console.log(nucleotide);
 
 colorFigure4.domain(["Data", "FF", "CHX"]);
 legendSpace =  widthFigure4/nucleotide.length;
 
 
  xFigure4.domain(d3.extent(datanew, function(d) { return d.Position; }));
// 
//   		yFigure4.domain([
//   			0,
//     		d3.max(nucleotide, function(c) { return d3.max(c.values, function(v) { return v.Mean+v.SD; }); })+
//     		d3.max(nucleotide, function(c) { return d3.max(c.values, function(v) { return v.Mean+v.SD; }); })/5
//   		]);
  		
 yFigure4.domain([0, 6]);

   		 		
var freqFigure4 = svgFigure4.selectAll(".nucleotide")
      		.data(nucleotide);
      		
	
	    freqFigure4.select("g .nucleotide text")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(500);
      		
      		freqFigure4.select("g .nucleotide .area")
      		.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(500)
      		.attr("class", "area")
      		.attr("id", function(d) { return "tag2"+d.key; }) // assign ID
      		.attr("d", function(d) { return area(d.values); })
      		.style("stroke-fill", function(d) { return colorFigure4(d.key); })
      		.style("stroke-opacity", 0.5);
        

 			freqFigure4.select("g .nucleotide .line")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(500)
      		.attr("class", "line")
      		.attr("id", function(d) { return "tag"+d.key; }) // assign ID
      		.attr("d", function(d) { return svglineFigure4(d.values); })
      		.style("stroke", function(d) { return colorFigure4(d.key); });
      	


	var freqgroup=freqFigure4.enter().append("g").attr("class", "nucleotide");
	   
      
          	freqgroup.append("path")
      		.attr("class", "area")
      		.attr("id", function(d) { return "tag2"+d.key; }) // assign ID
      		.attr("d",  function(d) { return area(d.values); })
      		.style("fill", function(d) { return colorFigure4(d.key); })
      		.style("fill-opacity", 0.5);
         
  			freqgroup.append("path")
      		.attr("class", "line")
      		.attr("id", function(d) { return "tag"+d.key; }) // assign ID
      		.attr("d", function(d) { console.log(d.values);return svglineFigure4(d.values); })
      		.style("stroke", function(d) { return colorFigure4(d.key); });
  
   

 legendSpace = heightFigure4/nucleotide.length;
    
  
     freqgroup.append("text")
     .datum(function(d) { return {name: d.key, value: d.values[d.values.length - 1]}; })
			.attr("transform", function(d, i) { return "translate(" + (widthFigure4) + "," + (i*legendSpace/2) + ")"; })
			//.attr("transform", "translate(" + (widthFigure4) + " ," + (heightFigure4 -paddingFigure4*5) + ")")
      		.attr("x", paddingFigure4-20)
      		.attr("y", 135)
      		.attr("class", "legend")
      		.text(function(d) {
      		if(d.name=="Data" & value1=="RPF"){return "Data"};
      		if(d.name=="FF" & value1=="RPF"){return "FF"};
      		if(d.name=="CHX" & value1=="RPF"){return "CHX"};
      		if(d.name==1 & value1=="mRNA"){return "Other"};
      		//if(d.name==1& value1=="RPF"){return "CHX"};
      		if(d.name==100 & value1=="mRNA"){return "Data"};
      		 })
      		.style("stroke", function(d) { return colorFigure4(d.name); })
      		.style("fill", function(d) { return colorFigure4(d.name); })
      		.style("font-size","20px")
      		.on("click", function(name){
      				console.log(name.name);
                 
                 //Determine if current line is visible 
                 var active   = name.active ? false : true,
                 newOpacity = active ? 0 : 1;
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
                 
 
 freqgroup.transition()
      		.transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
      		.duration(750);
      		
     		
      	// EXIT
  // Remove old elements as needed.	
  freqFigure4.exit()
  .transition()
    .ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
      .duration(750)
      .style("fill-opacity", 1)
      .remove();
      
      
      


       //Update X axis
	svgFigure4.select(".x.axis")
		.transition()
					.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
		.duration(550)
		.call(xAxisFigure4);
					
	//Update Y axis
	svgFigure4.select(".y.axis")
			.transition()
						.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(550)
			.call(yAxisFigure4);
						
		svgFigure4.select(".xaxis_label")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(550);


};



	updateFigure4(data, "RPF",3, ["FF", "CHX", "Data"]);

	d3.selectAll(".fig4radio")
      .on("change", changeit4);
    
    function changeit4() {

		var prime= d3.select('input[name="inputsrc51"]:checked').node().value;

		updateFigure4(data, "RPF", prime, ["FF", "CHX", "Data"]);
	};
      



 	

});
