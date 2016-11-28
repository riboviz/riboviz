var marginFigure5 = {top: 20, right: 70, bottom: 70, left: 100},
		paddingFigure5=35,
    	widthFigure5 =700 - marginFigure5.left - marginFigure5.right,
    	heightFigure5 = 350 - marginFigure5.top - marginFigure5.bottom;


var xFigure5 = d3.scale.linear()
    	.range([0, widthFigure5]);

	var yFigure5 = d3.scale.linear()
    	.range([heightFigure5, paddingFigure5]).nice();

	var colorFigure5 = d3.scale.category20()
					.domain(["Data","FF", "CHX"])
  					.range(["#a1d99b", "#9ecae1", "#bcbddc"] ); 
	
	var xAxisFigure5 = d3.svg.axis()
    	.scale(xFigure5)
    	.orient("bottom").ticks(15);

	var yAxisFigure5 = d3.svg.axis()
    	.scale(yFigure5)
    	.orient("left").ticks(10);

    	
	var area = d3.svg.area()
    .x(function(d) { return xFigure5(d.Position); })
    .y0(function(d) { return yFigure5(d.Reads+d.SD); })
    .y1(function(d) { return yFigure5(d.Reads-d.SD); });
    
	var svglineFigure5 = d3.svg.line()
    	.x(function(d) { return xFigure5(d.Position); })
    	.y(function(d) { return yFigure5(d.Reads); });
    	

	var svgFigure5 = d3.select("#CodonBias").append("svg")
    	.attr("width", widthFigure5 + marginFigure5.left + marginFigure5.right)
    	.attr("height", heightFigure5 + marginFigure5.top + marginFigure5.bottom)
  		.append("g")
    		.attr("transform", "translate(" + marginFigure5.left + "," + marginFigure5.top + ")");
    		
	
	d3.tsv("../../Data/F5_Year_2016_Author_Weinberg_Dataset_RPF_data.tsv", function(error, data) {
  		
  		if (error) throw error;
  		data.forEach(function(d) {
  			d.Year = d.Year;
  			d.Author = d.Author;
  			d.Dataset = d.Dataset;
    		d.Position = +d.Position;
    		d.X5prime = +d.X5prime;
    		d.Reads = +d.Reads;
    		d.SD = +d.SD;
    		d.Type = +d.Type;
    		d.DataType = d.DataType;
  		});
      	
		//attach x axis
  		svgFigure5.append("g")
      		.attr("class", "x axis")
      		.attr("transform", "translate(0," + heightFigure5 + ")")
      		.call(xAxisFigure5);
      		
      		//x-axis line
    			svgFigure5.selectAll(".x.axis")	
  					.append("rect")
  	   				.attr("width", widthFigure5)
  	   				.attr("height",1)
  	   				.attr("fill","#000")
  
    	svgFigure5.append("text")      // text label for the x axis
      				.attr("class", "xaxis_label")
      				.attr("transform", "translate(" + (widthFigure5 / 2) + " ," + (heightFigure5 + paddingFigure5*8) + ")")
        			.style("text-anchor", "middle")
        			.style("font-size","13px");
  	   
  		svgFigure5.append("g")
      		.attr("class", "y axis")
      		.call(yAxisFigure5);
      		
      		// now add titles to the axes
        svgFigure5.append("text")
            .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
            .attr("transform", "translate("+(0-paddingFigure5)+","+(heightFigure5/2)+")rotate(-90)")  // text is drawn off the screen top left, move down and out and rotate
            .text("Relative normalized reads");

       svgFigure5.append("text")
            .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
            .attr("transform", "translate("+ (widthFigure5/2) +","+(heightFigure5+paddingFigure5)+")")  // centre below axis
            .text("Codon position");
      


function updateFigure5(data, value1, value2, value3) {

var datanewfirst = data;

var datanew = datanewfirst.filter(function(d, key) { 
            return (d.X5prime==value2 & d.DataType==value1 & (d.Type==value3[0] | d.Type==value3[1] | d.Type==value3[2])); 
    
    });
   


var nucleotide = d3.nest()
        .key(function(d) {return d.Type;})
        .entries(datanew);
  		
 	
 
 colorFigure5.domain(["Data", "FF", "CHX", "Other"]);
 legendSpace =  widthFigure5/nucleotide.length;
 
 
 xFigure5.domain(d3.extent(datanew, function(d) { return d.Position; }));

  		yFigure5.domain([
  			0,
    		d3.max(nucleotide, function(c) { return d3.max(c.values, function(v) { return v.Reads+v.SD; }); })+
    		d3.max(nucleotide, function(c) { return d3.max(c.values, function(v) { return v.Reads+v.SD; }); })/5
  		]);
  		
       		 		
var freqFigure5 = svgFigure5.selectAll(".nucleotide")
      		.data(nucleotide);
      		
	
	    freqFigure5.select("g .nucleotide text")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(500);
      		
      		freqFigure5.select("g .nucleotide .area")
      		.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(500)
      		.attr("class", "area")
      		.attr("id", function(d) { return "tag2"+d.key; }) // assign ID
      		.attr("d", function(d) { return area(d.values); })
      		.style("stroke-fill", function(d) { return colorFigure5(d.key); })
      		.style("stroke-opacity", 0.5);
        

 			freqFigure5.select("g .nucleotide .line")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(500)
      		.attr("class", "line")
      		.attr("id", function(d) { return "tag"+d.key; }) // assign ID
      		.attr("d", function(d) { return svglineFigure5(d.values); })
      		.style("stroke", function(d) { return colorFigure5(d.key); });
      	


	var freqgroup=freqFigure5.enter().append("g").attr("class", "nucleotide");
	   
      
          	freqgroup.append("path")
      		.attr("class", "area")
      		.attr("id", function(d) { return "tag2"+d.key; }) // assign ID
      		.attr("d",  function(d) { return area(d.values); })
      		.style("fill", function(d) { return colorFigure5(d.key); })
      		.style("fill-opacity", 0.5);
         
  			freqgroup.append("path")
      		.attr("class", "line")
      		.attr("id", function(d) { return "tag"+d.key; }) // assign ID
      		.attr("d", function(d) { return svglineFigure5(d.values); })
      		.style("stroke", function(d) { return colorFigure5(d.key); });
  
   

 legendSpace = heightFigure5/nucleotide.length;
    
  
     freqgroup.append("text")
     .datum(function(d) { return {name: d.key, value: d.values[d.values.length - 1]}; })
			.attr("transform", function(d, i) { return "translate(" + (widthFigure5) + "," + (i*legendSpace/2) + ")"; })
			//.attr("transform", "translate(" + (widthFigure5) + " ," + (heightFigure5 -paddingFigure5*5) + ")")
      		.attr("x", paddingFigure5-20)
      		.attr("y", 135)
      		.attr("class", "legend")
      		.text(function(d) {
      		if(d.name==100 & value1=="RPF"){return "Data"};
      		if(d.name==0 & value1=="RPF"){return "FF"};
      		if(d.name==1 & value1=="RPF"){return "CHX"};
      		if(d.name==1 & value1=="mRNA"){return "Other"};
      		//if(d.name==1& value1=="RPF"){return "CHX"};
      		if(d.name==100 & value1=="mRNA"){return "Data"};
      		 })
      		.style("stroke", function(d) { return colorFigure5(d.name); })
      		.style("fill", function(d) { return colorFigure5(d.name); })
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
                if (name.name=="0"){d3.select("#tag0")
                     .transition().duration(100) 
                     .style("opacity", newOpacity);} 
                if (name.name=="1"){d3.select("#tag1")
                     .transition().duration(100) 
                     .style("opacity", newOpacity);} 
                     
                if (name.name=="100"){d3.select("#tag100")
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
  freqFigure5.exit()
  .transition()
    .ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
      .duration(750)
      .style("fill-opacity", 1)
      .remove();
      
      
      


       //Update X axis
	svgFigure5.select(".x.axis")
		.transition()
					.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
		.duration(550)
		.call(xAxisFigure5);
					
	//Update Y axis
	svgFigure5.select(".y.axis")
			.transition()
						.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(550)
			.call(yAxisFigure5);
						
		svgFigure5.select(".xaxis_label")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(550);


};



	updateFigure5(data, "RPF",0, [0, 1, 100]);

	d3.selectAll(".fig5radio")
      .on("change", changeit5);
    
    function changeit5() {

		var prime= d3.select('input[name="inputsrc51"]:checked').node().value;

		updateFigure5(data, "RPF", prime, [0, 1, 100]);
	};
      



 	

});
