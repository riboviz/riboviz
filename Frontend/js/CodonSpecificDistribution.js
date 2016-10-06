setTimeout(function() {
var marginFigure5 = {top: 20, right: 200, bottom: 70, left: 100},
		paddingFigure5=35,
    	widthFigure5 =800 - marginFigure5.left - marginFigure5.right,
    	heightFigure5 = 400 - marginFigure5.top - marginFigure5.bottom;


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

	var svglineFigure5 = d3.svg.line()
    	.x(function(d) { return xFigure5(d.Position); })
    	.y(function(d) { return yFigure5(d.Reads); });
    	
    	

	var svgFigure5 = d3.select("#CodonBias").append("svg")
    	.attr("width", widthFigure5 + marginFigure5.left + marginFigure5.right)
    	.attr("height", heightFigure5 + marginFigure5.top + marginFigure5.bottom)
  		.append("g")
    		.attr("transform", "translate(" + marginFigure5.left + "," + marginFigure5.top + ")");
    	
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
        			//.text("Ribosome Profile")
        			.style("font-size","13px");
  	   
  		svgFigure5.append("g")
      		.attr("class", "y axis")
      		.call(yAxisFigure5);
      		
      		
      			
d3.selectAll(".form-control")
.on("change.5", change5);

function change5() {
	var year5=d3.select("#yearform").node().value;
	var author5=d3.select("#authorform").node().value;
	var thedataset5=d3.select("#dataform").node().value;
	var string="../Data/";
	
	var thefile=string.concat("F5_Year_", year5, "_Author_", author5,"_Dataset_",thedataset5,"_data.tsv");

	       
	d3.tsv(thefile, function(error, data) {
  		
  		if (error) throw error;
  		data.forEach(function(d) {
  			d.Year = d.Year;
  			d.Author = d.Author;
  			d.Dataset = d.Dataset;
    		d.Position = +d.Position;
    		d.X5prime = +d.X5prime;
    		d.Reads = +d.Reads;
    		d.Type = +d.Type;
    		d.DataType = d.DataType;
  		});
    


function updateFigure5(data, value2, value3) {

var datanewfirst = data;

var datanew = datanewfirst.filter(function(d, key) { 
            return (d.X5prime==value2  & (d.Type==value3[0] | d.Type==value3[1] | d.Type==value3[2])); 
    
    });
   
var value1=d3.map(datanew, function(d){return d.DataType;}).keys();


var nucleotide = d3.nest()
        .key(function(d) {return d.Type;})
        .entries(datanew);
  		
 	
 
 colorFigure5.domain(["Data", "FF", "CHX"]);
 legendSpace =  widthFigure5/nucleotide.length;
 
 
 xFigure5.domain(d3.extent(datanew, function(d) { return d.Position; }));

  		yFigure5.domain([
  			0,
    		d3.max(nucleotide, function(c) { return d3.max(c.values, function(v) { return v.Reads; }); })+
    		d3.max(nucleotide, function(c) { return d3.max(c.values, function(v) { return v.Reads; }); })/5
  		]);
  		
       		 		
var freqFigure5 = svgFigure5.selectAll(".nucleotide")
      		.data(nucleotide);
	

 
  freqFigure5.select("g .nucleotide path")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(750)
      		.attr("class", "line")
      		.attr("id", function(d) { return "tag"+d.key; }) // assign ID
      		.attr("d", function(d) { return svglineFigure5(d.values); })
      		.style("stroke", function(d) { return colorFigure5(d.key); });
	
	 freqFigure5.selectAll("g .nucleotide text")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(500)
			.attr("x", paddingFigure5-20)
      		.attr("y", 135)
      		.attr("class", "legend")
      		.text(function(d) {
      		if(d.name==100 & value1=="RPF"){return " RPF Data"};
      		if(d.name==0 & value1=="RPF"){return "FF"};
      		if(d.name==1 & value1=="RPF"){return "CHX"};
      		if(d.name==1 & value1=="mRNA"){return "BG"};
      		if(d.name==100 & value1=="mRNA"){return " mRNA Data"};
      		 })
      		.style("stroke", function(d) { return colorFigure5(d.name); })
      		.style("fill", function(d) { return colorFigure5(d.name); })
      		.style("font-size","20px");
      		
      		
	  

	var freqgroup=freqFigure5.enter().append("g").attr("class", "nucleotide");
	   
      
      //add path to line group		
  			freqgroup.append("path")
      		.attr("class", "line")
      		.attr("id", function(d) { return "tag"+d.key; }) // assign ID
      		.attr("d", function(d) { return svglineFigure5(d.values); })
      		.style("stroke", function(d) { return colorFigure5(d.key); });
  
 legendSpace = heightFigure5/nucleotide.length;
    
  
     freqgroup.append("text")
     .datum(function(d) { return {name: d.key, value: d.values[d.values.length - 1]}; })
			.attr("transform", function(d, i) { return "translate(" + (widthFigure5) + "," + (i*legendSpace/2) + ")"; })
      		.attr("x", paddingFigure5-20)
      		.attr("y", 135)
      		.attr("class", "legend")
      		.text(function(d) {
      		if(d.name==100 & value1=="RPF"){return " RPF Data"};
      		if(d.name==0 & value1=="RPF"){return "FF"};
      		if(d.name==1 & value1=="RPF"){return "CHX"};
      		if(d.name==1 & value1=="mRNA"){return "BG"};
      		//if(d.name==1& value1=="RPF"){return "CHX"};
      		if(d.name==100 & value1=="mRNA"){return " mRNA Data"};
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



updateFigure5(data, d3.select('input[name="inputsrc51"]:checked').node().value, [0, 1, 100]);

	d3.selectAll(".fig5radio")
      .on("change", changeit5);
    
    function changeit5() {
		var prime= d3.select('input[name="inputsrc51"]:checked').node().value;
		updateFigure5(data,  prime, [0, 1, 100]);
	};
      
d3.select("#download3")
		.on("click", function (){
			window.open(thefile );
		});
		
}); //data

}; //change form

}, 30); //timeout
