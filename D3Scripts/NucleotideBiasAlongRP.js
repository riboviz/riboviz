setTimeout(function() {

var marginFigure3 = {top: 0, right: 100, bottom: 70, left: 70},
		paddingFigure3=35,
    	widthFigure3 =750 - marginFigure3.left - marginFigure3.right,
    	heightFigure3 = 400 - marginFigure3.top - marginFigure3.bottom;


	var xFigure3 = d3.scale.linear()
    	.range([0, widthFigure3]);

	var yFigure3 = d3.scale.linear()
    	.range([heightFigure3, paddingFigure3]).nice();

	var colorFigure3 = d3.scale.category20()
					.domain(["A", "T", "C", "G"])
  					.range(["#d6604d", "#bababa", "rgb(239,138,98)" , "#878787"]);
  					
  					
  	var colorFigure32 = d3.scale.category20()
					.domain(["non-CHX", "CHX", "Data"])
  					.range(["red", "black", "grey"]);
  					
  					
	var bisectDateFigure3 = d3.bisector(function(d) { return d.Position; }).right;
	
	var xAxisFigure3 = d3.svg.axis()
    	.scale(xFigure3)
    	.orient("bottom").ticks(15);

	var yAxisFigure3 = d3.svg.axis()
    	.scale(yFigure3)
    	.orient("left").ticks(10);

	var svglineFigure3 = d3.svg.line()
    	.x(function(d) { return xFigure3(d.Position); })
    	.y(function(d) { return yFigure3(d.frequencies); });
    	
    	var svglineFigure33 = d3.svg.line()
    	.x(function(d) { return xFigure3(d.Position); })
    	.y(function(d) { return yFigure3(d.frequencies); });

	var svgFigure3 = d3.select("#Nucleotide").append("svg")
    	.attr("width", widthFigure3 + marginFigure3.left + marginFigure3.right)
    	.attr("height", heightFigure3 + marginFigure3.top + marginFigure3.bottom)
  		.append("g")
    		.attr("transform", "translate(" + marginFigure3.left + "," + marginFigure3.top + ")");
    		
//attach x axis
  		svgFigure3.append("g")
      		.attr("class", "x axis")
      		.attr("transform", "translate(0," + heightFigure3 + ")")
      		.call(xAxisFigure3);
      		
      		//x-axis line
    			svgFigure3.selectAll(".x.axis")	
  					.append("rect")
  	   				.attr("width", widthFigure3)
  	   				.attr("height",1)
  	   				.attr("fill","#000")
  
    	svgFigure3.append("text")      // text label for the x axis
      				.attr("class", "xaxis_label")
      				.attr("transform", "translate(" + (widthFigure3 / 2) + " ," + (heightFigure3 + paddingFigure3*8) + ")")
        			.style("text-anchor", "middle")
        			.text("Ribosome Profile")
        			.style("font-size","13px");
  	   
  		svgFigure3.append("g")
      		.attr("class", "y axis")
      		.call(yAxisFigure3);
      		


d3.selectAll(".form-control")
.on("change.3", change3);


function change3() {
	console.log("testingthisshit");
	var year3=d3.select("#yearform").node().value;
	console.log(year3);
	var author3=d3.select("#authorform").node().value;
	console.log(author3);
	var thedataset3=d3.select("#dataform").node().value;
	console.log(thedataset3);
	


    //var focusFigure3 = svgFigure3.append("g") 
    //				.style("display", "none");
    


	var string="../../VizData/";
	console.log(string);
	var thefile=string.concat("F3_Temp_Year_",year3,"_Author_",author3, "_Dataset_", thedataset3,"_data.tsv");
	console.log(thefile);
	d3.tsv(thefile, function(error, data) {
  		if (error) throw error;
  		data.forEach(function(d) {
    		d.Position = +d.Position;
    		d.Length= +d.Length;
    		d.Frame= +d.Frame;
    		d.CHX= d.CHX;
  		});
    
console.log(data);
console.log("here");
var nrLength=data[data.length-1].Position;

      		
		

function updateFigure3(data, value1, value2, value3) {
console.log("value3");
var datanew = data.filter(function(d) { 
            return (d.Frame==value1 & d.Length==value2 & d.CHX==value3); //& d.CHX==value3
    
    });
 console.log(datanew);   
legendSpace = heightFigure3/datanew.length;


colorFigure3.domain(d3.keys(datanew[0]).filter(function(key) { return (key !== "Length" & key !=="Position" & key !== "Frame" & key !== "SD_A" & key !== "SD_T" & key !== "SD_G" & key !== "SD_C" & key !== "CHX"); }));


  	var nucleotide = colorFigure3.domain().map(function(name) {
    		return {
      			name: name,
      			values: datanew.map(function(d) {
        			return {Position: d.Position, frequencies: +d[name], length: d.Length, Frame: d.Frame, CHX:d.CHX};
      			})
    		};
		});
  		
  		
	xFigure3.domain(d3.extent(datanew, function(d) { return d.Position; }));

  		yFigure3.domain([
  			0,
    		//d3.min(freq, function(c) { return d3.min(c.values, function(v) { return v.frequencies; }); }),
    		d3.max(nucleotide, function(c) { return d3.max(c.values, function(v) { return v.frequencies; }); })+
    		d3.max(nucleotide, function(c) { return d3.max(c.values, function(v) { return v.frequencies; }); })/5
  		]);
  		
 
      		 		
var freqFigure3 = svgFigure3.selectAll(".nucleotide")
      		.data(nucleotide);


  		freqFigure3.select("g .nucleotide path")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(750)
      		.attr("class", "line")
      		.attr("id", function(d) { return "tag"+d.name; }) // assign ID
      		.attr("d", function(d) { return svglineFigure3(d.values); })
      		.style("stroke", function(d) { return colorFigure3(d.name); });
      		
      		freqFigure3.select("g .nucleotide text")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(750)
      		.attr("x", paddingFigure3)
      		.attr("y", 30)
      		.attr("class", "legend")
      		.text(function(d) { return d.name; })
      		.style("stroke", function(d) { return colorFigure3(d.name); })
      		.style("fill", function(d) { return colorFigure3(d.name); })
      		.style("font-size","20px");

    		
      		// ENTER
  // add a line group
      		var freqgroup=freqFigure3.enter().append("g").attr("class", "nucleotide");
      		 console.log(freqgroup);
      //add path to line group		
  			freqgroup.append("path")
      		.attr("class", "line")
      		.attr("id", function(d) { return "tag"+d.name; }) // assign ID
      		.attr("d", function(d) { return svglineFigure3(d.values); })
      		.style("stroke", function(d) { return colorFigure3(d.name); });
      		
      		
            
      	freqgroup.append("text")
  		.datum(function(d) { return {name: d.name, value: d.values[d.values.length - 1]}; })
      		 
      		.attr("transform", function(d, i) { return "translate(" + (xFigure3(d.value.Position)) + "," + ((7*legendSpace)+3*i*legendSpace) + ")"; })
      		.attr("x", paddingFigure3)
      		.attr("y", 30)
      		.attr("class", "legend")
      		.text(function(d) { return d.name; })
      		.style("stroke", function(d) { return colorFigure3(d.name); })
      		.style("fill", function(d) { return colorFigure3(d.name); })
      		.style("font-size","20px")
      		.on("click", function(name){
      				console.log(name.name);
                 //Determine if current line is visible 
                 var active   = name.active ? false : true,
                 newOpacity = active ? 0 : 1;
                 console.log(active);
                 console.log(newOpacity); 
                 //Hide or show the elements based on the ID
                 d3.select("#tag"+name.name)
                     .transition().duration(100) 
                     .style("opacity", newOpacity);
                if (name.name=="A"){d3.select("#tagA")
                     .transition().duration(100) 
                     .style("opacity", newOpacity);} 
                if (name.name=="T"){d3.select("#tagT")
                     .transition().duration(100) 
                     .style("opacity", newOpacity);} 
                     
                if (name.name=="C"){d3.select("#tagC")
                     .transition().duration(100) 
                     .style("opacity", newOpacity);} 
                if (name.name=="G"){d3.select("#tagG")
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
  freqFigure3.exit()
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
	svgFigure3.select(".x.axis")
		.transition()
					.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
		.duration(550)
		.call(xAxisFigure3);
					
	//Update Y axis
	svgFigure3.select(".y.axis")
			.transition()
						.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(550)
			.call(yAxisFigure3);
						
		svgFigure3.select(".xaxis_label")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 500;
			})
			.duration(550);
      		  		

 }; 


function updateFigure3CHX(data, value1, value2, value3) {
      

var datanew2 = data.filter(function(d, key) { 
            return (d.Frame==value1 & d.Length==value2 ); //& d.CHX==value3
    
    });


var thenucleo=value3;
console.log(thenucleo);

//filter before the nesting
console.log(datanew2);

var  datanew3 = datanew2.map( function (d) { 
	//console.log(d.A);
	//console.log(thenucleo);
	if(thenucleo=="A"){thefreqtest=d.A};
	if(thenucleo=="T"){thefreqtest=d.T};
	if(thenucleo=="G"){thefreqtest=d.G};
	if(thenucleo=="C"){thefreqtest=d.C};
    return { 
      frequencies: +thefreqtest,
      Frame: d.Frame,
      Length: d.Length,
      Position: d.Position,
      CHX: d.CHX}; 
});
    
  
//console.log(datanew3); 
var nucleotide = d3.nest()
		
        .key(function(d) {return d.CHX;})
        .entries(datanew3);
        
      //console.log(dataNest);
      
legendSpace2 = heightFigure3/value2;
 console.log(nucleotide);
 
 

//console.log(d3.keys(datanew2));
 
// colorFigure32.domain(d3.keys(dataNest[0]).filter(function(key) { return (key);}));
//  
//    console.log(dataNest.key);     
// console.log(d3.keys(datanew2[0]).filter(function(key) { return (key == "CHX" ); }));
colorFigure32.domain(["non-CHX", "CHX", "Data"]);

       xFigure3.domain(d3.extent(datanew3, function(d) { return d.Position; }));

  		yFigure3.domain([
  			0,
    		//d3.min(freq, function(c) { return d3.min(c.values, function(v) { return v.frequencies; }); }),
    		d3.max(nucleotide, function(c) { return d3.max(c.values, function(v) { return v.frequencies; }); })+
    		d3.max(nucleotide, function(c) { return d3.max(c.values, function(v) { return v.frequencies; }); })/5
  		]);
  		
       		 		
var freqFigure3 = svgFigure3.selectAll(".nucleotide")
      		.data(nucleotide);
			//.enter().append("g")
			//.attr("class", "freq");

  
  freqFigure3.select("g .nucleotide path")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / value2 * 500;
			})
			.duration(750)
      		.attr("class", "line")
      		.attr("id", function(d) { return "tag"+d.key; }) // assign ID
      		.attr("d", function(d) { return svglineFigure33(d.values); })
      		.style("stroke", function(d) { return colorFigure32(d.key); });
    			
	
	    freqFigure3.select("g .nucleotide text")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / value2 * 500;
			})
			.duration(500)
      		.attr("x", paddingFigure3)
      		.attr("y", 30)
      		.attr("class", "legend")
      		.text(function(d) { 
      		if(d.key==0){return "non-CHX"};
      		if(d.key==1){return "CHX"};
      		if(d.key==100){return "Data"};
      		 })
      		.style("stroke", function(d) { return colorFigure32(d.key); })
      		.style("fill", function(d) { return colorFigure32(d.key); })
      		.style("font-size","20px");		
	

//  freqFigure3.selectAll('g .nucleotide text')
//  		.transition()
//   			.ease("linear")
// 			.delay(10)
// 			.remove();
//   	
	
	
	var freqgroup2=freqFigure3.enter().append("g").attr("class", "nucleotide");
    //  console.log(freqgroup2);	
      	

      
      
      //add path to line group		
  			freqgroup2.append("path")
      		.attr("class", "line")
      		.attr("id", function(d) { return "tag"+d.key; }) // assign ID
      		.attr("d", function(d) { return svglineFigure33(d.values); })
      		.style("stroke", function(d) { return colorFigure32(d.key); });
 
        console.log(nucleotide); 
 
     freqgroup2.append("text")
     .data(nucleotide)	 
			.attr("transform", function(d, i) { return "translate(" + (paddingFigure3+xFigure3(d.values.Length)) + "," + ((8*legendSpace2)+5*i*legendSpace2) + ")"; })
			.attr("transform", "translate(" + (widthFigure3/2) + " ," + (heightFigure3 + paddingFigure3*5) + ")")
      		.attr("x", paddingFigure3)
      		.attr("y", 30)
      		.attr("class", "legend")
      		.text(function(d) { 
      		if(d.key==0){return "non-CHX"};
      		if(d.key==1){return "CHX"};
      		if(d.key==100){return "Data"};
      		 })
      		.style("stroke", function(d) { return colorFigure32(["non-CHX", "CHX", "Data"]); })
      		.style("fill", function(d) { return colorFigure32(["non-CHX", "CHX", "Data"]); })
      		.style("font-size","20px")
      		// .on("click", function(key){
//       				console.log(nucleotide.key);
//       				
//                  //Determine if current line is visible 
//                  var active   = name.active ? false : true,
//                  newOpacity = active ? 0 : 1;
//                  console.log(active);
//                  console.log(newOpacity); 
//                  //Hide or show the elements based on the ID
//                  d3.select("#tag"+name.key)
//                      .transition().duration(100) 
//                      .style("opacity", newOpacity);
//                 if (name.key=="0"){d3.select("#tag0")
//                      .transition().duration(100) 
//                      .style("opacity", newOpacity);} 
//                 if (name.key=="1"){d3.select("#tag1")
//                      .transition().duration(100) 
//                      .style("opacity", newOpacity);} 
//                      
//                 if (name.key=="100"){d3.select("#tag100")
//                      .transition().duration(100) 
//                      .style("opacity", newOpacity);}  
//                  //Update whether or not the elements are active
//                  name.active = active;
//                  }) ;
                 
 
 freqgroup2.transition()
      		.transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / value2 * 500;
			})
      		.duration(750);
      		
     		
      	// EXIT
  // Remove old elements as needed.	
  freqFigure3.exit()
  .transition()
    .ease("linear")
			.delay(function(d, i) {
					return i / value2 * 500;
			})
      .duration(750)
      .style("fill-opacity", 1)
      .remove();
      
      
      

       		// now should really update the axes
       //Update X axis
	svgFigure3.select(".x.axis")
		.transition()
					.ease("linear")
			.delay(function(d, i) {
					return i / value2 * 500;
			})
		.duration(550)
		.call(xAxisFigure3);
					
	//Update Y axis
	svgFigure3.select(".y.axis")
			.transition()
						.ease("linear")
			.delay(function(d, i) {
					return i / nucleotide.length * 500;
			})
			.duration(550)
			.call(yAxisFigure3);
						
		svgFigure3.select(".xaxis_label")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / value2 * 500;
			})
			.duration(550);
  //  console.log(datanew3.length);  		  		
             
}; //updateCHX


//	updateFigure3(data, 0, 28, 100);

    
 updateFigure3(data, 0, 28, 100);   
    d3.selectAll(".fig3radio")
      .on("change", changeit3);
	
function changeit3() {

	var nval=d3.select("#nValue").node().value;
	d3.select("#nValue").attr("max",nrLength);
	var frames= d3.select('input[name="inputsrc1"]:checked').node().value;
	var nucleo= d3.select('input[name="inputsrc3"]:checked').node();
	if ( nucleo == null){
		updateFigure3(data, frames, nval, 100);
	}
	 
  			
  	d3.selectAll(".fig3radio2")
      				.on("change.32", change32);
    
    function change32() {
  		var nucleo= d3.select('input[name="inputsrc3"]:checked').node().value;
  		updateFigure3CHX(data, frames, nval, nucleo);
  		
  	};
  	
  d3.select("#Compare")
  				.on("click",  function() {
  				console.log("blah");
  				
  			});
  						
  d3.select("#reset")
  				.on("click",  function() {
  				updateFigure3(data, 0, 28, 100);
  			});
//when the input range changes update value 

};
// 
// function handleClick(event){
//                 console.log(document.getElementById("myVal").value)
//                 return false;
//             };
            

}); //data

}; //change form

}, 10); //timeout
	



