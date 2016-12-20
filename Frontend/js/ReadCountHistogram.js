setTimeout(function() {
            //Width and height
			var marginFigure2 = {top: 0, right: 20, bottom: 70, left: 100},
    			widthFigure2 = 700 - marginFigure2.left - marginFigure2.right,
    			paddingFigure2=30,
    			paddingxFigure2=5,
    			heightFigure2 = 300 - marginFigure2.top - marginFigure2.bottom;


			var xFigure2 = d3.scale.ordinal().rangeRoundBands([0, widthFigure2], .05);
			var yFigure2 = d3.scale.linear().range([heightFigure2-5,paddingFigure2]);

			var xAxisFigure2 = d3.svg.axis()
    			.scale(xFigure2)
    			.orient("bottom");

			var yAxisFigure2 = d3.svg.axis()
    			.scale(yFigure2)
    			.orient("left")
    			.ticks(5);

			var svgFigure2 = d3.select("#ReadCountHistogram").append("svg")
    			.attr("width", widthFigure2 + marginFigure2.left + marginFigure2.right)
    			.attr("height", heightFigure2 + marginFigure2.top + marginFigure2.bottom)
  				.append("g")
    			.attr("transform", "translate(" + marginFigure2.left + "," + marginFigure2.top + ")");

				svgFigure2.append("g")
      				.attr("class", "x axis")
      				.attr("transform", "translate(0," + heightFigure2 + ")")
      				.call(xAxisFigure2)
    				.selectAll("text")
      					.style("text-anchor", "end")
      					.attr("dx", "-.8em")
      					.attr("dy", "-.55em")
      					.attr("transform", "rotate(0)" )
      					.attr("transform", "translate(17," +15+ ")")
      					.style("font-size","16px");
      					
      						svgFigure2.append("g")
      				.attr("class", "y axis")
      				.call(yAxisFigure2);
      				
      				// now add titles to the axes
        svgFigure2.append("text")
            .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
            .attr("transform", "translate("+(0-paddingFigure2*2.2)+","+(heightFigure2/2)+")rotate(-90)")  // text is drawn off the screen top left, move down and out and rotate
            .text("Read count").style("font-size","16px").style("fill","#777777");

       svgFigure2.append("text")
            .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
            .attr("transform", "translate("+ (widthFigure2/2) +","+(heightFigure2+paddingFigure2*1.3)+")")  // centre below axis
            .text("Read length").style("font-size","16px").style("fill","#777777");
			
d3.selectAll(".form-control")
.on("change.2",change2);

function change2() {
	var year2=d3.select("#yearform").node().value;
	var author2=d3.select("#authorform").node().value;
	var thedataset2=d3.select("#dataform").node().value;
	var string="../Data/";
	
	var thefile=string.concat("F2_", year2, "_", author2,"_",thedataset2,".tsv");

queue()
  .defer(d3.tsv,thefile)
  .await(analyze);

function analyze(error, data) {
  if(error) { console.log(error); }	   
  
    			data.forEach(function(d) {
        			d.key = d.Length;
        			d.value = +d.Counts;
    			});
	
						
	
	function updateFigure2(data) {


		var datanew = data;
 


  				xFigure2.domain(datanew.map(function(d) { return d.key; }));
  				yFigure2.domain([0, d3.max(datanew, function(d) { return d.value; })]);
				var key = function(d) {
					return d.key;
				};
  	
   				var freqFigure2 = svgFigure2.selectAll(".datanew2")
      		.data(datanew);
      		
    freqFigure2.select("g .datanew2 rect") 
						.transition()
						 .ease("linear")
 						// .delay(function(d, i) {
//  						   return i / datanew.length * 1000;
//  						})
						.duration(600)
						.style("fill", "rgb(166,206,227)")
      				.attr("x", function(d) { return xFigure2(d.key); })
      				.attr("width", xFigure2.rangeBand())
      				.attr("y", function(d) { return yFigure2(d.value); })
      				.attr("height", function(d) { return heightFigure2 - yFigure2(d.value)-paddingxFigure2; });
					
					
  			
        									
 // ENTER
  // add a line group
var freqgroup2=freqFigure2.enter().append("g").attr("class", "datanew2"); 


freqgroup2.append("rect")
      				.style("fill", "rgb(166,206,227)")
      				.attr("x", function(d) { return xFigure2(d.key); })
      				.attr("width", xFigure2.rangeBand())
      				.attr("y", function(d) { return yFigure2(d.value); })
      				.attr("height", function(d) { return heightFigure2 - yFigure2(d.value)-paddingxFigure2; });
             		
             	
    		 
freqgroup2.transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 1000;
			})
      		.duration(600);
      		
freqFigure2.exit()
  .transition()
    .ease("linear")
			.delay(function(d, i) {
					return i / datanew.length * 1000;
			})
      .duration(600)
      .style("fill-opacity", 1)
      .remove();     		
      		      		      		

			
			

						svgFigure2.select(".x.axis")
							.transition()
							.ease("linear")
							.delay(function(d, i) {
								return i / datanew.length * 500;
							})
							.duration(600)
							.call(xAxisFigure2);
							
							
							//Upposition Y axis
						svgFigure2.select(".y.axis")
							.transition()
							.ease("linear")
							.delay(function(d, i) {
								return i / datanew.length * 500;
							})
							.duration(600)
							.call(yAxisFigure2);
										
      			
};


updateFigure2(data);
d3.select("#download2")
		.on("click", function (){
			window.open(thefile );
		});
		
}; //data

}; //change form

}, 500); //timeout