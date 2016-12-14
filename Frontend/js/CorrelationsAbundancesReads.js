setTimeout(function() {

var MyApp = {
    		thecorr: null
		};


var marginFigure6 = {top: 100, right: 200, bottom: 50, left: 100},
		paddingFigure6=35,
    	widthFigure6 = 800 - marginFigure6.left - marginFigure6.right,
    	heightFigure6 = 500 - marginFigure6.top - marginFigure6.bottom;
    	
var marginFigure6C = {top: 20, right: 60, bottom: 70, left:100},
		paddingFigure6C=35,
		paddingxFigure6C=10,
    	widthFigure6C = 350 - marginFigure6C.left - marginFigure6C.right,
    	heightFigure6C = 350 - marginFigure6C.top - marginFigure6C.bottom;
    
var xFigure6 = d3.scale.linear()
    .range([paddingFigure6, widthFigure6]);

var yFigure6 = d3.scale.linear()
    .range([heightFigure6-paddingFigure6, 0]);

var colorFigure6 = d3.scale.category20()
  			.range(["#9ecae1"]);

var xAxisFigure6 = d3.svg.axis()
    .scale(xFigure6)
    .orient("bottom");
    //.ticks(7);

var yAxisFigure6 = d3.svg.axis()
    .scale(yFigure6)
    .orient("left");//.ticks(7);


var svgFigure6 = d3.select("#CorrelationsAbundancesReads").append("svg")
    .attr("width", widthFigure6 + marginFigure6.left + marginFigure6.right)
    .attr("height", heightFigure6+2*marginFigure6.top + marginFigure6.bottom)
  .append("g")
    .attr("transform", "translate(" + marginFigure6.left + "," + marginFigure6.top + ")");
 
			
var xFigure6C = d3.scale.ordinal().rangeRoundBands([0, widthFigure6C], .05);
var yFigure6C = d3.scale.linear().range([heightFigure6C-paddingFigure6C,paddingFigure6C]);

var xAxisFigure6C = d3.svg.axis()
    			.scale(xFigure6C)
    			.orient("bottom");

var yAxisFigure6C = d3.svg.axis()
    			.scale(yFigure6C)
    			.orient("left")
    			.ticks(10);

var svgFigure6C = d3.select("#HistogramsCorrsnonCHX").append("svg")
    			.attr("width", widthFigure6C + marginFigure6C.left + marginFigure6C.right)
    			.attr("height", heightFigure6C + marginFigure6C.top + marginFigure6C.bottom)
  				.append("g")
    			.attr("transform", "translate(" + marginFigure6C.left + "," + marginFigure6C.top + ")");

var xFigure6nC = d3.scale.ordinal().rangeRoundBands([0, widthFigure6C], .05);
var yFigure6nC = d3.scale.linear().range([heightFigure6C-paddingFigure6C,paddingFigure6C]);

var xAxisFigure6nC = d3.svg.axis()
    			.scale(xFigure6nC)
    			.orient("bottom");

var yAxisFigure6nC = d3.svg.axis()
    			.scale(yFigure6nC)
    			.orient("left")
    			.ticks(10);

var svgFigure6nC = d3.select("#HistogramsCorrsCHX").append("svg")
    			.attr("width", widthFigure6C + marginFigure6C.left + marginFigure6C.right)
    			.attr("height", heightFigure6C + marginFigure6C.top + marginFigure6C.bottom)
  				.append("g")
    			.attr("transform", "translate(" + marginFigure6C.left + "," + marginFigure6C.top + ")");


					svgFigure6.append("g")
      				.attr("class", "x axis")
      				.attr("transform", "translate(0," + heightFigure6 + ")")
      				.call(xAxisFigure6)
      				.selectAll("text")
      					.style("text-anchor", "end")
      					.attr("dx", "-.8em")
      					.attr("dy", "-.55em")
      					.attr("transform", "rotate(0)" )
      					.attr("transform", "translate(17," +15+ ")")
      					.style("font-size","12px");
      				
      				
      				
      					
    				    				
				// svgFigure6.append("text")      // text label for the x axis
// 					.attr("class", "xaxis_label1")
// 					.attr("transform", "translate(" + (widthFigure6 / 2) + " ," + (heightFigure6 + paddingFigure6*8) + ")")
// 					.style("text-anchor", "middle")
// 					//.text("Length")
// 					.style("font-size","13px");
        				
				//y-axis
  				svgFigure6.append("g")
      				.attr("class", "y axis")
      				//.attr("transform", "translate("  paddingFigure6 + ",0)")
      				.call(yAxisFigure6);
	
      			svgFigure6.append("text")      // text label for the y axis
      				.attr("class", "yaxis_label")
        			.attr("y", heightFigure6 /2 )
        			.attr("x", -paddingFigure6*6)
        			.style("text-anchor", "middle")
        			.attr("transform", "rotate(0)")
        			//.text("FEatg")
        			.style("font-size","13px");
        		
        		        			
        		svgFigure6.append("text")
					.attr("class", "r-label")
					.attr("y", -1*paddingFigure6 )
        			.attr("x", 2.6*paddingFigure6)
        			.style("text-anchor", "middle")
        			.attr("transform", "rotate(0)")
        			//.text("FEatg")
        			.style("font-size","13px");
        			
        			svgFigure6.append("line")
					.attr("class", "r-line")
 					.attr("stroke", "black")
 					.style("stroke-width",0.9);
 						
 						
 						
 				
        svgFigure6.append("text")
            .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
            .attr("transform", "translate("+(0-paddingFigure6)+","+(heightFigure6/2)+")rotate(-90)")  // text is drawn off the screen top left, move down and out and rotate
            .text("");

       svgFigure6.append("text")
            .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
            .attr("transform", "translate("+ (widthFigure6/2) +","+(heightFigure6+paddingFigure6)+")")  // centre below axis
            .text("A, P or E").style("font-size","16px").style("fill","#777777");
        
        svgFigure6C.append("text")
            .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
            .attr("transform", "translate("+ (widthFigure6C/2) +","+(heightFigure6C+paddingFigure6C)+")")  // centre below axis
            .text("background correlations CHX pre-treatment").style("fill","#777777");
            
            svgFigure6nC.append("text")
            .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
            .attr("transform", "translate("+ (widthFigure6C/2) +","+(heightFigure6C+paddingFigure6C)+")")  // centre below axis
            .text("background correlations FF pre-treatment").style("fill","#777777");
            
            
svgFigure6C.append("g").append("line")
		.attr("class", "r-line").attr("stroke", "grey")
 					.style("stroke-width",0.9);		
    					
svgFigure6C.append("g")
      				.attr("class", "x axis")
      				.attr("transform", "translate(0," + heightFigure6C + ")")
      				.call(xAxisFigure6C);
    
svgFigure6C.append("text")
    .attr("class", "xaxis_label")
      					.style("text-anchor", "end")
      					.attr("dx", "-.8em")
      					.attr("dy", "-.55em")
      					.attr("transform", "rotate(0)" )
      					.attr("transform", "translate(17," +15+ ")")
      					.style("font-size","12px");
      					
 
 svgFigure6nC.append("g").append("line")
		.attr("class", "r-line2").attr("stroke", "grey")
 					.style("stroke-width",0.9);		
    					
svgFigure6nC.append("g")
      				.attr("class", "x axis")
      				.attr("transform", "translate(0," + heightFigure6C + ")")
      				.call(xAxisFigure6nC);
    
svgFigure6nC.append("text")
    .attr("class", "xaxislabel")
      					.style("text-anchor", "end")
      					.attr("dx", "-.8em")
      					.attr("dy", "-.55em")
      					.attr("transform", "rotate(0)" )
      					.attr("transform", "translate(17," +15+ ")")
      					.style("font-size","12px");
      					
      					        			
d3.selectAll(".form-control")
.on("change.6", change6);


function change6() {
	var year6=d3.select("#yearform").node().value;
	var author6=d3.select("#authorform").node().value;
	var thedataset6=d3.select("#dataform").node().value;
	var string6="../Data/";
	var thefile6=string6.concat("F6_",year6,"_",author6, "_", thedataset6,".tsv");	
	
queue()
  .defer(d3.tsv,thefile6)
  .defer(d3.tsv, "../Data/F6_Counts.tsv")
  .await(analyze);

function analyze(error, data, data1) {
  if(error) { console.log(error); }

  
  	data.forEach(function(d) {
    d.A = +d.A;
    d.P = +d.P;
    d.E = +d.E;
    d.GCN = +d.tRNA;
    d.tAI = +d.tAI ;
    d.RNAseq = +Math.log10(d.RNAseq);
    d.Microarray = +Math.log10(d.Microarray);
  });

function updateFigure6(data, value1, value2) {


var  data = data.map( function (d) { 
	if(value2=="A"){choosethetwo=d.A};
	if(value2=="E"){choosethetwo=d.E};
	if(value2=="P"){choosethetwo=d.P};
	if(value1=="GCN"){choosetheone=d.GCN};
	if(value1=="tAI"){choosetheone=d.tAI}; //.toFixed(9)
	if(value1=="RNAseq"){choosetheone=d.RNAseq};
	if(value1=="Microarray"){choosetheone=d.Microarray};
	
    return { 
      d1: +choosetheone,
      d2: +choosethetwo}; 
});
    
    
    
	xFigure6.domain(d3.extent(data, function(d) { return d.d1; }));
  	yFigure6.domain(d3.extent(data, function(d) { return d.d2; }));
  
  // DATA JOIN
  // Join new data with old elements, if any.
  var textFigure6 = svgFigure6.selectAll("circle")
      .data(data);

  // UPDATE
  // Update old elements as needed.
  textFigure6.attr("class", "update")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      .attr("cx", function(d) { return xFigure6(d.d1); })
      .attr("cy", function(d) { return yFigure6(d.d2); })
      .filter(function(d) { return !isNaN(d.d1) && !isNaN(d.d2)});

  // ENTER
  // Create new elements as needed.
  textFigure6.enter().append("circle")
      .attr("class", "enter")
      .attr("cx", function(d) { return xFigure6(d.d1); })
      .attr("cy", function(d) { return yFigure6(d.d2); })
      .style("fill-opacity", 0.2)
      .style("fill","#bcbddc")
      .attr("r", 5)
      .filter(function(d) { return !isNaN(d.d1) && !isNaN(d.d2)})
      //.text(function(d) { return d; })
    .transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      .attr("cx", function(d) { return xFigure6(d.d1); })
      .attr("cy", function(d) { return yFigure6(d.d2); })
      .style("fill-opacity", 1)
      .filter(function(d) { return !isNaN(d.d1) && !isNaN(d.d2)});


textFigure6.exit().attr("class", "exit")
    .transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      .style("fill-opacity", 1e-6)
      .remove();
      
	  
		var xSeries = data.map(function(d) { return d.d1; });
		xSeriesreg=xSeries.sort(function(a, b) {
  				return Number(a) - Number(b);
				});


	var xSeriesreal = data.map(function(d) { return d.d1; });
	var ySeries = data.map(function(d) { return d.d2; });
	thecorv=[xSeriesreal, ySeries];
	var thecor=pearsonCorrelation(thecorv,0,1);
	var leastSquaresCoeff = leastSquares(xSeriesreal, ySeries);
	var x1 = xSeriesreg[0];
	var y1 = leastSquaresCoeff[0] *x1+ leastSquaresCoeff[1];
	var x2 = xSeriesreg[xSeriesreg.length - 1];
	var y2 =leastSquaresCoeff[0] *x2+ leastSquaresCoeff[1]; //leastSquaresCoeff[0] * xSeries[xSeries.length-1] + leastSquaresCoeff[1];
	var trendData = [[x1,y1,x2,y2]];
    MyApp.thecorr=thecor;
      
    svgFigure6.select(".r-label")			
			.text("Correlation coefficient: r=" + thecor.toFixed(2))//leastSquaresCoeff[2])
			.style("font-size","16px").style("fill","#777777")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(750);

			
	svgFigure6.select(".r-line")			
			.attr("x1", xFigure6(trendData[0][0]))
  			.attr("y1", yFigure6(trendData[0][1]))
  			.attr("x2",  xFigure6(trendData[0][2]))
  			.attr("y2", yFigure6(trendData[0][3]))
  			.attr("stroke", "grey")
 			.style("stroke-width",1.1)
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(750);
			
			
      // now should really update the axes
      //Update X axis
	svgFigure6.select(".x.axis")
		.transition()
					.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
		.duration(250)
		.call(xAxisFigure6);
					
	//Update Y axis
	svgFigure6.select(".y.axis")
			.transition()
						.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(250)
			.call(yAxisFigure6);
						
		svgFigure6.select(".xaxis_label1")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(250);
            // .text("d1");
						
		svgFigure6.select(".yaxis_label")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(250);
             //.text("d2");
};

function pearsonCorrelation(prefs, p1, p2) {
  var si = [];

  for (var key in prefs[p1]) {
    if (prefs[p2][key]) si.push(key);
  }

  var n = si.length;

  if (n == 0) return 0;

  var sum1 = 0;
  for (var i = 0; i < si.length; i++) sum1 += prefs[p1][si[i]];

  var sum2 = 0;
  for (var i = 0; i < si.length; i++) sum2 += prefs[p2][si[i]];

  var sum1Sq = 0;
  for (var i = 0; i < si.length; i++) {
    sum1Sq += Math.pow(prefs[p1][si[i]], 2);
  }

  var sum2Sq = 0;
  for (var i = 0; i < si.length; i++) {
    sum2Sq += Math.pow(prefs[p2][si[i]], 2);
  }

  var pSum = 0;
  for (var i = 0; i < si.length; i++) {
    pSum += prefs[p1][si[i]] * prefs[p2][si[i]];
  }

  var num = pSum - (sum1 * sum2 / n);
  var den = Math.sqrt((sum1Sq - Math.pow(sum1, 2) / n) *
      (sum2Sq - Math.pow(sum2, 2) / n));

  if (den == 0) return 0;

  return num / den;
};



// returns slope, intercept and r-square of the line
	function leastSquares(xSeries, ySeries) {
		var reduceSumFunc = function(prev, cur) { return prev + cur; };
		
		var xBar = xSeries.reduce(reduceSumFunc) * 1.0 / xSeries.length;
		var yBar = ySeries.reduce(reduceSumFunc) * 1.0 / ySeries.length;

		var ssXX = xSeries.map(function(d) { return Math.pow(d - xBar, 2); })
			.reduce(reduceSumFunc);
		
		var ssYY = ySeries.map(function(d) { return Math.pow(d - yBar, 2); })
			.reduce(reduceSumFunc);
			
		var ssXY = xSeries.map(function(d, i) { return (d - xBar) * (ySeries[i] - yBar); })
			.reduce(reduceSumFunc);
			
		var slope = ssXY / ssXX;
		var intercept = yBar - (xBar * slope);
		var rSquare = Math.pow(ssXY, 2) / (ssXX * ssYY);
		
		return [slope, intercept, rSquare];
	};





//now move to the next panel   			


data1.forEach(function(d) {
        			d.Site = d.Site;
        			d.Abundance = d.Abundance;
        			d.Counts=+d.Counts;
        			d.Values=+d.Values;
        			d.CHX=d.CHX
    			});  			
      					
      					    					
//update function for the second panel     			
function updatehistogram1(data, value1, value2) {

			
var data = data.filter(function(d) { 
            return (d.Site==value2 & d.Abundance==value1 & d.CHX=="non-CHX"); //& d.CHX==value3
    
    });

xFigure6C.domain(data.map(function(d) { return d.Values; }));
yFigure6C.domain([0, d3.max(data, function(d) { return d.Counts; })]);

var values = function(d) {
					return d.Values;
				};

function closest (num, arr) {
                var curr = arr[0];
                var diff = Math.abs (num - curr);
                for (var val = 0; val < arr.length; val++) {
                    var newdiff = Math.abs (num - arr[val]);
                    if (newdiff < diff) {
                        diff = newdiff;
                        curr = arr[val];
                    }
                }
                return curr;
            }
            
var testnumber=MyApp.thecorr;

var testarray=data.map(function(d) { return d.Values; });
var closestnr=closest(testnumber, testarray);
    					
      					
var textFigure6 = svgFigure6C.selectAll("rect")
      .data(data);
      
      	 textFigure6.attr("class", "update")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      				.attr("x", function(d) { return xFigure6C(d.Values); })
      				.attr("width", xFigure6C.rangeBand())
      				.attr("y", function(d) { return yFigure6C(d.Counts); })
      				.attr("height", function(d) { return heightFigure6C - yFigure6C(d.Counts)-paddingxFigure6C; }).style("fill", "#bcbddc");

  textFigure6.enter().append("rect")
      .attr("class", "enter")
      .attr("x", function(d) { return xFigure6C(d.Values); })
      				.attr("width", xFigure6C.rangeBand())
      				.attr("y", function(d) { return yFigure6C(d.Counts); })
      				.attr("height", function(d) { return heightFigure6C - yFigure6C(d.Counts)-paddingxFigure6C; }).style("fill", "#bcbddc").style("opacity", .65)
      //.text(function(d) { return d; })
    .transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      .attr("x", function(d) { return xFigure6C(d.Values); })
      				.attr("width", xFigure6C.rangeBand())
      				.attr("y", function(d) { return yFigure6C(d.Counts); })
      				.attr("height", function(d) { return heightFigure6C - yFigure6C(d.Counts)-paddingxFigure6C; }).style("fill", "#bcbddc").style("opacity", .65);


textFigure6.exit().attr("class", "exit")
    .transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      //.style("fill", "rgb(153,216,201)")
      .remove();
      

 svgFigure6C.select(".x.axis")
		.transition()
					.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
		.duration(500)
		.call(xAxisFigure6C);
		
svgFigure6C.select(".xaxis_label")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(500);

		
svgFigure6C.select(".r-line")
           // .attr("class", "line")
        .attr("x1", xFigure6C(closestnr) + xFigure6C.rangeBand()/2)
            .attr("x2", xFigure6C(closestnr) + xFigure6C.rangeBand()/2)
            .attr("y1", 0)
            .attr("y2", heightFigure6C).style("stroke-width",2)
      		.style("stroke", "grey").transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(750).attr("x1", xFigure6C(closestnr) + xFigure6C.rangeBand()/2)
            .attr("x2", xFigure6C(closestnr) + xFigure6C.rangeBand()/2)
            .attr("y1", 0)
            .attr("y2", heightFigure6C).style("stroke-width",2)
      		.style("stroke", "grey");	     


      					
}; //update function 2


//update function for the second panel     			
function updatehistogram2(data, value1, value2) {

			
var data = data.filter(function(d) { 
            return (d.Site==value2 & d.Abundance==value1 & d.CHX=="CHX"); //& d.CHX==value3
    
    });

xFigure6nC.domain(data.map(function(d) { return d.Values; }));
yFigure6nC.domain([0, d3.max(data, function(d) { return d.Counts; })]);

var values = function(d) {
					return d.Values;
				};

function closest (num, arr) {
                var curr = arr[0];
                var diff = Math.abs (num - curr);
                for (var val = 0; val < arr.length; val++) {
                    var newdiff = Math.abs (num - arr[val]);
                    if (newdiff < diff) {
                        diff = newdiff;
                        curr = arr[val];
                    }
                }
                return curr;
            }
            
var testnumber=MyApp.thecorr;

var testarray=data.map(function(d) { return d.Values; });
var closestnr=closest(testnumber, testarray);
      					
      					
var textFigure6 = svgFigure6nC.selectAll("rect")
      .data(data);
      
      	 textFigure6.attr("class", "update")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      				.attr("x", function(d) { return xFigure6nC(d.Values); })
      				.attr("width", xFigure6nC.rangeBand())
      				.attr("y", function(d) { return yFigure6nC(d.Counts); })
      				.attr("height", function(d) { return heightFigure6C - yFigure6nC(d.Counts)-paddingxFigure6C; }).style("fill", "#bcbddc");

  textFigure6.enter().append("rect")
      .attr("class", "enter")
      .attr("x", function(d) { return xFigure6nC(d.Values); })
      				.attr("width", xFigure6nC.rangeBand())
      				.attr("y", function(d) { return yFigure6nC(d.Counts); })
      				.attr("height", function(d) { return heightFigure6C - yFigure6nC(d.Counts)-paddingxFigure6C; }).style("fill", "#bcbddc").style("opacity", .65)
      //.text(function(d) { return d; })
    .transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      .attr("x", function(d) { return xFigure6nC(d.Values); })
      				.attr("width", xFigure6nC.rangeBand())
      				.attr("y", function(d) { return yFigure6nC(d.Counts); })
      				.attr("height", function(d) { return heightFigure6C - yFigure6nC(d.Counts)-paddingxFigure6C; }).style("fill", "#bcbddc").style("opacity", .65);


textFigure6.exit().attr("class", "exit")
    .transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      //.style("fill", "rgb(153,216,201)")
      .remove();
      

 svgFigure6nC.select(".x.axis")
		.transition()
					.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
		.duration(250)
		.call(xAxisFigure6nC);
		
svgFigure6nC.select(".xaxislabel")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(250);

		
svgFigure6nC.select(".r-line2")
           // .attr("class", "line")
        .attr("x1", xFigure6nC(closestnr) + xFigure6nC.rangeBand()/2)
            .attr("x2", xFigure6nC(closestnr) + xFigure6nC.rangeBand()/2)
            .attr("y1", 0)
            .attr("y2", heightFigure6C).style("stroke-width",2)
      		.style("stroke", "grey").transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(750).attr("x1", xFigure6nC(closestnr) + xFigure6nC.rangeBand()/2)
            .attr("x2", xFigure6nC(closestnr) + xFigure6nC.rangeBand()/2)
            .attr("y1", 0)
            .attr("y2", heightFigure6C).style("stroke-width",2)
      		.style("stroke", "grey");	     
					
}; //update function 3
		   

//now render the graphs
updateFigure6(data, "GCN", "A"); 
updatehistogram1(data1, "GCN", "A");
updatehistogram2(data1, "GCN", "A");	
d3.selectAll(".Thecorrs")
      .on("change", changeit6);

function changeit6() {
  
  	var testabundance= d3.select('input[name="abundance"]:checked').node().value;
    var testsite= d3.select('input[name="site"]:checked').node().value;
    updateFigure6(data, testabundance, testsite);
   	updatehistogram1(data1, testabundance, testsite);	
   	updatehistogram2(data1, testabundance, testsite);	
};
d3.select("#download6")
		.on("click", function (){
			window.open(thefile6);
		});

}; //end of analyze function

 }; //change form
}, 500); //timeout