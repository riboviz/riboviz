setTimeout(function() {

var MyApp = {
    		thecorr: null
		};

	
var marginFigure7 = {top: 60, right: 20, bottom: 20, left: 100},
		paddingFigure7=5,
    	widthFigure7 = 700 - marginFigure7.left - marginFigure7.right,
    	heightFigure7 = 300 - marginFigure7.top - marginFigure7.bottom;
    	
var marginFigure7C = {top: 20, right: 60, bottom: 20, left:50},
		paddingFigure7C=5,
		paddingxFigure7C=10,
    	widthFigure7C = 400 - marginFigure7C.left - marginFigure7C.right,
    	heightFigure7C = 400 - marginFigure7C.top - marginFigure7C.bottom;
    
var xFigure7 = d3.scale.linear()
    .range([paddingFigure7, widthFigure7]);

var yFigure7 = d3.scale.linear()
    .range([heightFigure7-paddingFigure7, 0]);

var colorFigure7 = d3.scale.category20()
  			.range(["#9ecae1"]);

var xAxisFigure7 = d3.svg.axis()
    .scale(xFigure7)
    .orient("bottom");
    //.ticks(7);

var yAxisFigure7 = d3.svg.axis()
    .scale(yFigure7)
    .orient("left");//.ticks(7);


var svgFigure7 = d3.select("#CorrelationsAbundancesReads").append("svg")
    .attr("width", widthFigure7 + marginFigure7.left + marginFigure7.right)
    .attr("height", heightFigure7+2*marginFigure7.top + marginFigure7.bottom)
  .append("g")
    .attr("transform", "translate(" + marginFigure7.left + "," + marginFigure7.top + ")");
 
			
var xFigure7C = d3.scale.ordinal().rangeRoundBands([0, widthFigure7C], .05);
var yFigure7C = d3.scale.linear().range([heightFigure7C-paddingFigure7C,paddingFigure7C]);

var xAxisFigure7C = d3.svg.axis()
    			.scale(xFigure7C)
    			.orient("bottom");

var yAxisFigure7C = d3.svg.axis()
    			.scale(yFigure7C)
    			.orient("left")
    			.ticks(10);

var svgFigure7C = d3.select("#HistogramsCorrsnonCHX").append("svg")
    			.attr("width", widthFigure7C + marginFigure7C.left + marginFigure7C.right)
    			.attr("height", heightFigure7C + marginFigure7C.top + marginFigure7C.bottom)
  				.append("g")
    			.attr("transform", "translate(" + marginFigure7C.left + "," + marginFigure7C.top + ")");

var xFigure7nC = d3.scale.ordinal().rangeRoundBands([0, widthFigure7C], .05);
var yFigure7nC = d3.scale.linear().range([heightFigure7C-paddingFigure7C,paddingFigure7C]);

var xAxisFigure7nC = d3.svg.axis()
    			.scale(xFigure7nC)
    			.orient("bottom");

var yAxisFigure7nC = d3.svg.axis()
    			.scale(yFigure7nC)
    			.orient("left")
    			.ticks(10);

var svgFigure7nC = d3.select("#HistogramsCorrsCHX").append("svg")
    			.attr("width", widthFigure7C + marginFigure7C.left + marginFigure7C.right)
    			.attr("height", heightFigure7C + marginFigure7C.top + marginFigure7C.bottom)
  				.append("g")
    			.attr("transform", "translate(" + marginFigure7C.left + "," + marginFigure7C.top + ")");


svgFigure7.append("g")
      				.attr("class", "x axis")
      				.attr("transform", "translate(0," + heightFigure7 + ")")
      				.call(xAxisFigure7);
    				    				
				svgFigure7.append("text")      // text label for the x axis
					.attr("class", "xaxis_label1")
					.attr("transform", "translate(" + (widthFigure7 / 2) + " ," + (heightFigure7 + paddingFigure7*8) + ")")
					.style("text-anchor", "middle")
					//.text("Length")
					.style("font-size","13px");
        				
				//y-axis
  				svgFigure7.append("g")
      				.attr("class", "y axis")
      				//.attr("transform", "translate("  paddingFigure7 + ",0)")
      				.call(yAxisFigure7);
	
      			svgFigure7.append("text")      // text label for the y axis
      				.attr("class", "yaxis_label")
        			.attr("y", heightFigure7 /2 )
        			.attr("x", -paddingFigure7*6)
        			.style("text-anchor", "middle")
        			.attr("transform", "rotate(0)")
        			//.text("FEatg")
        			.style("font-size","13px");
        		
        		        			
        		svgFigure7.append("text")
					.attr("class", "r-label")
					.attr("y", -5*paddingFigure7 )
        			.attr("x", paddingFigure7*30)
        			.style("text-anchor", "middle")
        			.attr("transform", "rotate(0)")
        			//.text("FEatg")
        			.style("font-size","13px");
        			
        			svgFigure7.append("line")
					.attr("class", "r-line")
 					.attr("stroke", "black")
 					.style("stroke-width",0.9);
 						
        			
d3.selectAll(".form-control")
.on("change.7", change7);


function change7() {
	console.log("testingthisshit");
	var year7=d3.select("#yearform").node().value;
	console.log(year7);
	var author7=d3.select("#authorform").node().value;
	console.log(author7);
	var thedataset7=d3.select("#dataform").node().value;
	console.log(thedataset7);
	var string7="../../VizData/";
	console.log(string7);
	var thefile7=string7.concat("F7_Temp_Year_",year7,"_Author_",author7, "_Dataset_", thedataset7,"_data.tsv");
console.log(thefile7);

queue()
  .defer(d3.tsv,thefile7)
  .defer(d3.tsv, "../VizData/Figure7HistCounts.tsv")
  .await(analyze);

function analyze(error, data, data1) {
  if(error) { console.log(error); }

  	data.forEach(function(d) {
    d.A = +d.A;
    d.P = +d.P;
    d.E = +d.E;
    d.GCN = +d.tRNA;
    d.tAI = +d.tAI;
    d.RNAseq = +d.RNAseq;
    d.Microarray = +d.Microarray;
  });

	
console.log(data);
xFigure7.domain(d3.extent(data, function(d) { return d.A; }));
yFigure7.domain(d3.extent(data, function(d) { return d.GCN; }));

	
  			
        			

function updateFigure7(data, value1, value2) {


var  data = data.map( function (d) { 
	if(value2=="A"){choosethetwo=d.A};
	if(value2=="E"){choosethetwo=d.E};
	if(value2=="P"){choosethetwo=d.P};
	if(value1=="GCN"){choosetheone=1/d.GCN};
	if(value1=="tAI"){choosetheone=1/d.tAI.toFixed(9)};
	if(value1=="RNAseq"){choosetheone=1/d.RNAseq.toFixed(9)};
	if(value1=="Microarray"){choosetheone=1/d.Microarray.toFixed(9)};
	
    return { 
      d1: +choosetheone,
      d2: +choosethetwo}; 
});
    
    
    
	xFigure7.domain(d3.extent(data, function(d) { return d.d1; }));
  	yFigure7.domain(d3.extent(data, function(d) { return d.d2; }));
  
  // DATA JOIN
  // Join new data with old elements, if any.
  var textFigure7 = svgFigure7.selectAll("circle")
      .data(data);

  // UPDATE
  // Update old elements as needed.
  textFigure7.attr("class", "update")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      .attr("cx", function(d) { return xFigure7(d.d1); })
      .attr("cy", function(d) { return yFigure7(d.d2); })
      .filter(function(d) { return !isNaN(d.d1) && !isNaN(d.d2)});

  // ENTER
  // Create new elements as needed.
  textFigure7.enter().append("circle")
      .attr("class", "enter")
      .attr("cx", function(d) { return xFigure7(d.d1); })
      .attr("cy", function(d) { return yFigure7(d.d2); })
      .style("fill-opacity", 1e-6)
      .attr("r", 4)
      .filter(function(d) { return !isNaN(d.d1) && !isNaN(d.d2)})
      //.text(function(d) { return d; })
    .transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      .attr("cx", function(d) { return xFigure7(d.d1); })
      .attr("cy", function(d) { return yFigure7(d.d2); })
      .style("fill-opacity", 1)
      .filter(function(d) { return !isNaN(d.d1) && !isNaN(d.d2)});


textFigure7.exit().attr("class", "exit")
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
      
    svgFigure7.select(".r-label")			
			.text("corr: " + thecor)//leastSquaresCoeff[2])
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(750);

			
	svgFigure7.select(".r-line")			
			.attr("x1", xFigure7(trendData[0][0]))
  			.attr("y1", yFigure7(trendData[0][1]))
  			.attr("x2",  xFigure7(trendData[0][2]))
  			.attr("y2", yFigure7(trendData[0][3]))
  			.attr("stroke", "green")
 			.style("stroke-width",0.9)
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(750);
			
			
      // now should really update the axes
      //Update X axis
	svgFigure7.select(".x.axis")
		.transition()
					.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
		.duration(250)
		.call(xAxisFigure7);
					
	//Update Y axis
	svgFigure7.select(".y.axis")
			.transition()
						.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(250)
			.call(yAxisFigure7);
						
		svgFigure7.select(".xaxis_label1")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(250);
            // .text("d1");
						
		svgFigure7.select(".yaxis_label")
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
	console.log("testthismf");
//now move to the next panel   			
data1.forEach(function(d) {
        			d.Site = d.Site;
        			d.Abundance = d.Abundance;
        			d.Counts=+d.Counts;
        			d.Values=+d.Values;
        			d.CHX=d.CHX
    			});
//console.log(data1);     			
svgFigure7C.append("g").append("line")
		.attr("class", "r-line").attr("stroke", "black")
 					.style("stroke-width",0.9);		
    					
svgFigure7C.append("g")
      				.attr("class", "x axis")
      				.attr("transform", "translate(0," + heightFigure7C + ")")
      				.call(xAxisFigure7C);
    
svgFigure7C.append("text")
    .attr("class", "xaxis_label")
      					.style("text-anchor", "end")
      					.attr("dx", "-.8em")
      					.attr("dy", "-.55em")
      					.attr("transform", "rotate(0)" )
      					.attr("transform", "translate(17," +15+ ")")
      					.style("font-size","12px");
      					
 
 svgFigure7nC.append("g").append("line")
		.attr("class", "r-line2").attr("stroke", "black")
 					.style("stroke-width",0.9);		
    					
svgFigure7nC.append("g")
      				.attr("class", "x axis")
      				.attr("transform", "translate(0," + heightFigure7C + ")")
      				.call(xAxisFigure7nC);
    
svgFigure7nC.append("text")
    .attr("class", "xaxislabel")
      					.style("text-anchor", "end")
      					.attr("dx", "-.8em")
      					.attr("dy", "-.55em")
      					.attr("transform", "rotate(0)" )
      					.attr("transform", "translate(17," +15+ ")")
      					.style("font-size","12px");
      					
      					    					
//update function for the second panel     			
function updatehistogram1(data, value1, value2) {

			
var data = data.filter(function(d) { 
            return (d.Site==value2 & d.Abundance==value1 & d.CHX=="non-CHX"); //& d.CHX==value3
    
    });
//console.log(data);
xFigure7C.domain(data.map(function(d) { return d.Values; }));
yFigure7C.domain([0, d3.max(data, function(d) { return d.Counts; })]);

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
//console.log(testnumber);
var testarray=data.map(function(d) { return d.Values; });
var closestnr=closest(testnumber, testarray);
//console.log(testarray);
//console.log(closestnr);      					
      					
var textFigure7 = svgFigure7C.selectAll("rect")
      .data(data);
      
      	 textFigure7.attr("class", "update")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      				.attr("x", function(d) { return xFigure7C(d.Values); })
      				.attr("width", xFigure7C.rangeBand())
      				.attr("y", function(d) { return yFigure7C(d.Counts); })
      				.attr("height", function(d) { return heightFigure7C - yFigure7C(d.Counts)-paddingxFigure7C; }).style("fill", "rgb(153,216,201)");

  textFigure7.enter().append("rect")
      .attr("class", "enter")
      .attr("x", function(d) { return xFigure7C(d.Values); })
      				.attr("width", xFigure7C.rangeBand())
      				.attr("y", function(d) { return yFigure7C(d.Counts); })
      				.attr("height", function(d) { return heightFigure7C - yFigure7C(d.Counts)-paddingxFigure7C; }).style("fill", "rgb(153,216,201)")
      //.text(function(d) { return d; })
    .transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      .attr("x", function(d) { return xFigure7C(d.Values); })
      				.attr("width", xFigure7C.rangeBand())
      				.attr("y", function(d) { return yFigure7C(d.Counts); })
      				.attr("height", function(d) { return heightFigure7C - yFigure7C(d.Counts)-paddingxFigure7C; }).style("fill", "rgb(153,216,201)");


textFigure7.exit().attr("class", "exit")
    .transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      //.style("fill", "rgb(153,216,201)")
      .remove();
      

 svgFigure7C.select(".x.axis")
		.transition()
					.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
		.duration(250)
		.call(xAxisFigure7C);
		
svgFigure7C.select(".xaxis_label")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(250);

		
svgFigure7C.select(".r-line")
           // .attr("class", "line")
        .attr("x1", xFigure7C(closestnr) + xFigure7C.rangeBand()/2)
            .attr("x2", xFigure7C(closestnr) + xFigure7C.rangeBand()/2)
            .attr("y1", 0)
            .attr("y2", heightFigure7C).style("stroke-width",2)
      		.style("stroke", "grey").transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(750).attr("x1", xFigure7C(closestnr) + xFigure7C.rangeBand()/2)
            .attr("x2", xFigure7C(closestnr) + xFigure7C.rangeBand()/2)
            .attr("y1", 0)
            .attr("y2", heightFigure7C).style("stroke-width",2)
      		.style("stroke", "grey");	     
					
}; //update function 2


//update function for the second panel     			
function updatehistogram2(data, value1, value2) {

			
var data = data.filter(function(d) { 
            return (d.Site==value2 & d.Abundance==value1 & d.CHX=="CHX"); //& d.CHX==value3
    
    });
//console.log(data);
xFigure7nC.domain(data.map(function(d) { return d.Values; }));
yFigure7nC.domain([0, d3.max(data, function(d) { return d.Counts; })]);

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
//console.log(testnumber);
var testarray=data.map(function(d) { return d.Values; });
var closestnr=closest(testnumber, testarray);
//console.log(testarray);
//console.log(closestnr);      					
      					
var textFigure7 = svgFigure7nC.selectAll("rect")
      .data(data);
      
      	 textFigure7.attr("class", "update")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      				.attr("x", function(d) { return xFigure7nC(d.Values); })
      				.attr("width", xFigure7nC.rangeBand())
      				.attr("y", function(d) { return yFigure7nC(d.Counts); })
      				.attr("height", function(d) { return heightFigure7C - yFigure7nC(d.Counts)-paddingxFigure7C; }).style("fill", "rgb(153,216,201)");

  textFigure7.enter().append("rect")
      .attr("class", "enter")
      .attr("x", function(d) { return xFigure7nC(d.Values); })
      				.attr("width", xFigure7nC.rangeBand())
      				.attr("y", function(d) { return yFigure7nC(d.Counts); })
      				.attr("height", function(d) { return heightFigure7C - yFigure7nC(d.Counts)-paddingxFigure7C; }).style("fill", "rgb(153,216,201)")
      //.text(function(d) { return d; })
    .transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      .attr("x", function(d) { return xFigure7nC(d.Values); })
      				.attr("width", xFigure7nC.rangeBand())
      				.attr("y", function(d) { return yFigure7nC(d.Counts); })
      				.attr("height", function(d) { return heightFigure7C - yFigure7nC(d.Counts)-paddingxFigure7C; }).style("fill", "rgb(153,216,201)");


textFigure7.exit().attr("class", "exit")
    .transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      //.style("fill", "rgb(153,216,201)")
      .remove();
      

 svgFigure7nC.select(".x.axis")
		.transition()
					.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
		.duration(250)
		.call(xAxisFigure7nC);
		
svgFigure7nC.select(".xaxislabel")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(250);

		
svgFigure7nC.select(".r-line2")
           // .attr("class", "line")
        .attr("x1", xFigure7nC(closestnr) + xFigure7nC.rangeBand()/2)
            .attr("x2", xFigure7nC(closestnr) + xFigure7nC.rangeBand()/2)
            .attr("y1", 0)
            .attr("y2", heightFigure7C).style("stroke-width",2)
      		.style("stroke", "grey").transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(750).attr("x1", xFigure7nC(closestnr) + xFigure7nC.rangeBand()/2)
            .attr("x2", xFigure7nC(closestnr) + xFigure7nC.rangeBand()/2)
            .attr("y1", 0)
            .attr("y2", heightFigure7C).style("stroke-width",2)
      		.style("stroke", "grey");	     
					
}; //update function 3
		   

//now render the graphs
updateFigure7(data, "GCN", "A"); 
//updatehistogram1(data1, "GCN", "A");
//updatehistogram2(data1, "GCN", "A");	
d3.selectAll(".Thecorrs")
      .on("change", changeit7);

function changeit7() {
  
  	var testabundance= d3.select('input[name="abundance"]:checked').node().value;
    var testsite= d3.select('input[name="site"]:checked').node().value;
    updateFigure7(data, testabundance, testsite);
   	//updatehistogram1(data1, testabundance, testsite);	
   //	updatehistogram2(data1, testabundance, testsite);	
};


}; //end of analyze function

}; //change form
}, 10); //timeout