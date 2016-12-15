


var MyApp2 = {
    		thecorr: null
		};

	
var marginFigure7 = {top: 60, right: 20, bottom: 50, left: 100},
		paddingFigure7=5,
    	widthFigure7 = 700 - marginFigure7.left - marginFigure7.right,
    	heightFigure7 = 400 - marginFigure7.top - marginFigure7.bottom;
    	
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


var svgFigure7 = d3.select("#CorrelationsReadsandFeatures").append("svg")
    .attr("width", widthFigure7 + marginFigure7.left + marginFigure7.right)
    .attr("height", heightFigure7+2*marginFigure7.top + marginFigure7.bottom)
  .append("g")
    .attr("transform", "translate(" + marginFigure7.left + "," + marginFigure7.top + ")");
 
			
// var xFigure7C = d3.scale.ordinal().rangeRoundBands([0, widthFigure7C], .05);
// var yFigure7C = d3.scale.linear().range([heightFigure7C-paddingFigure7C,paddingFigure7C]);


				svgFigure7.append("g")
      				.attr("class", "x axis")
      				.attr("transform", "translate(0," + heightFigure7 + ")")
      				.call(xAxisFigure7);
    				    				
				svgFigure7.append("text")     
					.attr("class", "xaxis_label1")
					.attr("transform", "translate(" + (widthFigure7 / 2) + " ," + (heightFigure7 + paddingFigure7*8) + ")")
					.style("text-anchor", "middle")
					//.text("Length")
					.style("font-size","16px");
        				
				//y-axis
  				svgFigure7.append("g")
      				.attr("class", "y axis")
      				.call(yAxisFigure7);
	
      			svgFigure7.append("text")    
      				.attr("class", "yaxis_label")
        			.attr("y", heightFigure7 /2 )
        			.attr("x", -paddingFigure7*6)
        			.style("text-anchor", "middle")
        			.attr("transform", "rotate(0)")
        			//.text("FEatg")
        			.style("font-size","16px");
        			
        		svgFigure7.append("text")
					.attr("class", "r-label")
					.attr("y", -5*paddingFigure7 )
        			.attr("x", paddingFigure7*15)
        			.style("text-anchor", "middle")
        			.attr("transform", "rotate(0)")
        			.style("font-size","16px");
        			
        		// now add titles to the axes
        		svgFigure7.append("text")
            		.attr("text-anchor", "middle")  
            		.attr("transform", "translate("+(0-paddingFigure7*8)+","+(heightFigure7/2)+")rotate(-90)")  
            		.text("Selected feature").style("font-size","16px").style("fill","#777777");

       			svgFigure7.append("text")
            		.attr("text-anchor", "middle")  
            		.attr("transform", "translate("+ (widthFigure7/2) +","+(heightFigure7+paddingFigure7*6.5)+")")  
            		.text("Log 10 read counts").style("font-size","16px").style("fill","#777777");
 
	var string="../../Data/";

	var thefile=string.concat("F7_2016_Weinberg_RPF.tsv");


queue()
  .defer(d3.tsv, thefile)
  .await(analyze);

function analyze(error, data, data1) {
  if(error) { console.log(error); }



  	data.forEach(function(d) {
    d.RPF = +d.RPF;
    d.Length = +d.Length;
    d.FEatg = +d.FE_atg;
    d.FEcap = +d.FE_cap;
    d.uATG = +d.uATGs;
    d.utr = +d.utr;
    d.gc = +d.utr_gc;
    d.polyA = +d.polyA;
  });

	

 xFigure7.domain(d3.extent(data, function(d) { return d.RPF; }));
 yFigure7.domain(d3.extent(data, function(d) { return d.Length; }));


function updateFigure7(data, value1, value2) {


var  data = data.map( function (d) { 
	if(value1=="RPF"){choosetheone=d.RPF};
	if(value2=="Length"){choosethetwo=d.Length};
	if(value2=="FEatg"){choosethetwo=d.FEatg};
	if(value2=="FEcap"){choosethetwo=d.FEcap};
	if(value2=="uATG"){choosethetwo=d.uATG};
	if(value2=="utr"){choosethetwo=d.utr};
	if(value2=="gc"){choosethetwo=d.gc};
    if(value2=="polyA"){choosethetwo=d.polyA};
    return { 
      d1: +choosetheone,
      d2: +choosethetwo}; 
});
    
   
    
	xFigure7.domain(d3.extent(data, function(d) { return d.d1; })).nice();
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
      .style("fill-opacity", 0.2)
      .style("fill","rgb(166,206,227)")
      .attr("r", 3)
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
      .style("fill-opacity", 1);
      //.filter(function(d) { return !isNaN(d.d1) && !isNaN(d.d2)});


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
	var y2 =leastSquaresCoeff[0] *x2+ leastSquaresCoeff[1]; 
    MyApp2.thecorr=thecor;
      
    svgFigure7.select(".r-label")			
			.text("Correlation coefficient: r=" + thecor.toFixed(2))
			.style("font-size","16px").style("fill","#777777")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(750);

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

updateFigure7(data, "RPF", "Length"); 

d3.selectAll(".Thecorrs2")
      .on("change", changeit7);

function changeit7() {
    var feature= d3.select('input[name="features"]:checked').node().value;
    updateFigure7(data, "RPF", feature);	
};


}; //end of analyze function




	