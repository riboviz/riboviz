setTimeout(function() {


var MyApp2 = {
    		thecorr: null
		};

	
var marginFigure8 = {top: 60, right: 20, bottom: 20, left: 100},
		paddingFigure8=5,
    	widthFigure8 = 700 - marginFigure8.left - marginFigure8.right,
    	heightFigure8 = 400 - marginFigure8.top - marginFigure8.bottom;
    	
var marginFigure8C = {top: 20, right: 60, bottom: 20, left:50},
		paddingFigure8C=5,
		paddingxFigure8C=10,
    	widthFigure8C = 400 - marginFigure8C.left - marginFigure8C.right,
    	heightFigure8C = 400 - marginFigure8C.top - marginFigure8C.bottom;
    
var xFigure8 = d3.scale.linear()
    .range([paddingFigure8, widthFigure8]);

var yFigure8 = d3.scale.linear()
    .range([heightFigure8-paddingFigure8, 0]);

var colorFigure8 = d3.scale.category20()
  			.range(["#9ecae1"]);

var xAxisFigure8 = d3.svg.axis()
    .scale(xFigure8)
    .orient("bottom");
    //.ticks(7);

var yAxisFigure8 = d3.svg.axis()
    .scale(yFigure8)
    .orient("left");//.ticks(7);


var svgFigure8 = d3.select("#CorrelationsReadsandFeatures").append("svg")
    .attr("width", widthFigure8 + marginFigure8.left + marginFigure8.right)
    .attr("height", heightFigure8+2*marginFigure8.top + marginFigure8.bottom)
  .append("g")
    .attr("transform", "translate(" + marginFigure8.left + "," + marginFigure8.top + ")");
 
			
var xFigure8C = d3.scale.ordinal().rangeRoundBands([0, widthFigure8C], .05);
var yFigure8C = d3.scale.linear().range([heightFigure8C-paddingFigure8C,paddingFigure8C]);

var xAxisFigure8C = d3.svg.axis()
    			.scale(xFigure8C)
    			.orient("bottom");

var yAxisFigure8C = d3.svg.axis()
    			.scale(yFigure8C)
    			.orient("left")
    			.ticks(10);

var svgFigure8C = d3.select("#HistogramsCorrsnonCHX").append("svg")
    			.attr("width", widthFigure8C + marginFigure8C.left + marginFigure8C.right)
    			.attr("height", heightFigure8C + marginFigure8C.top + marginFigure8C.bottom)
  				.append("g")
    			.attr("transform", "translate(" + marginFigure8C.left + "," + marginFigure8C.top + ")");

svgFigure8.append("g")
      				.attr("class", "x axis")
      				.attr("transform", "translate(0," + heightFigure8 + ")")
      				.call(xAxisFigure8);
    				    				
				svgFigure8.append("text")      // text label for the x axis
					.attr("class", "xaxis_label1")
					.attr("transform", "translate(" + (widthFigure8 / 2) + " ," + (heightFigure8 + paddingFigure8*8) + ")")
					.style("text-anchor", "middle")
					//.text("Length")
					.style("font-size","13px");
        				
				//y-axis
  				svgFigure8.append("g")
      				.attr("class", "y axis")
      				//.attr("transform", "translate("  paddingFigure8 + ",0)")
      				.call(yAxisFigure8);
	
      			svgFigure8.append("text")      // text label for the y axis
      				.attr("class", "yaxis_label")
        			.attr("y", heightFigure8 /2 )
        			.attr("x", -paddingFigure8*6)
        			.style("text-anchor", "middle")
        			.attr("transform", "rotate(0)")
        			//.text("FEatg")
        			.style("font-size","13px");
        			
        		svgFigure8.append("text")
					.attr("class", "r-label")
					.attr("y", -5*paddingFigure8 )
        			.attr("x", paddingFigure8*30)
        			.style("text-anchor", "middle")
        			.attr("transform", "rotate(0)")
        			//.text("FEatg")
        			.style("font-size","13px");
        			
        		// 	svgFigure8.append("line")
// 					.attr("class", "r-line")
//  					.attr("stroke", "black")
//  					.style("stroke-width",0.9);


d3.selectAll(".form-control")
.on("change.8", change8);


function change8() {

	var year8=d3.select("#yearform").node().value;
	var author8=d3.select("#authorform").node().value;
	var thedataset8=d3.select("#dataform").node().value;
	var string="../Data/";

	var thefile=string.concat("F8_Year_",year8,"_Author_",author8, "_Dataset_", thedataset8,"_data.tsv");


queue()
  .defer(d3.tsv, thefile)
  .await(analyze);

function analyze(error, data, data1) {
  if(error) { console.log(error); }


  	data.forEach(function(d) {
    d.RPF = +d.RPF;
    d.mRNA = +d.mRNA;
    d.Length = +d.Length;
    d.FEatg = +d.FE_atg;
    d.FEcap = +d.FE_cap;
    d.uATG = +d.uATGs;
    d.utr = +d.utr;
    d.gc = +d.utr_gc;
    d.polyA = +d.polyA;
  });

console.log(data);	

xFigure8.domain(d3.extent(data, function(d) { return d.RPF; }));
yFigure8.domain(d3.extent(data, function(d) { return d.Length; }));

	

function updateFigure8(data, value1, value2) {


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
    
   
    
	xFigure8.domain(d3.extent(data, function(d) { return d.d1; }));
  	yFigure8.domain(d3.extent(data, function(d) { return d.d2; }));
  
  // DATA JOIN
  // Join new data with old elements, if any.
  var textFigure8 = svgFigure8.selectAll("circle")
      .data(data);

  // UPDATE
  // Update old elements as needed.
  textFigure8.attr("class", "update")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      .attr("cx", function(d) { return xFigure8(d.d1); })
      .attr("cy", function(d) { return yFigure8(d.d2); })
      .filter(function(d) { return !isNaN(d.d1) && !isNaN(d.d2)});

  // ENTER
  // Create new elements as needed.
  textFigure8.enter().append("circle")
      .attr("class", "enter")
      .attr("cx", function(d) { return xFigure8(d.d1); })
      .attr("cy", function(d) { return yFigure8(d.d2); })
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
      .attr("cx", function(d) { return xFigure8(d.d1); })
      .attr("cy", function(d) { return yFigure8(d.d2); })
      .style("fill-opacity", 1);
      //.filter(function(d) { return !isNaN(d.d1) && !isNaN(d.d2)});


textFigure8.exit().attr("class", "exit")
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
      
    svgFigure8.select(".r-label")			
			.text("Correlation coefficient: r=" + thecor.toFixed(2))
			.style("font-size","16px")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(750);

      //Update X axis
	svgFigure8.select(".x.axis")
		.transition()
					.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
		.duration(250)
		.call(xAxisFigure8);
					
	//Update Y axis
	svgFigure8.select(".y.axis")
			.transition()
						.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(250)
			.call(yAxisFigure8);
						
		svgFigure8.select(".xaxis_label1")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(250);
            // .text("d1");
						
		svgFigure8.select(".yaxis_label")
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

updateFigure8(data, "RPF", d3.select('input[name="features"]:checked').node().value); 

d3.selectAll(".Thecorrs2")
      .on("change", changeit8);

function changeit8() {
    var feature= d3.select('input[name="features"]:checked').node().value;
    updateFigure8(data, "RPF", feature);	
};
d3.select("#download6")
		.on("click", function (){
			window.open(thefile );
		});
		
}; //end of analyze function

}; //change form

}, 30); //timeout
	