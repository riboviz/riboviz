setTimeout(function() {


var MyApp2 = {
    		thecorr: null
		};

	
var marginFigureE = {top: 60, right: 20, bottom: 20, left: 100},
		paddingFigureE=5,
    	widthFigureE = 700 - marginFigureE.left - marginFigureE.right,
    	heightFigureE = 400 - marginFigureE.top - marginFigureE.bottom;
    	
var marginFigureEC = {top: 20, right: 60, bottom: 20, left:50},
		paddingFigureEC=5,
		paddingxFigureEC=10,
    	widthFigureEC = 400 - marginFigureEC.left - marginFigureEC.right,
    	heightFigureEC = 400 - marginFigureEC.top - marginFigureEC.bottom;
    
var xFigureE = d3.scale.linear()
     .range([paddingFigureE, widthFigureE]);

var yFigureE = d3.scale.linear()
    .range([heightFigureE-paddingFigureE, 0]);

var colorFigureE = d3.scale.category20()
  			.range(["#9ecae1"]);

var xAxisFigureE = d3.svg.axis()
    .scale(xFigureE)
    .orient("bottom");
    //.ticks(7);

var yAxisFigureE = d3.svg.axis()
    .scale(yFigureE)
    .orient("left");//.ticks(7);

var tipFigureE = d3.select("#CorrelationsRNAandRPF")
      .append('div')
      .attr('class', 'tooltip')
      .style("opacity", 0);


var svgFigureE = d3.select("#CorrelationsRNAandRPF").append("svg")
    .attr("width", widthFigureE + marginFigureE.left + marginFigureE.right)
    .attr("height", heightFigureE+2*marginFigureE.top + marginFigureE.bottom)
  .append("g")
    .attr("transform", "translate(" + marginFigureE.left + "," + marginFigureE.top + ")");
 
		

svgFigureE.append("g")
      				.attr("class", "x axis")
      				.attr("transform", "translate(0," + heightFigureE + ")")
      				.call(xAxisFigureE);
    				
    				
    		    //x-axis line
    			svgFigureE.selectAll(".x.axis")	
  					.append("rect")
  	   				.attr("width", widthFigureE)
  	   				.attr("height",1)
  	   				.attr("fill","#000")
    				    				
				svgFigureE.append("text")      // text label for the x axis
					.attr("class", "xaxis_label1")
					.attr("transform", "translate(" + (widthFigureE / 2) + " ," + (heightFigureE + paddingFigureE*8) + ")")
					.style("text-anchor", "middle")
					//.text("Length")
					.style("font-size","16px");
        				
				//y-axis
  				svgFigureE.append("g")
      				.attr("class", "y axis")
      				//.attr("transform", "translate("  paddingFigureE + ",0)")
      				.call(yAxisFigureE);
	
      			svgFigureE.append("text")      // text label for the y axis
      				.attr("class", "yaxis_label")
        			.attr("y", heightFigureE /2 )
        			.attr("x", -paddingFigureE*6)
        			.style("text-anchor", "middle")
        			.attr("transform", "rotate(0)")
        			//.text("FEatg")
        			.style("font-size","16px");
        			
        		svgFigureE.append("text")
					.attr("class", "r-label")
					.attr("y", -5*paddingFigureE )
        			.attr("x", paddingFigureE*15)
        			.style("text-anchor", "middle")
        			.attr("transform", "rotate(0)")
        			//.text("FEatg")
        			.style("font-size","16px");
        			
// now add titles to the axes
        		svgFigureE.append("text")
            		.attr("text-anchor", "middle")  
            		.attr("transform", "translate("+(0-paddingFigureE*8.5)+","+(heightFigureE/2)+")rotate(-90)")  
            		.text("Log10 RPF").style("font-size","16px").style("fill","#777777");

       			svgFigureE.append("text")
            		.attr("text-anchor", "middle")  
            		.attr("transform", "translate("+ (widthFigureE/2) +","+(heightFigureE+paddingFigureE*6.9)+")")  
            		.text("Log10 RNA").style("font-size","16px").style("fill","#777777");


d3.selectAll(".form-control")
.on("change.1", change1);


function change1() {

	var year1=d3.select("#yearform").node().value;
	var author1=d3.select("#authorform").node().value;
	var thedataset1=d3.select("#dataform").node().value;
	var year1RNA=d3.select("#yearform2").node().value;
	var author1RNA=d3.select("#authorform2").node().value;
	var thedataset1RNA=d3.select("#dataform2").node().value;
	var string="../Data/";

	var thefile=string.concat("F7_",year1,"_",author1, "_", thedataset1,".tsv");
	var thefileRNA=string.concat("F7_",year1RNA,"_",author1RNA, "_", thedataset1RNA,".tsv");
	

queue()
  .defer(d3.tsv, thefile)
  .defer(d3.tsv, thefileRNA)
  .await(analyze);

function analyze(error, data, data1) {
  if(error) {  }


  	data.forEach(function(d) {
    d.RPF = +d.Data;
    d.Length = +d.Length;
    d.FEatg = +d.FE_atg;
    d.FEcap = +d.FE_cap;
    d.uATG = +d.uATGs;
    d.utr = +d.utr;
    d.gc = +d.utr_gc;
    d.polyA = +d.polyA;
    d.ORF=d.ORF;
  });

	data1.forEach(function(d) {
    d.RNA = +d.Data;
    d.Length = +d.Length;
    d.FEatg = +d.FE_atg;
    d.FEcap = +d.FE_cap;
    d.uATG = +d.uATGs;
    d.utr = +d.utr;
    d.gc = +d.utr_gc;
    d.polyA = +d.polyA;
    d.Corr=+d.Corr;
    d.ORF=d.ORF;
  });
  
  
  var  dataRPF = data.map(function(d) { return d.RPF; });
  var  dataRNA = data1.map(function(d) { return d.RNA; }); 
  //var te =dataRPF.map(function(n, i) { return n / dataRNA[i]; });
 
 
 // xFigureE.domain(d3.extent(te));
//  yFigureE.domain(d3.extent(data, function(d) { return d.Length; }));
// 
// 	



var  data = data.map( function (d, i) { 
	choosetheone=dataRNA[i];
	choosethetwo=dataRPF[i];
    return { 
      d1: +choosetheone,
      d2: +choosethetwo,
      ORF: d.ORF,
      Corr: +d.Corr}; 
});

    var opa = d3.scale.linear()
            .domain([0, 1])
            .range([0.4, 1]);
    var opaColor = d3.scale.linear()
            .domain([0, 1])
            .range(["rgb(239,138,98)", "rgb(166,206,227)"]);
            
	xFigureE.domain(d3.extent(data, function(d) { return d.d1; }));
	//xFigureE.domain([0, d3.max(data, function(d) { return d.d1; })]);
  	yFigureE.domain(d3.extent(data, function(d) { return d.d2; }));
  	//yFigureE.domain([0, d3.max(data, function(d) { return d.d2; })]);
  
  // DATA JOIN
  // Join new data with old elements, if any.
  var textFigureE = svgFigureE.selectAll("circle")
      .data(data);

  // UPDATE
  // Update old elements as needed.
  textFigureE.attr("class", "update")
  			.transition()
  			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      .attr("cx", function(d) { return xFigureE(d.d1); })
      .attr("cy", function(d) { return yFigureE(d.d2); })
      .filter(function(d) { return !isNaN(d.d1) && !isNaN(d.d2)});

  // ENTER
  // Create new elements as needed.
  textFigureE.enter().append("circle")
      .attr("class", "enter")
      .attr("cx", function(d) { return xFigureE(d.d1); })
      .attr("cy", function(d) { return yFigureE(d.d2); })
      .style("fill-opacity", function(d) { return opa(d.Corr); })
      .style("fill", function(d) { return opaColor(d.Corr); })
      .attr("r", 3)
      .filter(function(d) { return !isNaN(d.d1) && !isNaN(d.d2)})
      .on("mouseover", function(d) {		
            tipFigureE.transition()		
                .duration(20)		
                .style("opacity", .9);		
            tipFigureE.html(d.ORF)	
                .style("left",xFigureE(d.d1) +70+ 'px')		
                .style("top", yFigureE(d.d2) +20 + 'px');	
            })					
        .on("mouseout", function(d) {		
            tipFigureE.transition()		
                .delay(50)		
                .style("opacity", 0);	
        })
    .transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      .attr("cx", function(d) { return xFigureE(d.d1); })
      .attr("cy", function(d) { return yFigureE(d.d2); })
      .style("fill-opacity", function(d) { return opa(d.Corr); })
      .style("fill", function(d) { return opaColor(d.Corr); });

      //.filter(function(d) { return !isNaN(d.d1) && !isNaN(d.d2)});


textFigureE.exit().attr("class", "exit")
    .transition()
    			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
      .duration(750)
      .style("fill-opacity", 1e-6)
      .remove();
      
	  
	var datafiltered = data.filter(function(d, key) { 
            return (d.Corr==1 )}); 
    var datafiltered = datafiltered.filter(function(d, key) { 
            return (!isNaN(d.d1) && !isNaN(d.d2))});
	  var xSeries = datafiltered.map(function(d) { return d.d1; });
	  xSeriesreg=xSeries.sort(function(a, b) {
  				return Number(a) - Number(b);
				});
	var xSeriesreal = datafiltered.map(function(d) { return d.d1; });
	var ySeries = datafiltered.map(function(d) { return d.d2; });
	thecorv=[xSeriesreal, ySeries];
	var thecor=pearsonCorrelation(thecorv,0,1);
	var leastSquaresCoeff = leastSquares(xSeriesreal, ySeries);
	var x1 = xSeriesreg[0];
	var y1 = leastSquaresCoeff[0] *x1+ leastSquaresCoeff[1];
	var x2 = xSeriesreg[xSeriesreg.length - 1];
	var y2 =leastSquaresCoeff[0] *x2+ leastSquaresCoeff[1]; 
    MyApp2.thecorr=thecor;

    svgFigureE.select(".r-label")			
			.text("Pearson correlation coefficient: r=" + thecor.toFixed(2))
			.style("font-size","16px").style("fill","#777777")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(750);

      //Update X axis
	svgFigureE.select(".x.axis")
		.transition()
					.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
		.duration(250)
		.call(xAxisFigureE);
					
	//Update Y axis
	svgFigureE.select(".y.axis")
			.transition()
						.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(250)
			.call(yAxisFigureE);
						
		svgFigureE.select(".xaxis_label1")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(250);
            // .text("d1");
						
		svgFigureE.select(".yaxis_label")
			.transition()
			.ease("linear")
			.delay(function(d, i) {
					return i / data.length * 500;
			})
			.duration(250);
             //.text("d2");
};

function cleaner(arr) {
   return arr.filter(function(item){ 
      return typeof item == "string" || (typeof item == "number" && item);
              /** Any string**/        /** Numbers without NaN & 0 **/
   });
}


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


}; //end of analyze function



}, 500); //timeout
	