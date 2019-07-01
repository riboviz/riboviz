$(document).ready(function() {
toggleFields();

function toggleFields() {
     if ($("#yearform").val()==0){
         $("#author").hide();
         $("#dataset").hide();
         $("#hideallatfirst").hide();
     } else {
         $("#author").show();
         $("#dataset").show();
         $("#hideallatfirst").show()
         }
  }
function isEmpty(str) {
    return (!str || 0 === str.length);
}

d3.tsv("../../Data/AllData.tsv", function(error, data) {
  	
  		if (error) throw error;
  		data.forEach(function(d) {
  			d.Year = +d.Year;
  			d.Author = d.Author;
  			d.Dataset = d.Dataset;
  			d.mRNA= +d.mRNA;
  		});

		$("#yearform").change(function () {
         toggleFields();
         var yval = $(this).val();
         
         var datay = data.filter(function(d, key) { 
            return ((d.Year==yval) );
    
    	});
          
    	var authors =d3.map(datay, function(d){return d.Author;}).keys();
        var dsets =d3.map(datay, function(d){return d.Dataset;}).keys();
         	var acc="";
         	for (i=0; i<authors.length; i++){
         		acc=acc+"<option value="+authors[i]+">"+authors[i]+"</option>";
         		$("#authorform").html(acc);
         	}

		var dataya1 = datay.filter(function(d, key) { 
            return ((d.Author==authors[0]) );
    
    	});

        var dsets =d3.map(dataya1, function(d){return d.Dataset;}).keys();	
         var acc2="";
         	for (i=0; i<dsets.length; i++){
         		acc2=acc2+"<option value="+dsets[i]+">"+dsets[i]+"</option>";
         		$("#dataform").html(acc2);
         	}
        var str=$("#dataform").val();
  		if (isEmpty(str.match(/RNA/g)) && isEmpty(str.match(/Dynabeads/g))){
  		$("#hideformRNA").show();
  		} else {
  		$("#hideformRNA").hide();
  		}	

   });
   
  		$("#authorform").change(function () {
        var aval = $(this).val();
         
         var dataya = data.filter(function(d, key) { 
            return ((d.Author==aval) );
    
    		});
         
         var datasets =d3.map(dataya, function(d){return d.Dataset;}).keys();

         	var acc2="";
         	for (i=0; i<datasets.length; i++){
         		acc2=acc2+"<option value="+datasets[i]+">"+datasets[i]+"</option>";
         		$("#dataform").html(acc2);
         	}
        
  		});
  		
  		$("#dataform").change(function () {
  		var str=$("#dataform").val();
  		if (isEmpty(str.match(/RNA/g)) && isEmpty(str.match(/Dynabeads/g))){
  		$("#hideformRNA").show();
  		} else {
  		$("#hideformRNA").hide();
  		}
  		
  		});	
  		
    		
 });  
 
});

