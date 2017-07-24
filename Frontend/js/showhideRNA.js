$(document).ready(function() {
toggleFields2();

function toggleFields2() {
     if ($("#yearform2").val()==0){
         $("#author2").hide();
         $("#dataset2").hide();
         $("#hideallatfirst").hide();
     } else {
         $("#author2").show();
         $("#dataset2").show();
         $("#hideallatfirst").show()
         }
  }
function isEmpty(str) {
    return (!str || 0 === str.length);
}

d3.tsv("../../Data/AllDataRNA.tsv", function(error, data) {
  	
  		if (error) throw error;
  		data.forEach(function(d) {
  			d.Year = +d.Year;
  			d.Author = d.Author;
  			d.Dataset = d.Dataset;
  			d.mRNA= +d.mRNA;
  		});

		$("#yearform2").change(function () {
         toggleFields2();
         var yval2 = $(this).val();
         
         var datay2 = data.filter(function(d, key) { 
            return ((d.Year==yval2) );
    
    	});
          
    	var authors2 =d3.map(datay2, function(d){return d.Author;}).keys();
        var dsets2 =d3.map(datay2, function(d){return d.Dataset;}).keys();
         	var acc="";
         	for (i=0; i<authors2.length; i++){
         		acc=acc+"<option value="+authors2[i]+">"+authors2[i]+"</option>";
         		$("#authorform2").html(acc);
         	}

		var dataya12 = datay2.filter(function(d, key) { 
            return ((d.Author==authors2[0]) );
    
    	});

        var dsets2 =d3.map(dataya12, function(d){return d.Dataset;}).keys();	
         var acc2="";
         	for (i=0; i<dsets2.length; i++){
         		acc2=acc2+"<option value="+dsets2[i]+">"+dsets2[i]+"</option>";
         		$("#dataform2").html(acc2);
         	}
        var str=$("#dataform2").val();
  		if (isEmpty(str.match(/RNA/g)) && isEmpty(str.match(/Dynabeads/g))){
  		$("#hideformRNA").show();
  		} else {
  		$("#hideformRNA").hide();
  		}	

   });
   
  		$("#authorform2").change(function () {
        var aval2 = $(this).val();
         
         var dataya2 = data.filter(function(d, key) { 
            return ((d.Author==aval2) );
    
    		});
         
         var datasets2 =d3.map(dataya2, function(d){return d.Dataset;}).keys();

         	var acc2="";
         	for (i=0; i<datasets2.length; i++){
         		acc2=acc2+"<option value="+datasets2[i]+">"+datasets2[i]+"</option>";
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

