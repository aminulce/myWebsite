	$(document).ready(function(){
	var max_level = 5; //maximum input boxes allowed
	var visc_holder = $("#visc-holder"); //Fields wrapper
	var add_visc = $("#add_field"); //Add button ID
	const unique = (value, index, self) => {//Deletes duplicate values from arrays
		return self.indexOf(value) === index;
	}
	function add_field(holder, field, level, max_level, placeholder) {
		
		if(level < max_level)
		{
			$(holder).append('<input name="'+field+x+'" type="text" id="'+field+x+'" placeholder="'+placeholder+'" class = "calc" value="">' ); //add input box
		}
	}
	
	function arrstr2posnum (num1)//string array to positive numeric values
	{
		for(var i=0; i<num1.length; i++) 
		{ 
			if($.isNumeric(num1[i]))
			{
				num1[i] = parseFloat(num1[i]); 
				if (num1[i]<=0)
					return -1;
			}
			else{
				return 0;
			}
			
		} 
		return num1;
	}
	function arrstr2negnum (num1)//string array to positive numeric values
	{
		if(!num1)
			return -2;
		else
		{
			for(var i=0; i<num1.length; i++) 
			{ 
				if($.isNumeric(num1[i]))
				{
					num1[i] = parseFloat(num1[i]); 
					if (num1[i]>0)
						return -1;
				}
				else{
					return 0;
				}
				
			} 
		}
		return num1;
	}
	function dsrvalidate(temps,vals,name)//Validation function for dsr values and temperatures
	{
		var errormsg = '';
		if (temps==0)
			errormsg += '<li>'+name+' test temperatures should be numeric</li>';
		if (temps==-1)
			errormsg += '<li>'+name+' test temperatures should be positive</li>';
		if (vals==0)
			errormsg += '<li>'+name+' DSR value should be numeric</li>';
		if (vals==-1)
			errormsg += '<li>'+name+' DSR value should be positive</li>';
		if (vals.length != temps.length && vals!=0 && vals !=-1 && temps!=0 && temps !=-1 && typeof temps.length !== "undefined" && typeof vals.length !== "undefined")
			errormsg += '<li>'+'Number of '+name+' test temperature and number of '+name+' DSR values should be same\nYou provided '+temps.length+' temperatures and '+vals.length+' DSR values</li>';
		return errormsg;
	}
	function validateSingleValue (param,name)
	{
		var errormsg = '';
		if (param==-1)
			errormsg += '<li>'+name+' should be a positive number</li>';
		else if (param==0)
			errormsg += '<li>'+name+' should be a numeric value</li>';
		return errormsg;
	}
	function checkVals(values,ls,checkparam)//checks values depending on the checkparams and ls determines largeer-then or smaller-than direction
	{
		var index_array = []; //this variable saves the indices of the passed values
		var index = 0; //index of the index_array
		if (ls=='lt')//less than
		{
			for(var i=0; i<values.length; i++) 
			{ 
				if(values[i]<checkparam)
				{
					index_array[index]=i; 
					index++;
				}
			}
		}
		else if (ls=='gt')//greater than
		{
			for(var i=0; i<values.length; i++) 
			{ 
				if(values[i]>checkparam)
				{
					index_array[index]=i; 
					index++;
				}
			}
		}
		return index_array;
	}
	function sm_indices(svalues,mvalues,cstemps, dttemps, dtvals)//checks values depending on the checkparams and ls determines largeer-then or smaller-than direction
	{
		var sindex_array = [], mindex_array = []; //this variable saves the indices of the passed values
		var sindex = 0, mindex = 0; //index of the index_array
		for(var i=0; i<svalues.length; i++) 
		{ 
			if(svalues[i]<300 && svalues[i]>0.3)
			{
				sindex_array[sindex]=i; 
				sindex++;
			}
			else if(svalues[i]>300 && svalues[i]<600 && mvalues[i]>0.3 && dttemps!=-2 && dtvals!=-2)
			{
				
				for(var j=0; j<dtvals.length; j++) 
				{
					if (dtvals[j]>1)
					{
						sindex_array[sindex]=i;
						sindex++;
					}
				}
			}
		}
		
		return sindex_array.filter(unique);
	}

	function pg_inter_temp(a,b) //get intermediate temperatures for passed high and low temps
	{
		var index = 0;
		var intermediate = [];
		for (var i=0; i<a.length; i++)
			for (var j=0; j<b.length; j++)
			{
				intermediate[index] = ((a[i]+b[j]-10)/2)+4;
				index++;
			}
		return intermediate.filter(unique);
	}
	
	function index2val(index,arr)
	{
		var arr2 = [];
		for (var i=0; i<index.length; i++)
		{
			arr2[i]=arr[index[i]];
		}
		return arr2;
	}
	var x = 1; //initlal text box count
	$(add_visc).click(function(e){ //on add input button click
	e.preventDefault();
	add_field(visc_holder,'visc',x++,max_level,'Viscosity at 135&#8451');
	});
	
	$("#pg-calc").click(function(){
	$('#result').hide('fast');
	$('#result-error').hide('fast');
	var errormsg = 'Error List:<br/><ol>';
	var error = 0;
	var fpt = ($("#flp").val());
	var visc = ($("#visc1").val());
	var ds1temps = ($("#temp-ds1-1").val().split(',')); //Comma seperated valued to string arrays
	var ds1vals = ($("#ds1-1").val().split(',')); //Comma seperated valued to string arrays
	var ds2temps = ($("#temp-ds2-1").val().split(',')); //Comma seperated valued to string arrays
	var ds2vals = ($("#ds2-1").val().split(',')); //Comma seperated valued to string arrays
	var ds3temps = ($("#temp-ds3-1").val().split(',')); //Comma seperated valued to string arrays
	var ds3vals = ($("#ds3-1").val().split(',')); //Comma seperated valued to string arrays
	var cstemps = ($("#temp-cs-1").val().split(',')); //Comma seperated valued to string arrays
	var s_vals = ($("#s-1").val().split(',')); //Comma seperated valued to string arrays
	var m_vals = ($("#m-1").val().split(',')); //Comma seperated valued to string arrays
	var dttemps = ($("#temp-dt-1").val().split(',')); //Comma seperated valued to string arrays
	var dtvals = ($("#dt-1").val().split(',')); //Comma seperated valued to string arrays
	fpt = arrstr2posnum (parseFloat(fpt));
	visc = arrstr2posnum (parseFloat(visc));
	ds1temps = arrstr2posnum (ds1temps);
	ds1vals = arrstr2posnum (ds1vals);
	ds2temps = arrstr2posnum (ds2temps);
	ds2vals = arrstr2posnum (ds2vals);
	ds3temps = arrstr2posnum (ds3temps);
	ds3vals = arrstr2posnum (ds3vals);
	cstemps = arrstr2negnum (cstemps);
	s_vals = arrstr2posnum (s_vals);
	m_vals = arrstr2posnum (m_vals);
	dttemps = arrstr2negnum (dttemps);
	dtvals = arrstr2posnum (dtvals);

	//Validating Flash point input
	errormsg +=validateSingleValue(fpt,'Flash Point Temperature');
	//Validating Viscosity input
	errormsg +=validateSingleValue(visc,'Viscosity');

	//Validating DS1 temperature and test values
	errormsg +=dsrvalidate(ds1temps,ds1vals,"Original Binder");
	//Validating DS2 temperature and test values
	errormsg +=dsrvalidate(ds2temps,ds2vals,"RTFO");
	//Validating DS3 temperature and test values
	errormsg +=dsrvalidate(ds3temps,ds3vals,"PAV");
	if (cstemps==0)
		errormsg += '<li>Creep Stiffness test temperatures should be numeric</li>';
	if (cstemps==-1)
		errormsg += '<li>Creep Stiffness test temperatures should be Negative</li>';
	if (s_vals==0)
		errormsg += '<li>S-values should be numeric</li>';
	if (s_vals==-1)
		errormsg += '<li>m-values should be positive</li>';
	if (m_vals==0)
		errormsg += '<li>m-values should be numeric</li>';
	if (m_vals==-1)
		errormsg += '<li>S-values should be positive</li>';
	if(s_vals!=0 && s_vals !=-1 && cstemps!=0 && cstemps !=-1 && m_vals!=0 && m_vals !=-1 && typeof m_vals.length !== "undefined" && typeof s_vals.length !== "undefined")
	{
		if (s_vals.length != cstemps.length ||  cstemps.length!= m_vals.length)
			errormsg += '<li>Number of Creep stiffness test temperatures and number of s and m values should be same\nYou provided '+cstemps.length+' temperatures and '+s_vals.length+' s-values and '+m_vals.length+' m-values</li>';
	}
	if (dttemps==0)
		errormsg += '<li>Direct Tension test temperatures should be numeric</li>';
	if (dttemps==-1)
		errormsg += '<li>Direct Tension test temperatures should be Negative</li>';
	if (dtvals==0)
		errormsg += '<li>Direct Tension test results should be numeric</li>';
	if (s_vals==-1)
		errormsg += '<li>Direct tension test values should be positive</li>';
	if(dtvals!=0 && dttemps!=0 && typeof dtvals.length !== "undefined" && typeof dttemps.length !== "undefined")
	{
		if (dtvals.length != dttemps.length)
			errormsg += 'Number of Direct Tension test temperatures and number of Direct Tension results should be same.\nYou provided '+dttemps.length+' temperatures and '+dtvals.length+' Direct tension values\n';
	}
	if (errormsg !='Error List:<br/><ol>')
	{
		$('#result-error').html('<h4>'+errormsg+'</ol></h4>');
		$('#result-error').show('fast');	
	}
	
	else{
		var errormsgfrmcalc='';
		if (fpt<230){
			errormsgfrmcalc += "Failed in Flash point temperature";
			$('#result-error').html('<h4>'+errormsgfrmcalc+'</h4>');
			$('#result-error').show('fast');
		}
		else if (visc>3){
			errormsgfrmcalc += "Failed in viscosity at 135&#8451";
			$('#result-error').html('<h4>'+errormsgfrmcalc+'</h4>');
			$('#result-error').show('fast');
		}
		else{
			var ds1_indices=checkVals(ds1vals,'gt',1);//get the indices where ds1vals is greater than 1
			var ds2_indices=checkVals(ds2vals,'gt',2.2);//get the indices where ds2vals is greater than 2.2
			var ds3_indices=checkVals(ds3vals,'lt',5000);//get the indices where ds3vals is less than 5000
			var s_indices=sm_indices(s_vals,m_vals,cstemps, dttemps, dtvals);//get the indices where s_vals is less than 300
			var valid_ds1_temps = index2val(ds1_indices,ds1temps);
			var valid_ds2_temps = index2val(ds2_indices,ds2temps);
			var valid_ds3_temps = index2val(ds3_indices,ds3temps);
			var valid_s_temps = index2val(s_indices,cstemps);
			
			var highTemps = _.intersection(valid_ds1_temps,valid_ds2_temps); //function from underscore.js
			var pg_intermediate_temps = pg_inter_temp(highTemps,valid_s_temps);
			var valid_intermediate_temps = _.intersection(pg_intermediate_temps,valid_ds3_temps); //function from underscore.js
			var ult_high = _.max(highTemps);
			var ult_inter = _.max(valid_intermediate_temps);
			var ult_low = ult_high - (ult_inter-4)*2;
			if(!ult_high){
				errormsgfrmcalc += 'Unsuccessful in determining PG high temp';
				$('#result-error').html('<h4>'+errormsgfrmcalc+'</ol></h4>');
				$('#result-error').show('fast');
			}
			else if(!ult_low){
				errormsgfrmcalc += 'Unsuccessful in determining PG low temp\n';
				$('#result-error').html('<h4>'+errormsgfrmcalc+'</h4>');
				$('#result-error').show('fast');
			}
			else {
			$('#result').html('<h4>Provided data matches best with <span class="result">PG '+ult_high+'-'+ult_low+'</span></h4>');
			$('#result').show('fast');
			}
		}
	}

	});
	
	});
	
	
