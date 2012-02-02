CONTIGO.float_precision = 4;

function xy_chart(x_cat, y_cat) {
    /*
      Plots XY-Chart: 
          - Allowed axis are specified in pulldown menus of HTML
	  - Ranges are REQUIRED and are specified in assembly.js
	  - Data points are specified in assembly.js
	  - Axis transforms and tick_labels are OPTIONALLY specified in assembly.js.
     */
    
    var d1 = [];
    var i = 0;
    
    // Get ranges, REQUIRED!
    var min_x = CONTIGO.ranges.min[x_cat];
    var max_x = CONTIGO.ranges.max[x_cat];
    var min_y = CONTIGO.ranges.min[y_cat];
    var max_y = CONTIGO.ranges.max[y_cat];
    
    $("#xy_chart_detail").html("<span class=\"help\">Click point to get detailed info.</span>");
    
    // Optionally, transform data range to log
    if (CONTIGO.is_log[x_cat]) {
       min_x = Math.log(min_x)/Math.LN10;
       max_x = Math.log(max_x)/Math.LN10;
    }
    if (CONTIGO.is_log[y_cat]) {
       min_y = Math.log(min_y)/Math.LN10;
       max_y = Math.log(max_y)/Math.LN10;
    }
    
    // Round off data ranges
    min_x = Math.floor(min_x);
    min_y = Math.floor(min_y);
    max_x = Math.ceil(max_x);
    max_y = Math.ceil(max_y);
    
    // Set x tick labels
    var x_axis = { min: min_x, max: max_x, ticks: null };
    
    // If log x or y axis, push tick custom labels
    var x_ticks = [];
    if (CONTIGO.is_log[x_cat]){
	for (i = min_x; i <= max_x; i += 1){
            x_ticks.push([i,CONTIGO.tick_labels[x_cat][i]]);
	}
	x_axis.ticks = x_ticks;
    }
    
    var y_axis = { min: min_y, max: max_y, ticks: null };
    var y_ticks = [];
    if (CONTIGO.is_log[y_cat]){
	for (i = min_y; i <= max_y; i += 1){
            y_ticks.push([i,CONTIGO.tick_labels[y_cat][i]]);
	}
	y_axis.ticks = y_ticks;
    }    
    
    
    // Get data points
    var x_val = null;
    var y_val = null;
    for (i = 0; i < CONTIGO.assembly.length; i += 1) {
	x_val = CONTIGO.assembly[i][x_cat];
        y_val = CONTIGO.assembly[i][y_cat];
        if (CONTIGO.is_log[x_cat]) {
            x_val = Math.log(CONTIGO.assembly[i][x_cat])/Math.LN10;
            
        }
        
        if (CONTIGO.is_log[y_cat]) {
            y_val = Math.log(CONTIGO.assembly[i][y_cat])/Math.LN10;
            
        }
        d1.push([x_val, y_val]);
    }
    
    // plot options
    var options = { xaxis: x_axis,                         
                    yaxis: y_axis,
                    colors: ["#000"],
                    selection: {
                                 mode: "xy",
                                 color: "#f00"
	                       },
                    grid: {
                           color: "#000",
                           backgroundColor: "#fff",
                           borderWidth: 1,
                           tickColor: "#ccc",
                           clickable: true,
			   autoHighlight: true
		    }
                  
		    
    };
    $.plot($("#xy_chart"), [{ data: d1, 
                              points: { show: true, radius: 0.5},
                             shadowSize: 0,
                             colors: ["#000"]
                           }                              
                          ],
                          options
      );
  
}

$("#xy_chart").bind("plotclick", function (event, pos, item) {
    var x_cat = $('#xy_chart_x').val();
    var y_cat = $('#xy_chart_y').val();
    var x_txt = $('#xy_chart_x :selected').text();
    var y_txt = $('#xy_chart_y :selected').text();
    if (item) {
	x_val = item.datapoint[0];
	y_val = item.datapoint[1];
	if (CONTIGO.is_log[x_cat]){
	    x_val = Math.pow(10, x_val);
	} 
	if (CONTIGO.is_log[y_cat]){
	    y_val = Math.pow(10, y_val);
	} 
	x_val = Math.round(x_val);
	y_val = Math.round(y_val);
	$("#xy_chart_detail").html(CONTIGO.assembly[item.dataIndex].name + ": " + CONTIGO.assembly[item.dataIndex].length + "nt, " + CONTIGO.assembly[item.dataIndex].avg_depth + "x, " + CONTIGO.assembly[item.dataIndex].avg_gc + "%GC");
	
    }
});


$("#xy_chart").bind("plotselected", function(event, ranges) {
    var x_cat = $('#xy_chart_x').val();
    var y_cat = $('#xy_chart_y').val();
    var x_txt = $('#xy_chart_x :selected').text();
    var y_txt = $('#xy_chart_y :selected').text();
    
    var x1_val = ranges.xaxis.from;
    var x2_val = ranges.xaxis.to;
    var y1_val = ranges.yaxis.from;
    var y2_val = ranges.yaxis.to;
    
    if (CONTIGO.is_log[x_cat]){
	x1_val = Math.pow(10, x1_val);
	x2_val = Math.pow(10, x2_val);
    } 
    if (CONTIGO.is_log[y_cat]){
	y1_val = Math.pow(10, y1_val);
	y2_val = Math.pow(10, y2_val);
    } 
    x1_val = Math.round(x1_val);
    y1_val = Math.round(y1_val);
    x2_val = Math.round(x2_val);
    y2_val = Math.round(y2_val);
    
    alert(x_txt + ": " + x1_val + "-" + x2_val + " (" + (x2_val-x1_val) + "), " + y_txt + ": " + y1_val + "-" + y2_val + " (" + (y2_val-y1_val) + ")");
});

function m_chart(cat){
       
    var d = [];
    var min = CONTIGO.ranges.min[cat];
    var max = CONTIGO.ranges.max[cat];
    var i = 0;
    $("#m_chart_detail").html("<span class=\"help\">Click point to get approx value</span>");
    if (CONTIGO.is_log[cat]) {
       min = Math.log(min)/Math.LN10;
       max = Math.log(max)/Math.LN10;
    }
    min = Math.floor(min);
    max = Math.ceil(max);
    var x_axis = { min: min,
		   max: max,
		   ticks: null
             };
    var ticks = [];
    if (CONTIGO.is_log[cat]){
       for (i = min; i <= max; i += 1){
           ticks.push([i,CONTIGO.tick_labels[cat][i]]);
       }
       x_axis.ticks = ticks;
    }

    var total_length = 0;
    for (i = 0; i < CONTIGO.assembly.length; i += 1) {
       d.push([CONTIGO.assembly[i].length, CONTIGO.assembly[i][cat]]);
       total_length += CONTIGO.assembly[i].length;
    }
    d.sort(function(a,b){return b[1] - a[1];});
    var sum = 0;
    var d2 = [];
    var x_val = null;
    var y_val = null;
    var x_95 = 0;
    for (i = 0; i < (d.length); i += 1) {
        sum += d[i][0];
        x_val = d[i][1];
        y_val = (sum/total_length)*100;
        if (CONTIGO.is_log[cat]){
           x_val = Math.log(d[i][1])/Math.LN10;
        }
        if (y_val < 95){
	    x_95 = x_val;
	}
        d2.push([x_val, y_val]);
        
    }
    var options = { xaxis: x_axis,
                    colors: ["#000"],
                    selection: {
                                 mode: "xy",
                                 color: "#f00"
	                       },
                    grid: {
                           color: "#000",
                           backgroundColor: "#fff",
                           borderWidth: 1,
                           tickColor: "#ccc",
                           clickable: true,
                           autoHighlight: true
			   
		    }
		    

    };
    $.plot($("#m_chart"), [{ data: d2, 
                             points: { show: true, radius: 0.5 },
                             shadowSize: 0,
                             colors: ["#000"],
                             fillColor: ["#000"]
                            }                              
                          ],
                          options
                          
      );
   
}

$("#m_chart").bind("plotclick", function (event, pos, item) {
    var cat = $('#m_chart_cat').val();
    var txt = $('#m_chart_cat :selected').text();
    if (item) {
	x_val = item.datapoint[0];
	y_val = item.datapoint[1];
	if (CONTIGO.is_log[cat]){
	    x_val = Math.pow(10, x_val);
	}
	x_val = Math.round(x_val);
	y_val = y_val.toFixed(2);
	$("#m_chart_detail").html(y_val + "% with " + txt + " >= " + x_val);
    }
});

$("#m_chart").bind("plotselected", function(event, ranges) {
    var x_cat = $('#m_chart_cat').val();
    var x_txt = $('#m_chart_cat :selected').text();
    var y_txt = "Percent";
    
    var x1_val = ranges.xaxis.from;
    var x2_val = ranges.xaxis.to;
    var y1_val = ranges.yaxis.from;
    var y2_val = ranges.yaxis.to;
    
    if (CONTIGO.is_log[x_cat]){
	x1_val = Math.pow(10, x1_val);
	x2_val = Math.pow(10, x2_val);
    } 
    x1_val = Math.round(x1_val);
    y1_val = Math.round(y1_val);
    x2_val = Math.round(x2_val);
    y2_val = Math.round(y2_val);
    
    alert(x_txt + ": " + x1_val + "-" + x2_val + " (" + (x2_val-x1_val) + "), " + y_txt + ": " + y1_val + "-" + y2_val + " (" + (y2_val-y1_val) + ")");
    
});
    

function h_chart(cat) {
/*
  Plots Histogram 
*/

    var d3 = [];
    var sum = 0;
    var min = CONTIGO.ranges.min[cat];
    var max = CONTIGO.ranges.max[cat];
    var i = 0;
    $("#h_chart_detail").html("<span class=\"help\">Click point to get approx value</span>");
    if (CONTIGO.is_log[cat]) {
	min = Math.log(min)/Math.LN10;
	max = Math.log(max)/Math.LN10;
    }
    min = Math.floor(min);
    max = Math.ceil(max);
    var x_axis = { min: min,
		   max: max,
		   ticks: null
             };
    var ticks = [];
    if (CONTIGO.is_log[cat]){
       for (i = min; i <= max; i += 1){
           ticks.push([i,CONTIGO.tick_labels[cat][i]]);
       }
       x_axis.ticks = ticks;
    }
    for (i = 0; i < CONTIGO.hist_data[cat].length; i += 1) {
	sum += CONTIGO.hist_data[cat][i];
    }
    var y_min = 1;
    var y_max = 0;
    for (i = 0; i < CONTIGO.hist_data[cat].length; i += 1) {
	y_val = CONTIGO.hist_data[cat][i];
	if (y_val > 0){
	    
	    x_val = i;
	    if (CONTIGO.is_log[cat]){
		x_val = Math.log(x_val)/Math.LN10;
	    }
	    y_val = Math.log(y_val)/Math.LN10;
	    if (y_val == -Infinity) {
		y_val = 0;
	    }
	    if (x_val == -Infinity) {
		x_val = 0;
	    }
	    if (y_val > y_max) {
		y_max = y_val;
	    }
	    if (y_val < y_min) {
		y_min = y_val;
	    }
	    d3.push([x_val, y_val]);
	}
	
    }
    y_min = Math.floor(y_min);
    y_max = Math.ceil(y_max);
    var y_axis = { min: y_min,
		   max: y_max,
		   ticks: null
             };
    
    ticks = [];
    for (i = y_min; i <= y_max; i += 1){
        ticks.push([i,Math.pow(10,i)]);
	
    }
    y_axis.ticks = ticks;
    
    var options = { xaxis: x_axis,
		    yaxis: y_axis,
                    colors: ["#000"],
                    selection: {
 	              mode: "xy",
 	              color: "#f00"
	            },
                    grid: {
                        color: "#000",
                        backgroundColor: "#fff",
			borderWidth: 1,
                        tickColor: "#ccc",
                        clickable: true,
			autoHighlight: true
                    }
                  };
    
    $.plot($("#h_chart"), [{ data: d3,
                             points: { show: true, radius: 0.5 },
                             shadowSize: 0,
                             colors: ["#000"],
                             fillColor: ["#000"]
                           }                              
                          ],
	                  options
          );
}

$("#h_chart").bind("plotclick", function (event, pos, item) {
    var cat = $('#h_chart_cat').val();
    var txt = $('#h_chart_cat :selected').text();
    var x_val = null;
    var y_val = null;
    if (item) {
	x_val = item.datapoint[0];
	y_val = item.datapoint[1];
	if (CONTIGO.is_log[cat]){
	    x_val = Math.pow(10, x_val);
	}
	
	x_val = Math.round(x_val);
	y_val = Math.pow(10,y_val).toFixed(0);
	$("#h_chart_detail").html(y_val + " with " + txt + " = " + x_val);
    }
});



$("#h_chart").bind("plotselected", function(event, ranges) {
    var x_cat = $('#h_chart_cat').val();
    var x_txt = $('#h_chart_cat :selected').text();
    var y_txt = "Counts";
    
    var x1_val = ranges.xaxis.from;
    var x2_val = ranges.xaxis.to;
    var y1_val = ranges.yaxis.from;
    var y2_val = ranges.yaxis.to;
    var y_sum = 0;
    var sum = 0;
    var perc = 0.0;
    
    if (CONTIGO.is_log[x_cat]){
	x1_val = Math.pow(10, x1_val);
	x2_val = Math.pow(10, x2_val);
    }
    x1_val = Math.round(x1_val);
    x2_val = Math.round(x2_val);
    y1_val = Math.pow(10,y1_val).toFixed(0);
    y2_val = Math.pow(10,y2_val).toFixed(0);
    for (i = 1; i < CONTIGO.hist_data[x_cat].length; i += 1) {
	sum += CONTIGO.hist_data[x_cat][i];
	    if ((i >= x1_val) && (i <= x2_val) && (CONTIGO.hist_data[x_cat][i] >= y1_val) && (CONTIGO.hist_data[x_cat][i] <= y2_val)){
		y_sum += CONTIGO.hist_data[x_cat][i];
	    }
    }
    
    perc = ((y_sum/sum)*100).toFixed(2);
    alert(y_sum + " (" + perc + "%) w/ " + x_txt + " from " + x1_val + "-" + x2_val + " and " + y_txt + " from " + y1_val + "-" + y2_val);
    
});
