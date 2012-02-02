// var PRECISION = 5;
var viewer = null;
var fastaSeq = null;
var point = null;
var zoom = null;

function init_seadragon() {
    //CONTIGO.current_contig_length = '';
    document.getElementById("zoomtig").innerHTML = "";    // for CMS
    viewer = new Seadragon.Viewer("zoomtig");
    viewer.addEventListener("open", showViewport);
    viewer.addEventListener("animation", showViewport);
    Seadragon.Utils.addEvent(viewer.elmt, "mousemove", showMouse);
    viewer.addControl(makeFetcherControl(), Seadragon.ControlAnchor.TOP_RIGHT);
    viewer.addControl(makeAssemblyControl(), Seadragon.ControlAnchor.TOP_LEFT);
    viewer.addControl(makeReadnameControl(), Seadragon.ControlAnchor.TOP_LEFT);
    viewer.addControl(makeTemplateControl(), Seadragon.ControlAnchor.TOP_LEFT);
    viewer.addEventListener("open", onViewerOpen);
}


function makeFetcherControl() {
    
    var control = document.createElement("a");
    var controlText = document.createTextNode("Text Reads");
    control.href = "contigo.html";
    control.className = "control";
    control.setAttribute('style', 'font-size:10px; color: #000; background-color: white; border: 1px solid #000; padding: 2px; margin: 2px;');
    control.setAttribute('title', 'View FASTA reads contained within the viewable area');
    control.appendChild(controlText);
    Seadragon.Utils.addEvent(control, "click", onControlClick);
    return control;
}

function makeReadnameControl() {
    
    //var container = document.createElement("div");
    //container.setAttribute('style', 'font-size:10px; color: #800; background-color: white;');
    //container.innerHTML += "View: ";
    
    var control = document.createElement("a");
    var controlText = document.createTextNode("Names/Status");
    control.href = "contigo.html";
    control.className = "control";
    control.setAttribute('style', 'font-size:10px; color: #000; background-color: white; border: 1px solid #000; padding: 2px; margin: 2px;');
    control.setAttribute('title', 'Switch to alignment with read names');
	
    control.appendChild(controlText);
    Seadragon.Utils.addEvent(control, "click", onReadnameControlClick);
    
    //container.appendChild(control);
    //container.innerHTML += " ";
    return control;
}

function makeAssemblyControl() {
    
    var control = document.createElement("a");
    var controlText = document.createTextNode("Alignment");
    control.href = "contigo.html";
    control.className = "control";
    control.setAttribute('style', 'font-size:10px; color: #000; background-color: white; border: 1px solid #000; padding: 2px; margin: 2px;');
    control.setAttribute('title', 'Switch to alignment with read sequences');
    control.appendChild(controlText);
    Seadragon.Utils.addEvent(control, "click", onAssemblyControlClick);
    return control;
}
function makeTemplateControl() {
    
    var control = document.createElement("a");
    var controlText = document.createTextNode("Templates");
    control.href = "contigo.html";
    control.className = "control";
    control.setAttribute('style', 'font-size:10px; color: #000; background-color: white; border: 1px solid #000; padding: 2px; margin: 2px;');
    control.setAttribute('title', 'Switch to alignment of templates');
    control.appendChild(controlText);
    Seadragon.Utils.addEvent(control, "click", onTemplateControlClick);
    return control;
}


function openFastaWindow(title, sequence) {
    top.consoleRef=window.open('','fastaWindow',
			       'width=700,height=700'
			       +',menubar=0'
			       +',toolbar=0'
			       +',status=0'
			       +',scrollbars=1'
			       +',resizable=1')
	top.consoleRef.document.writeln(
					'<html><head><title>'+ title +'</title></head>'
					+'<body onLoad="self.focus()" style="background-color: #FFF; color: #000;"><pre style="font-family: Courier New, fixed; font-size: 0.8em;">'
					+sequence
					+'</pre></body></html>')
	
	// top.consoleRef.close();
	
	}

function fetchReads(event) {
    var pixel = Seadragon.Utils.getMousePosition(event).minus(
							      Seadragon.Utils.getElementPosition(viewer.elmt));
    var win = viewer.viewport.getBounds();
    var span = (win.width * CONTIGO.current_contig_assembly_length);
    var ws = (((win.x * CONTIGO.current_contig_assembly_length) / 1) - (Math.abs(CONTIGO.current_contig_padded_start) / 1));
    var we = ws + span;
    ws = ws.toFixed(0);
    we = we.toFixed(0);
    var ss = parseInt((ws - CONTIGO.current_contig_padded_start) / CONTIGO.default_segment_width);
    var se = parseInt((we - CONTIGO.current_contig_padded_start) / CONTIGO.default_segment_width);
    var visited = {};
    top.consoleRef=window.open('','fastaWindow',
			       'width=700,height=700'
			       +',menubar=0'
			       +',toolbar=0'
			       +',status=0'
			       +',scrollbars=1'
			       +',resizable=1');
    top.consoleRef.document.writeln('<html><head><title>Reads</title></head>'
				    +'<body onLoad="self.focus()" style="background-color: #FFF; color: #000;"><pre style="font-family: Courier New, fixed; font-size: 0.8em;">');
    var count = 0;
    for (var i = ss; i <= se; i += 1) {
	var filename = 'json/' + CONTIGO.current_contig_name + '_' + i + '_reads.json';
	$.ajax({ url: filename,
		    dataType: "json",
		    async: false,
		    success: function(data)
		    {
			
			var reads = new TAFFY(data);
			var res = "";
			reads.forEach(function (f,n) {
				var rs = parseInt(f.padded_start);
				var re = rs + parseInt(f.padded_length);
				if ((rs >= ws) && (re <= we)){
				    if (!visited[f.name]) {
					res += '>' + f.name + "\n" + f.sequence + "\n";
					count += 1;
				    }
				    visited[f.name] = true;
				}
			    });
			top.consoleRef.document.writeln(res);
		    }
	    });
    }
    top.consoleRef.document.writeln('*** Fetched ' + count + ' reads ***</pre></body></html>');
}

function onControlClick(event) {
    Seadragon.Utils.cancelEvent(event);    // so link isn't processed
    if (!viewer.isOpen()) {
	return;
    }
    fetchReads(event);
}

function onReadnameControlClick(event) {
    var dzi = 'dzi_2\/' + CONTIGO.current_contig_name + '.dzi';
    if (viewer.isOpen()) {
	point = viewer.viewport.getCenter();
	zoom = viewer.viewport.getZoom();
    }
    viewer.openDzi(dzi);
    Seadragon.Utils.cancelEvent(event);
}

function onTemplateControlClick(event) {
    var dzi = 'dzi_3\/' + CONTIGO.current_contig_name + '.dzi';
    if (viewer.isOpen()) {
	point = viewer.viewport.getCenter();
	zoom = viewer.viewport.getZoom();
    }
    viewer.openDzi(dzi);
    Seadragon.Utils.cancelEvent(event);
}

function onAssemblyControlClick(event) {
    var dzi = 'dzi\/' + CONTIGO.current_contig_name + '.dzi';
    if (viewer.isOpen()) {
	point = viewer.viewport.getCenter();
	zoom = viewer.viewport.getZoom();
    }
    viewer.openDzi(dzi);
    Seadragon.Utils.cancelEvent(event);
}

function switchTo(event, dzi, contig_name, contig_length, contig_padded_length, padded_start, padded_end, max_depth, height, ruler_height) {
    Seadragon.Utils.cancelEvent(event);
    if (dzi) {
	CONTIGO.current_contig_name = contig_name;
	CONTIGO.current_contig_assembly_length = padded_end - padded_start;
	CONTIGO.current_contig_padded_start = padded_start;
	CONTIGO.current_contig_padded_end = padded_end;
	CONTIGO.current_contig_length = contig_length;
	CONTIGO.current_contig_padded_length = contig_padded_length;
	CONTIGO.current_contig_max_depth = max_depth;
	CONTIGO.current_contig_height = height;
	CONTIGO.current_contig_ruler_height = ruler_height;
	$('#contig_viewed').html("<b>" + contig_name +                                      "</b> [" + CONTIGO.current_contig_assembly_length +                                      " cols, " + CONTIGO.current_contig_length +                                      " bp, " + (CONTIGO.current_contig_padded_length - CONTIGO.current_contig_length ) + " pads]");
	viewer.openDzi(dzi);
    } else {
	viewer.close();
	alert("Error with contig viewer. Please contact the developer.");
    }
    
}

function showViewport(viewer) {
    if (!viewer.isOpen()) {
	return;
    }
    
    var sizePoints = viewer.viewport.getBounds().getSize();
    var window_size = ((sizePoints.x * CONTIGO.current_contig_assembly_length)/1).toFixed(0);
    document.getElementById("contig_window_size").innerHTML = window_size;
}

function showMouse(event) {
    
    var pixel = Seadragon.Utils.getMousePosition(event).minus(
							      Seadragon.Utils.getElementPosition(viewer.elmt));
    if (!viewer.isOpen()) {
	return;
    }
    
    var point = viewer.viewport.pointFromPixel(pixel);
    var bp = (((point.x * CONTIGO.current_contig_assembly_length) / 1) - (Math.abs(CONTIGO.current_contig_padded_start) / 1)).toFixed(0);
    
    var depth = 0;
    if (point.y > CONTIGO.current_contig_ruler_height){
	var y1 = point.y - CONTIGO.current_contig_ruler_height;
	var y2 = CONTIGO.current_contig_height - CONTIGO.current_contig_ruler_height;
	
	depth = CONTIGO.current_contig_max_depth * (y1/y2);
    }
    document.getElementById("contig_pos").innerHTML = bp;
    document.getElementById("contig_depth").innerHTML = Math.ceil(depth);
}

function onViewerOpen(viewer) {
    if (point && zoom) {
	viewer.viewport.panTo(point, true);
	viewer.viewport.zoomTo(zoom, true);
    }
}

Seadragon.Utils.addEvent(window, "load", init_seadragon);
Seadragon.Config.immediateRender  = true;