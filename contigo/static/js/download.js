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

function displayFasta(contig_name, consensus){
    var fastaSeq = '';
    var i = 0;
    var cend = 0
	for(i = 0; i < consensus.length; i += 50) {
	    cend = i + 50;
	    if (cend >= consensus.length) {
		cend = consensus.length;
	    }
	    fastaSeq = fastaSeq + consensus.substring(i, cend) + '\n';
	}
    
    return "# 50 chars per line\n>" + contig_name + " length=" + consensus.length + "\n" + fastaSeq;
}

function displayQuality(contig_name, quality){
    var fastaSeq = '';
    var i = 0;
    var cend = 0;
    var width = 25;
    for(i = 0; i < quality.length; i += width){
	cend = i + width;
	if (cend >= quality.length) {
	    cend = quality.length;
	}
	fastaSeq = fastaSeq + quality.slice(i, cend).join(' ') + '\n';
    }
    return "# 25 values per line\n>" + contig_name + " length=" + quality.length + "\n" + fastaSeq;
}

function loadFasta(contig_name) {
    var seq = '';
    $.ajax(
	   {
	       url: 'json/' + contig_name + ".json",
		   dataType: "json",
		   async: false,
		   success: function(data)
		   {
		       seq = displayFasta(contig_name, data.consensus);
		       openFastaWindow(contig_name, seq);
		   } 
	   });
}

function loadQuality(contig_name) {
    var seq = '';
    $.ajax(
	   {
	       url: 'json/' + contig_name + ".json",
		   dataType: "json",
		   async: false,
		   success: function(data)
		   {
		       seq = displayQuality(contig_name, data.quality);
		       openFastaWindow(contig_name, seq);
		   } 
	   });
}

function downloadFasta(event) {
    var i = 0;
    var txt = '';
    var callback1 = function(){
	alert("still loading...");
    };
    top.consoleRef=window.open('','fastaWindow',
			       'width=700,height=700'
			       +',menubar=0'
			       +',toolbar=0'
			       +',status=0'
			       +',scrollbars=1'
			       +',resizable=1');
    top.consoleRef.document.writeln(
				    '<html><head><title>Assembly Contigs</title></head>'
				    +'<body onLoad="self.focus()" style="background-color: #FFF; color: #000;"><pre style="font-family: Courier New, fixed; font-size: 0.8em;">'
				    );
    for (var i = 0; i < CONTIGO.assembly.length; i += 1) {
	var contig_name = CONTIGO.assembly[i].name;
	var my_data = '';
	$.ajax({
		url: 'json/' + contig_name + ".json",
		    dataType: "json",
		    async: false,
		    success: function(data)
		    {
			top.consoleRef.document.writeln(displayFasta(contig_name, data.consensus));
		    }
	    });
    }
    
    top.consoleRef.document.writeln('</pre></body></html>');
    
}

function downloadQuality(event) {
    var i = 0;
    var txt = '';
    top.consoleRef=window.open('','fastaWindow',
			       'width=700,height=700'
			       +',menubar=0'
			       +',toolbar=0'
			       +',status=0'
			       +',scrollbars=1'
			       +',resizable=1');
    top.consoleRef.document.writeln(
				    '<html><head><title>Assembly Contigs</title></head>'
				    +'<body onLoad="self.focus()" style="background-color: #FFF; color: #000;"><pre style="font-family: Courier New, fixed; font-size: 0.8em;">'
				    );
    for (var i = 0; i < CONTIGO.assembly.length; i += 1) {
	var contig_name = CONTIGO.assembly[i].name;
	var my_data = '';
	$.ajax(
	       {
		   url: 'json/' + contig_name + ".json",
		       dataType: "json",
		       async: false,
		       success: function(data)
		       {
			   top.consoleRef.document.writeln(displayQuality(contig_name, data.quality));
		       }
	       });
    }
    
    top.consoleRef.document.writeln('</pre></body></html>');
    
}

function downloadContigTable(event) {
    var i = 0;
    top.consoleRef=window.open('','fastaWindow',
			       'width=700,height=700'
			       +',menubar=0'
			       +',toolbar=0'
			       +',status=0'
			       +',scrollbars=1'
			       +',resizable=1');
    top.consoleRef.document.writeln(
				    '<html><head><title>Assembly Contigs</title></head>'
				    +'<body onLoad="self.focus()" style="background-color: #FFF; color: #000;"><pre style="font-family: Courier New, fixed; font-size: 0.8em;">Name\tLength\tAvgDepth\t%LowQual\t%GC\tAvgPairedDepth\tAvgInsSize'
				    );
    for (var i = 0; i < CONTIGO.assembly.length; i += 1) {
	var contig = CONTIGO.assembly[i];
	
	top.consoleRef.document.writeln([contig.name, contig.length, contig.avg_depth, contig.perc_low_qual, contig.avg_gc, contig.avg_paired_depth, contig.avg_ins_size].join("\t"));
    }
    top.consoleRef.document.writeln('</pre></body></html>');
}


function loadReads(contig_name, num_segments) {
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
    for (var i = 0; i <= (num_segments - 1); i += 1) {
	var filename = 'json/' + contig_name + '_' + i + '_reads.json';
	$.ajax({ url: filename,
		    dataType: "json",
		    async: false,
		    success: function(data)
		    {
			
			var reads = new TAFFY(data);
			var res = "";
			reads.forEach(function (f,n) {
				if (!visited[f.name]) {
				    res += '>' + f.name + "\n" + f.sequence + "\n";
				    count += 1;
				}
				visited[f.name] = true;
				
			    });
			top.consoleRef.document.writeln(res);
		    }
	    });
    }
    top.consoleRef.document.writeln('*** Fetched ' + count + ' reads ***</pre></body></html>');
}


