const CLI = await new Aioli([
    "samtools/1.10", 
    "minimap2/2.22",
    // {
    //     tool: "faster2",
    //     version: "0.3.0",
    //     urlPrefix: "http://localhost:8000/faster2"
    // }
], {debug: true});

const MAX_POINTS = 300; 
const SPARKLINE_COLOR = '#4c9aff';

// **NEW GLOBAL VARIABLE** to store max length for proportional scaling
let maxContigLength = 1; // Initialize to 1 to avoid division by zero

// Function to reload the page
function resetPage() {
    window.location.reload();
}

// Aggressive cleaning function to ensure keys match
function cleanContigName(name) {
    // Trim whitespace
    return name.trim();
}

function subsampleDepthData(data) {
    if (data.length <= MAX_POINTS) {
        return data;
    }

    const step = Math.ceil(data.length / MAX_POINTS);
    const subsampled = [];
    
    for (let i = 0; i < data.length; i += step) {
        subsampled.push(data[i]);
    }
    
    if (subsampled.length === 0 && data.length > 0) {
         subsampled.push(data[0]);
    }
    
    return subsampled;
}

// **MODIFIED RENDER FUNCTION** to use proportional width
function renderSparklines() {
    $('#depthList .sparkline-container').each(function() {
        const $this = $(this);
        const depths = $this.data('depths'); 
        const contigLength = $this.data('length'); // Get length
        
        // Calculate proportional width (max 100%)
        const proportionalWidth = Math.max(20, (contigLength / maxContigLength) * 100); 
        
        // Apply proportional width to the container
        $this.css('width', `${proportionalWidth}%`);

        if (depths && depths.length > 0) {
            $this.sparkline(depths, {
                type: 'line',
                width: '100%', // Sparkline fills its container's proportional width
                height: '50px',
                lineColor: SPARKLINE_COLOR, 
                fillColor: SPARKLINE_COLOR + '33',
                spotColor: false,
                minSpotColor: false,
                maxSpotColor: false, 
                chartRangeMin: 0,
                highlightLineColor: '#6c757d', 
                highlightSpotColor: '#6c757d', 
                tooltipFormat: 'Depth: {{y}}' 
            });
        } else {
            $this.text('No depth data');
        }
    });
}

async function run() {
    const outputDiv = document.getElementById("output");
    const ulList = document.getElementById('depthList');
    const downloadButton = document.getElementById("download-bam"); 
    const initialMessage = document.getElementById("initial-message"); 
    
    // 1. Initial State Setup
    // Display 'Running' message and hide the initial instruction
    outputDiv.innerHTML = '<div class="alert alert-info" role="alert">Running analysis... This may take a moment.</div>';
    document.getElementById("btn").disabled = true;
    //downloadButton.style.display = 'none'; // Hide button at start of run
    ulList.innerHTML = ''; 
    
    // Reset maxContigLength at the start of a run
    maxContigLength = 1;

    try {
        const q = document.getElementById("query").files;
        const r = document.getElementById("ref").files;

        if (q.length === 0 || r.length === 0) {
            outputDiv.innerHTML = '<div class="alert alert-danger" role="alert">Please select both a **Reference Genome** and **Query Reads** file(s).</div>';
            document.getElementById("btn").disabled = false;
            return;
        }

        const query = await CLI.mount(q); 
        const ref = await CLI.mount(r);   

        // --- PIPELINE EXECUTION ---
        
        
        const mmap_out = await CLI.exec("minimap2", [
            "-a",
            "-I", "1G",
            "-K", "100M",
            "-o", "output.sam",
            "-x", "map-ont",
            ref[0], 
            ...query
        ]);
        await CLI.exec("samtools view -S -b output.sam -o output.bam");
        
        await CLI.exec("samtools sort -o output.sorted.bam output.bam");
        await CLI.exec("samtools index output.sorted.bam");

        // **NEW STEP 1: Run samtools flagstat**
        const flagstatOutput = await CLI.exec("samtools flagstat output.sorted.bam");

       // --- DOWNLOAD FIX ---
        
        // 1. Read the file content from the Aioli filesystem as a Uint8Array
        const bamData = await CLI.fs.readFile("output.sorted.bam");

        // 2. Create a Blob object from the data
        const bamBlob = new Blob([bamData], { type: 'application/octet-stream' }); 

        // 3. Create a stable, revocable URL for the Blob
        // This URL remains valid as long as the Blob object exists
        const download_url = URL.createObjectURL(bamBlob);

        // **SUCCESS**: Set the download URL and show the button
        downloadButton.href = download_url;
        downloadButton.style.display = 'inline-block'; // Show the button
        
        // Clear the status/output div completely after a successful BAM generation
        outputDiv.innerHTML = '';

        
        // 1. Total Reads: Matches the first line (Total QC-passed reads)
        const totalReadsMatch = flagstatOutput.match(/^(\d+)\s*\+\s*\d+\s*in total/m);
        // 2. Secondary Reads
        const secondaryReadsMatch = flagstatOutput.match(/^(\d+)\s*\+\s*\d+\s*secondary/m);
        // 3. Supplementary Reads
        const supplementaryReadsMatch = flagstatOutput.match(/^(\d+)\s*\+\s*\d+\s*supplementary/m);
        // 4. Mapped Reads: To get the total mapped count and percentage
        const mappedReadsMatch = flagstatOutput.match(/^(\d+)\s*\+\s*\d+\s*mapped\s+\((\d+\.\d+)%/m);
        
        // --- DATA EXTRACTION & CALCULATION ---
        let primaryReads = 'N/A';
        let primaryMapped = 'N/A';
        let primaryMappedPercent = 'N/A';

        // Get the required raw counts (QC-passed)
        const total = totalReadsMatch ? parseInt(totalReadsMatch[1], 10) : 0;
        const secondary = secondaryReadsMatch ? parseInt(secondaryReadsMatch[1], 10) : 0;
        const supplementary = supplementaryReadsMatch ? parseInt(supplementaryReadsMatch[1], 10) : 0;
        
        const mappedTotal = mappedReadsMatch ? parseInt(mappedReadsMatch[1], 10) : 0;
        
        // 1. Calculate Primary Reads
        let calculatedPrimaryReads = total - secondary - supplementary;
        
        if (calculatedPrimaryReads > 0) {
            primaryReads = calculatedPrimaryReads.toLocaleString();
            
            const mappedRatio = total > 0 ? (mappedTotal / total) : 0;
            const calculatedPrimaryMapped = Math.round(calculatedPrimaryReads * mappedRatio);
            
            primaryMapped = calculatedPrimaryMapped.toLocaleString();
            
            // 3. Calculate Primary Mapped Percentage
            const calculatedPrimaryMappedPercent = calculatedPrimaryReads > 0 
                ? ((calculatedPrimaryMapped / calculatedPrimaryReads) * 100).toFixed(2)
                : '0.00';
                
            primaryMappedPercent = calculatedPrimaryMappedPercent;
        }

        const summaryHtml = `
            <div class="alert alert-info summary-box" role="alert">
                <span class="fw-bold">Total reads:</span> ${primaryReads} &nbsp;|&nbsp; 
                <span class="fw-bold">Mapped reads:</span> ${primaryMapped} 
                (${primaryMappedPercent}%)
            </div>
        `;
        document.getElementById("flagstat-summary").innerHTML = summaryHtml;
        
        // ... (rest of the code)


        const depthOutput = await CLI.exec("samtools depth -aa output.sorted.bam");
        
        // Capture samtools coverage output with -H (header)
        const covOutput = await CLI.exec("samtools coverage -H output.sorted.bam");
        
        // --- DATA AGGREGATION & LIST POPULATION ---
        
        // Robust parsing for samtools depth output
        const contigData = {}; 
        const lines = depthOutput.trim().split('\n').filter(line => line.length > 0);
        
        if (lines.length === 0) {
            outputDiv.innerHTML = '<div class="alert alert-warning" role="alert">Analysis completed, but no depth data was generated. Check your input files.</div>';
            downloadButton.style.display = 'none';
            document.getElementById("btn").disabled = false;
            return;
        }

        lines.forEach(line => {
            const fields = line.split('\t');
            if (fields.length === 3) {
                const [rname, pos, depthStr] = fields;
                
                // Use clean function
                const cleanRname = cleanContigName(rname); 
                const depth = parseInt(depthStr.trim(), 10);
                
                if (!contigData[cleanRname]) {
                    // **MODIFIED**: Added 'length' property
                    contigData[cleanRname] = { depths: [], sum: 0, count: 0, length: 0 };
                }
                
                contigData[cleanRname].depths.push(depth);
                contigData[cleanRname].sum += depth;
                contigData[cleanRname].count++;
                contigData[cleanRname].length++; // Increment length (number of depth points)
            }
        });
        
        // **NEW STEP**: Find the maximum contig length for scaling
        for (const rname in contigData) {
            maxContigLength = Math.max(maxContigLength, contigData[rname].length);
        }

        // CORRECTED: Parse samtools coverage output for 9 fields
        const covLines = covOutput.trim().split('\n').filter(line => line.length > 0);
        console.log("Coverage output lines:", covLines);
        const contigCoverageStats = {};

        if (covLines.length >= 1) {
            for (let i = 0; i < covLines.length; i++) {
                // Use split by any whitespace and filter out empty strings
                //const fields = covLines[i].split(/\s+/).filter(f => f.length > 0);
                const fields = covLines[i].split('\t');
                // Check for the expected 9 fields
                if (fields.length >= 8) {
                    try {
                        const rname = fields[0];
                        const cleanCovRname = cleanContigName(rname);
                        console.log("Processing coverage for contig:", cleanCovRname);
                        
                        // Skip the overall summary line ('*')
                        //if (cleanCovRname && cleanCovRname !== '*') {
                            // 9-Field Output Indices:
                            // 0:#rname, 1:startpos, 2:endpos, 3:numreads, 4:coveredbases, 5:coverage, 6:meandepth, 7:meanbaseq, 8:meanmapq
                            const numreads = parseInt(fields[3], 10); 
                            const contLength = parseInt(fields[2], 10);
                            const avgmapq = parseFloat(fields[8]).toFixed(0); 
                            const coverage = parseFloat(fields[5]).toFixed(0);
                            console.log(`Contig: ${cleanCovRname}, Length: ${contLength}, Reads: ${numreads}, Coverage: ${coverage}, Mean MAPQ: ${avgmapq}`);
                            
                            contigCoverageStats[cleanCovRname] = {
                                contLength: isNaN(contLength) ? 0 : contLength,
                                numreads: isNaN(numreads) ? 0 : numreads,
                                coverage: isNaN(parseFloat(coverage)) ? '0.0' : coverage,
                                meanmapq: isNaN(parseFloat(avgmapq)) ? '0' : avgmapq
                            };
                        //}
                    } catch (e) {
                        console.error("Error parsing coverage line:", covLines[i], e);
                    }
                } else {
                    console.warn("Skipping coverage line due to missing fields (expected 9+, got " + fields.length + "):", covLines[i]);
                }
            }
        }
        
        const listItemsHtml = [];
        // Iterate over the successfully parsed depth data (using the cleaned keys)
        for (const rname in contigData) {
            const data = contigData[rname];
            const meanDepth = data.count > 0 ? (data.sum / data.count).toFixed(0) : '0';
            
            const subsampledDepths = subsampleDepthData(data.depths);
            const depthsJSON = JSON.stringify(subsampledDepths);
            
            // **MODIFIED**: Pass the contig length using a data attribute
            const sparklineHtml = `<span class="sparkline-container" data-depths='${depthsJSON}' data-length='${data.length}'></span>`;

            // Lookup coverage stats using the same cleaned key
            const stats = contigCoverageStats[rname] || { numreads: 0, contLength: 0, coverage: '0.0', meanmapq: '0' };
            
            // Format output
            const numreads = stats.numreads.toLocaleString(); 
            const contLength = stats.contLength.toLocaleString();
            const coverage = stats.coverage + '%';
            const meanmapq = stats.meanmapq;

            // Insert new columns into the list item HTML
            const listItem = `
                <li class="depth-item">
                    <div class="contig-name">${rname}</div>
                    <div class="contig-length">${contLength}</div>
                    <div class="mean-depth">${meanDepth}</div>
                    
                    <div class="coverage">${coverage}</div>
                    <div class="mean-mapq">${meanmapq}</div>
                    <div class="depth-profile">${sparklineHtml}</div>
                </li>
            `;
            listItemsHtml.push(listItem);
        }
        
        ulList.innerHTML = listItemsHtml.join('');

        // --- SPARKLINE RENDERING ---
        renderSparklines(); 

    } catch (error) {
        console.error(error);
        // On error, display error message and hide the download button
        downloadButton.style.display = 'none';
        outputDiv.innerHTML = `<div class="alert alert-danger" role="alert">An error occurred during execution: **${error.message || String(error) || 'Unknown Aioli/CLI Error'}**</div>`;
    } finally {
        document.getElementById("btn").disabled = false;
    }
}

// Attach event listeners
document.getElementById("btn").addEventListener("click", run);
document.getElementById("reset-btn").addEventListener("click", resetPage);