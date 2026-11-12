const CLI = await new Aioli([
    "samtools/1.10", 
    "minimap2/2.22",
    "seqtk/1.4"
    // {
    //     tool: "faster2",
    //     version: "0.3.0",
    //     urlPrefix: "http://localhost:8000/faster2"
    // }
], {debug: true});

const NBINS = 300; // number of windows in each contig's depth profile
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

// use proportional width
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

// Progress bar: the wasm tools don't report internal progress, so the bar
// advances per pipeline step and an elapsed timer shows it's still working.
let progressTimer = null;

function showProgress(outputDiv) {
    outputDiv.innerHTML = `
        <div class="alert alert-info" role="alert">
            <div class="d-flex justify-content-between mb-2">
                <span id="progress-label">Starting…</span>
                <span id="progress-elapsed" class="text-muted">0s</span>
            </div>
            <div class="progress" role="progressbar" aria-label="Analysis progress" aria-valuemin="0" aria-valuemax="100">
                <div id="progress-bar" class="progress-bar progress-bar-striped progress-bar-animated" style="width: 0%"></div>
            </div>
        </div>`;
    const t0 = Date.now();
    clearInterval(progressTimer);
    progressTimer = setInterval(() => {
        const el = document.getElementById("progress-elapsed");
        if (el) el.textContent = `${Math.round((Date.now() - t0) / 1000)}s`;
    }, 1000);
}

function setProgress(pct, label) {
    const bar = document.getElementById("progress-bar");
    if (!bar) return;
    bar.style.width = `${pct}%`;
    bar.parentElement.setAttribute("aria-valuenow", pct);
    document.getElementById("progress-label").textContent = label;
}

function stopProgress() {
    clearInterval(progressTimer);
    progressTimer = null;
}

async function run() {
    const outputDiv = document.getElementById("output");
    const ulList = document.getElementById('depthList');
    const downloadButton = document.getElementById("download-bam");

    // 1. Initial State Setup
    showProgress(outputDiv);
    document.getElementById("btn").disabled = true;
    //downloadButton.style.display = 'none'; // Hide button at start of run
    ulList.innerHTML = ''; 
    document.getElementById("flagstat-summary").innerHTML = '';
    
    // Reset maxContigLength at the start of a run
    maxContigLength = 1;

    try {
        const q = document.getElementById("query").files;
        const r = document.getElementById("ref").files;

        if (q.length === 0 || r.length === 0) {
            outputDiv.innerHTML = '<div class="alert alert-danger" role="alert">Please select both a **Reference** and **Query** file(s).</div>';
            document.getElementById("btn").disabled = false;
            return;
        }

        setProgress(5, "Loading files…");
        const query = await CLI.mount(q);
        const ref = await CLI.mount(r);

        // --- PIPELINE EXECUTION ---
        setProgress(10, "Aligning reads with minimap2 (longest step)…");
        const mmap_out = await CLI.exec("minimap2", [
            "-a",
            "-K", "50M", // smaller query batches -> lower peak memory (default 500M)
            "-o", "output.sam",
            "-x", "map-ont",
            ref[0],
            ...query
        ]);

        // Sort SAM directly into BAM (skips the intermediate unsorted BAM).
        // -m caps the per-thread sort buffer so samtools spills to temp files.
        setProgress(65, "Sorting alignments…");
        await CLI.exec("samtools sort -m 200M -o output.sorted.bam output.sam");
        // Free the SAM from the in-memory filesystem as soon as it's no longer needed
        await CLI.fs.unlink("output.sam");
        setProgress(80, "Indexing BAM…");
        await CLI.exec("samtools index output.sorted.bam");

        // Run samtools flagstat
        setProgress(85, "Counting mapped reads…");
        const flagstatOutput = await CLI.exec("samtools flagstat output.sorted.bam");
        
        // BAM download bug fix
        // 1. Read the file content from the Aioli filesystem as a Uint8Array
        const bamData = await CLI.fs.readFile("output.sorted.bam");
        // 2. Create a Blob object from the data
        const bamBlob = new Blob([bamData], { type: 'application/octet-stream' }); 
        // 3. Create a stable, revocable URL for the Blob
        // This URL remains valid as long as the Blob object exists
        const download_url = URL.createObjectURL(bamBlob);
        
        downloadButton.href = download_url;
        downloadButton.style.display = 'inline-block'; 

        
        //Reads stats 
        const totalReadsMatch = flagstatOutput.match(/^(\d+)\s*\+\s*\d+\s*in total/m);
        const secondaryReadsMatch = flagstatOutput.match(/^(\d+)\s*\+\s*\d+\s*secondary/m);
        const supplementaryReadsMatch = flagstatOutput.match(/^(\d+)\s*\+\s*\d+\s*supplementary/m);
        const mappedReadsMatch = flagstatOutput.match(/^(\d+)\s*\+\s*\d+\s*mapped\s+\((\d+\.\d+)%/m);
        
        let primaryReads = 'N/A';
        let primaryMapped = 'N/A';
        let primaryMappedPercent = 'N/A';

        // Get the required raw counts (QC-passed)
        const total = totalReadsMatch ? parseInt(totalReadsMatch[1], 10) : 0;
        const secondary = secondaryReadsMatch ? parseInt(secondaryReadsMatch[1], 10) : 0;
        const supplementary = supplementaryReadsMatch ? parseInt(supplementaryReadsMatch[1], 10) : 0;
        
        const mappedTotal = mappedReadsMatch ? parseInt(mappedReadsMatch[1], 10) : 0;
        
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
                <span class="fw-bold">Total fastq reads:</span> ${primaryReads} &nbsp;|&nbsp; 
                <span class="fw-bold">Mapped reads:</span> ${primaryMapped} 
                (${primaryMappedPercent}%)
            </div>
        `;
        document.getElementById("flagstat-summary").innerHTML = summaryHtml;
        
        // --- PER-CONTIG STATS (samtools coverage) ---
        // Columns: 0:#rname 1:startpos 2:endpos 3:numreads 4:covbases 5:coverage 6:meandepth 7:meanbaseq 8:meanmapq
        setProgress(90, "Computing per-contig coverage…");
        const covOutput = await CLI.exec("samtools coverage -H output.sorted.bam");
        const contigStats = [];
        for (const line of covOutput.trim().split('\n')) {
            const f = line.split('\t');
            if (f.length < 9) continue;
            const num = (x) => { const v = parseFloat(x); return isNaN(v) ? 0 : v; };
            contigStats.push({
                name: cleanContigName(f[0]),
                length: num(f[2]),
                numreads: num(f[3]),
                coverage: num(f[5]).toFixed(0),
                meandepth: num(f[6]).toFixed(0),
                meanmapq: num(f[8]).toFixed(0),
            });
        }

        if (contigStats.length === 0) {
            outputDiv.innerHTML = '<div class="alert alert-warning" role="alert">Analysis completed, but no coverage data was generated. Check your input files.</div>';
            downloadButton.style.display = 'none';
            return;
        }

        // --- BINNED DEPTH PROFILE (samtools bedcov over ~NBINS windows per contig) ---
        // bedcov reports the summed per-base depth in each window; dividing by the
        // window length gives mean depth, so dips/spikes aren't lost to subsampling.
        const bed = [];
        for (const c of contigStats) {
            const w = Math.max(1, Math.ceil(c.length / NBINS));
            for (let st = 0; st < c.length; st += w) {
                bed.push(`${c.name}\t${st}\t${Math.min(st + w, c.length)}`);
            }
            maxContigLength = Math.max(maxContigLength, c.length);
        }
        await CLI.fs.writeFile("windows.bed", bed.join('\n') + '\n');
        setProgress(95, "Building depth profiles…");
        const bedcovOutput = await CLI.exec("samtools bedcov windows.bed output.sorted.bam");

        const profiles = {};
        for (const line of bedcovOutput.trim().split('\n')) {
            const f = line.split('\t');
            if (f.length < 4) continue;
            const len = parseInt(f[2], 10) - parseInt(f[1], 10);
            const mean = len > 0 ? parseFloat(f[3]) / len : 0;
            (profiles[cleanContigName(f[0])] ??= []).push(Math.round(mean * 10) / 10);
        }

        // --- LIST POPULATION ---
        const listItemsHtml = [];
        for (const c of contigStats) {
            const depthsJSON = JSON.stringify(profiles[c.name] || []);
            const sparklineHtml = `<span class="sparkline-container" data-depths='${depthsJSON}' data-length='${c.length}'></span>`;

            listItemsHtml.push(`
                <li class="depth-item">
                    <div class="contig-name">${c.name}</div>
                    <div class="contig-length">${c.length.toLocaleString()}</div>
                    <div class="mean-depth">${c.meandepth}</div>
                    
                    <div class="coverage">${c.coverage}%</div>
                    <div class="mean-mapq">${c.meanmapq}</div>
                    <div class="depth-profile">${sparklineHtml}</div>
                </li>
            `);
        }
        
        ulList.innerHTML = listItemsHtml.join('');

        // --- SPARKLINE RENDERING ---
        renderSparklines(); 
        outputDiv.innerHTML = ''; // done: remove progress bar

    } catch (error) {
        console.error(error);
        // On error, display error message and hide the download button
        downloadButton.style.display = 'none';
        outputDiv.innerHTML = `<div class="alert alert-danger" role="alert">An error occurred during execution: **${error.message || String(error) || 'Unknown Aioli/CLI Error'}**</div>`;
    } finally {
        stopProgress();
        document.getElementById("btn").disabled = false;
    }
}

// Attach event listeners
document.getElementById("btn").addEventListener("click", run);
document.getElementById("reset-btn").addEventListener("click", resetPage);