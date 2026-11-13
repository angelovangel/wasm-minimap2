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
const DEPTH_MAX_BP = 25e6; // up to this total reference size, profiles come from one samtools depth pass
const SPARKLINE_COLOR = '#4c9aff';
const MIN_LINE_COLOR = '#c2255c'; // lowest depth per bin

// **NEW GLOBAL VARIABLE** to store max length for proportional scaling
let maxContigLength = 1; // Initialize to 1 to avoid division by zero

// Per-contig results of the last run (stats + binned depth), used by the reads view
let contigResults = [];
let referenceFile = null; // the reference fasta File, read again by the reads view
const contigSeqCache = new Map();

// Function to reload the page
function resetPage() {
    window.location.reload();
}

// Aggressive cleaning function to ensure keys match
function cleanContigName(name) {
    // Trim whitespace
    return name.trim();
}

// Format a base-pair position compactly: 950 bp, 12.5 kb, 3 Mb
function formatBp(bp) {
    if (bp >= 1e6) return `${+(bp / 1e6).toFixed(2)} Mb`;
    if (bp >= 1e3) return `${+(bp / 1e3).toFixed(1)} kb`;
    return `${bp} bp`;
}

// "Nice" tick step (1, 2 or 5 x 10^k) giving roughly `target` ticks over `length`
function niceStep(length, target) {
    const raw = length / Math.max(1, target);
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const n = raw / mag;
    return (n < 1.5 ? 1 : n < 3.5 ? 2 : n < 7.5 ? 5 : 10) * mag;
}

// Thin contig line with tick marks at nice positions
function buildAxisHtml(length, widthPct) {
    const step = Math.max(1, niceStep(length, Math.max(2, Math.round(widthPct / 20))));
    let ticks = '';
    // nice ticks, skipping any too close to the end so they don't collide with the length label
    for (let pos = 0; pos <= length - 0.75 * step; pos += step) {
        const pct = (pos / length) * 100;
        // keep labels inside the line: left-align at start, centre elsewhere
        const shift = pct < 5 ? '0' : '-50%';
        ticks += `<span class="axis-tick" style="left:${pct}%"></span>` +
                 `<span class="axis-label" style="left:${pct}%;transform:translateX(${shift})">${formatBp(pos)}</span>`;
    }
    // final tick: exact contig length, right-aligned to the end of the line
    ticks += `<span class="axis-tick" style="left:100%"></span>` +
             `<span class="axis-label axis-end" style="left:100%;transform:translateX(-100%)">${formatBp(length)}</span>`;
    return `<div class="contig-axis">${ticks}</div>`;
}

// Sparkline + contig line, both scaled to contig length relative to the longest contig
function renderSparklines() {
    $('#depthList .profile-track').each(function() {
        const $track = $(this);
        const $spark = $track.find('.sparkline-container');
        const depths = $spark.data('depths');
        const mins = $spark.data('mins');
        const contigLength = $track.data('length');

        const proportionalWidth = Math.max(20, (contigLength / maxContigLength) * 100);
        $track.css('width', `${proportionalWidth}%`);
        $track.append(buildAxisHtml(contigLength, proportionalWidth));

        if (depths && depths.length > 0) {
            const binSize = $track.data('bin') || Math.max(1, Math.ceil(contigLength / NBINS));
            const yMax = Math.max(1, ...depths);
            $spark.sparkline(depths, {
                type: 'line',
                width: '100%',
                height: '50px',
                lineColor: SPARKLINE_COLOR,
                fillColor: SPARKLINE_COLOR + '33',
                spotColor: false,
                minSpotColor: false,
                maxSpotColor: false,
                chartRangeMin: 0,
                chartRangeMax: yMax,
                highlightLineColor: '#6c757d',
                highlightSpotColor: '#6c757d',
                // tooltip shows the window's position range and its mean depth
                tooltipFormatter: (sp, opts, fields) => {
                    const st = fields.x * binSize;
                    const en = Math.min(st + binSize, contigLength);
                    return `${formatBp(st + 1)}–${formatBp(en)}<br>Mean depth: ${fields.y}`;
                }
            });
            // Lowest depth within each bin, drawn as a thin line over the mean
            // so narrow dips (deletions, gaps) aren't averaged away
            if (mins && mins.length === depths.length) {
                $spark.sparkline(mins, {
                    composite: true,
                    type: 'line',
                    lineColor: MIN_LINE_COLOR,
                    lineWidth: 1,
                    fillColor: false,
                    spotColor: false,
                    minSpotColor: false,
                    maxSpotColor: false,
                    highlightSpotColor: MIN_LINE_COLOR,
                    chartRangeMin: 0,
                    chartRangeMax: yMax,
                    tooltipFormatter: (sp, opts, fields) => `<br>Lowest: ${fields.y}`
                });
            }
        } else {
            $spark.text('No depth data');
        }
    });
}

// Depth profiles from one `samtools depth` pass (fast; true per-base minimum per bin).
// Without -a only covered positions are printed, so bins with missing positions have min 0.
// Deletions are not counted as coverage (samtools depth default), so they show as dips.
async function depthProfilesFromDepth(contigStats) {
    const out = await CLI.exec("samtools depth -d 0 output.sorted.bam");
    const byName = Object.fromEntries(contigStats.map(c => [c.name, {
        c, sum: new Float64Array(Math.ceil(c.length / c.binSize)),
        n: new Int32Array(Math.ceil(c.length / c.binSize)),
        min: new Float64Array(Math.ceil(c.length / c.binSize)).fill(Infinity)
    }]));
    let cur = null, curName = null;
    // Walk the output line by line without splitting the whole (large) string
    for (let i = 0; i < out.length;) {
        let nl = out.indexOf('\n', i);
        if (nl === -1) nl = out.length;
        const t1 = out.indexOf('\t', i), t2 = t1 === -1 ? -1 : out.indexOf('\t', t1 + 1);
        if (t1 !== -1 && t2 !== -1 && t2 < nl) {
            const name = out.slice(i, t1);
            if (name !== curName) { curName = name; cur = byName[cleanContigName(name)] || null; }
            if (cur) {
                const pos = parseInt(out.slice(t1 + 1, t2), 10) - 1;
                const d = parseInt(out.slice(t2 + 1, nl), 10);
                const b = Math.floor(pos / cur.c.binSize);
                if (b >= 0 && b < cur.sum.length && !isNaN(d)) {
                    cur.sum[b] += d; cur.n[b]++;
                    if (d < cur.min[b]) cur.min[b] = d;
                }
            }
        }
        i = nl + 1;
    }
    const round1 = x => Math.round(x * 10) / 10;
    const profiles = {}, minProfiles = {};
    for (const { c, sum, n, min } of Object.values(byName)) {
        const mean = [], low = [];
        for (let b = 0; b < sum.length; b++) {
            const binLen = Math.min(c.binSize, c.length - b * c.binSize);
            mean.push(round1(sum[b] / binLen));
            low.push(n[b] < binLen ? 0 : round1(min[b])); // unreported positions have depth 0
        }
        profiles[c.name] = mean;
        minProfiles[c.name] = low;
    }
    return { profiles, minProfiles };
}

// Depth profiles for large references: samtools bedcov over NBINS windows per contig
// (mean only; per-base output would be too large). -j: don't count deletions as coverage.
async function depthProfilesFromBedcov(contigStats) {
    const bed = [];
    for (const c of contigStats) {
        for (let st = 0; st < c.length; st += c.binSize) {
            bed.push(`${c.name}\t${st}\t${Math.min(st + c.binSize, c.length)}`);
        }
    }
    await CLI.fs.writeFile("windows.bed", bed.join('\n') + '\n');
    let out = await CLI.exec("samtools bedcov -j windows.bed output.sorted.bam");
    if (!/^\S+\t\d+\t\d+\t\d+/.test(out)) {
        // older samtools without -j: fall back to counting deletions as covered
        console.warn("samtools bedcov -j failed, retrying without -j:", out.slice(0, 200));
        out = await CLI.exec("samtools bedcov windows.bed output.sorted.bam");
    }
    await CLI.fs.unlink("windows.bed");
    const profiles = {};
    for (const line of out.trim().split('\n')) {
        const f = line.split('\t');
        if (f.length < 4) continue;
        const len = parseInt(f[2], 10) - parseInt(f[1], 10);
        const mean = len > 0 ? parseFloat(f[3]) / len : 0;
        (profiles[cleanContigName(f[0])] ??= []).push(Math.round(mean * 10) / 10);
    }
    return { profiles, minProfiles: {} };
}

// Progress bar: the wasm tools don't report internal progress, so the bar
// advances per pipeline step and an elapsed timer shows it's still working.
let progressTimer = null;

function showProgress(outputDiv) {
    outputDiv.innerHTML = `
        <div class="run-progress">
            <div class="run-progress-text">
                <span id="progress-label">Starting…</span>
                <span id="progress-elapsed">0s</span>
            </div>
            <div class="run-progress-track" role="progressbar" aria-label="Analysis progress" aria-valuemin="0" aria-valuemax="100">
                <div id="progress-bar" class="run-progress-fill" style="width: 0%"></div>
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

        const subsamplePct = parseFloat(document.getElementById("subsample").value);
        if (!(subsamplePct > 0 && subsamplePct <= 100)) {
            outputDiv.innerHTML = '<div class="alert alert-danger" role="alert">Proportion of reads must be between 0 and 100%.</div>';
            return;
        }

        const preset = document.getElementById("preset").value;

        setProgress(5, "Loading files…");
        let query = await CLI.mount(q);
        const ref = await CLI.mount(r);
        referenceFile = r[0];
        contigSeqCache.clear();

        // --- OPTIONAL SUBSAMPLING (seqtk sample, fixed seed for reproducibility) ---
        // seqtk writes to stdout, so its output is written back into the
        // in-memory filesystem as a plain fastq that minimap2 can read.
        const subsampled = [];
        if (subsamplePct < 100) {
            const frac = String(subsamplePct / 100);
            for (let i = 0; i < query.length; i++) {
                setProgress(5, `Subsampling ${subsamplePct}% of reads (file ${i + 1}/${query.length})…`);
                const fq = await CLI.exec("seqtk", ["sample", "-s", "11", query[i], frac]);
                const name = `subsample_${i}.fastq`;
                await CLI.fs.writeFile(name, fq);
                subsampled.push(name);
            }
            query = subsampled;
        }

        // --- PIPELINE EXECUTION ---
        setProgress(10, "Aligning reads with minimap2 (longest step)…");
        const mmap_out = await CLI.exec("minimap2", [
            "-a",
            "-K", "50M", // smaller query batches -> lower peak memory (default 500M)
            "-o", "output.sam",
            "-x", preset,
            ref[0],
            ...query
        ]);
        // Subsampled fastqs are no longer needed once aligned
        for (const f of subsampled) await CLI.fs.unlink(f);

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

        // Unmapped records are always primary, so total - mapped is exact
        const unmapped = Math.max(0, total - mappedTotal);
        const summaryHtml = `
            <div class="alert alert-info summary-box" role="alert">
                <span class="fw-bold">${subsamplePct < 100 ? `Reads used (${subsamplePct}% subsample)` : "Total fastq reads"}:</span> ${primaryReads} &nbsp;|&nbsp; 
                <span class="fw-bold">Mapped reads:</span> ${primaryMapped} 
                (${primaryMappedPercent}%) &nbsp;|&nbsp;
                <span class="fw-bold">Unmapped:</span> ${unmapped.toLocaleString()}
                <span id="warn-summary"></span>
            </div>
        `;
        document.getElementById("flagstat-summary").innerHTML = summaryHtml;
        const unmappedBtn = document.getElementById("download-unmapped");
        unmappedBtn.style.display = unmapped > 0 ? 'inline-block' : 'none';
        
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

        // --- BINNED DEPTH PROFILE ---
        // Per contig, <= NBINS bins with the mean depth and the lowest depth in each bin,
        // so a dip narrower than a bin (e.g. a deletion) still shows up.
        setProgress(95, "Building depth profiles…");
        for (const c of contigStats) {
            c.binSize = Math.max(1, Math.ceil(c.length / NBINS));
            maxContigLength = Math.max(maxContigLength, c.length);
        }
        const totalRefBp = contigStats.reduce((t, c) => t + c.length, 0);
        const { profiles, minProfiles } = totalRefBp <= DEPTH_MAX_BP
            ? await depthProfilesFromDepth(contigStats)
            : await depthProfilesFromBedcov(contigStats);

        // --- CONTIG WARNINGS ---
        for (const c of contigStats) {
            c.warnings = contigWarnings(c, minProfiles[c.name]);
        }
        const nWarn = contigStats.filter(c => c.warnings.length).length;
        document.getElementById("warn-summary").innerHTML = nWarn
            ? `&nbsp;|&nbsp; <span class="warn-badge">!</span>${nWarn} of ${contigStats.length} contigs with warnings (hover the ! for details)`
            : '';

        // --- LIST POPULATION ---
        const listItemsHtml = [];
        contigResults = contigStats.map(c => ({ ...c, profile: profiles[c.name] || [] }));
        document.getElementById("download-csv").style.display = 'inline-block';

        for (const [i, c] of contigStats.entries()) {
            const depthsJSON = JSON.stringify(profiles[c.name] || []);
            const minsJSON = JSON.stringify(minProfiles[c.name] || []);
            const sparklineHtml = `<div class="profile-track" data-length='${c.length}' data-bin='${c.binSize}'><span class="sparkline-container" data-depths='${depthsJSON}' data-mins='${minsJSON}'></span></div>`;

            listItemsHtml.push(`
                <li class="depth-item clickable" data-index="${i}" title="Click to view reads for ${c.name}">
                    <div class="contig-name">${c.warnings.length ? `<span class="warn-badge" title="${c.warnings.join('\n')}">!</span>` : ''}${c.name}</div>
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
        // Results are shown: hide the inputs and Run button (RESET starts over)
        document.getElementById("inputs").style.display = 'none';
        document.getElementById("btn").style.display = 'none';

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

// --- CONTIG WARNINGS ---
const WARN_MIN_COVERAGE = 90; // % of bases covered
const WARN_MIN_DEPTH = 5;     // mean depth

// Returns a list of human-readable warnings for one contig
function contigWarnings(c, mins) {
    const w = [];
    const cov = parseFloat(c.coverage), depth = parseFloat(c.meandepth);
    if (c.numreads === 0) return ['no reads mapped'];
    if (cov < WARN_MIN_COVERAGE) w.push(`only ${cov}% of bases covered`);
    if (depth < WARN_MIN_DEPTH) w.push(`low mean depth (${depth}x)`);
    // zero-depth stretches away from the contig ends (first/last bins naturally thin out)
    if (mins && mins.length > 2) {
        let stretches = 0;
        for (let i = 1; i < mins.length - 1; i++) {
            if (mins[i] === 0 && mins[i - 1] !== 0) stretches++;
        }
        if (stretches) w.push(`${stretches} zero-depth stretch${stretches > 1 ? 'es' : ''} inside the contig`);
    }
    return w;
}

// --- DOWNLOADS ---
function downloadBlob(blob, filename) {
    const url = URL.createObjectURL(blob);
    const a = Object.assign(document.createElement('a'), { href: url, download: filename });
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 10000);
}

const csvCell = (v) => /[",\n]/.test(String(v)) ? `"${String(v).replace(/"/g, '""')}"` : String(v);

function downloadStatsCsv() {
    const header = ['contig', 'length', 'reads', 'covered_pct', 'mean_depth', 'mean_mapq', 'warnings'];
    const lines = [header.join(',')];
    for (const c of contigResults) {
        lines.push([c.name, c.length, c.numreads, c.coverage, c.meandepth, c.meanmapq, (c.warnings || []).join('; ')]
            .map(csvCell).join(','));
    }
    downloadBlob(new Blob([lines.join('\n') + '\n'], { type: 'text/csv' }), 'contig_stats.csv');
}

// Unmapped reads as FASTQ. samtools fastq writes reads to a file (-0) so its
// stderr messages can't end up mixed into the FASTQ.
async function downloadUnmapped() {
    const btn = document.getElementById("download-unmapped");
    const label = btn.textContent;
    btn.disabled = true;
    btn.textContent = 'Extracting…';
    try {
        await CLI.exec("samtools", ["fastq", "-f", "4", "-0", "unmapped.fastq", "output.sorted.bam"]);
        const data = await CLI.fs.readFile("unmapped.fastq");
        await CLI.fs.unlink("unmapped.fastq");
        downloadBlob(new Blob([data], { type: 'text/plain' }), 'unmapped.fastq');
    } catch (err) {
        console.error(err);
        alert(`Could not extract unmapped reads: ${err.message || err}`);
    } finally {
        btn.disabled = false;
        btn.textContent = label;
    }
}

document.getElementById("download-csv").addEventListener("click", downloadStatsCsv);
document.getElementById("download-unmapped").addEventListener("click", downloadUnmapped);

// --- READS VIEW (opens reads.html in a new tab) ---
const MAX_VIEW_READS = 3000; // above this, reads are randomly subsampled for display
const MIN_GAP_BP = 20;       // deletions/introns at least this long are drawn as gaps
const MIN_INS_BP = 20;       // insertions at least this long get a marker
const EXCLUDE_FLAGS = "0x704"; // unmapped, secondary, QC-fail, duplicate

const MAX_SEQ_BP = 30e6;     // longer contigs are shown without mismatches/consensus

// Stream the reference fasta and return one contig's sequence (uppercase), without
// holding the whole file in memory. Header name = text after '>' up to whitespace.
async function readContigSeq(file, name) {
    if (contigSeqCache.has(name)) return contigSeqCache.get(name);
    const reader = file.stream().pipeThrough(new TextDecoderStream()).getReader();
    const parts = [];
    let inTarget = false, done = false, rest = '';
    while (!done) {
        const { value, done: eof } = await reader.read();
        const text = rest + (value || '');
        const lines = text.split('\n');
        rest = eof ? '' : lines.pop();
        for (const raw of lines) {
            const line = raw.trimEnd();
            if (line.startsWith('>')) {
                if (inTarget) { done = true; break; }
                inTarget = line.slice(1).split(/\s/)[0] === name;
            } else if (inTarget) {
                parts.push(line);
            }
        }
        if (eof) break;
    }
    reader.cancel().catch(() => {});
    const seq = parts.join('').toUpperCase();
    contigSeqCache.set(name, seq);
    return seq;
}

// Parse one SAM line into a compact read record. Only the first 10 fields are split,
// so the long QUAL column is never copied. With the contig sequence `ref`, read bases
// are compared to the reference to find mismatches.
//   blocks: aligned segments [s0,e0,s1,e1,...], split at deletions >= MIN_GAP_BP
//   mmPos/mmBase: mismatch positions and read bases
//   dels: every deletion as [pos,len,...]; ins: every insertion as [pos,len,...]
function parseSamRead(line, ref) {
    const f = line.split('\t', 10);
    if (f.length < 6 || f[5] === '*') return null;
    const flag = parseInt(f[1], 10);
    const start = parseInt(f[3], 10) - 1; // SAM is 1-based
    if (isNaN(flag) || isNaN(start)) return null;
    const seq = ref && f[9] && f[9] !== '*' ? f[9] : null;

    const blocks = [], dels = [], ins = [], mmPos = [];
    let mmBase = '';
    let pos = start, qpos = 0, blockStart = start, clipL = 0, clipR = 0;
    for (const [, lenStr, op] of f[5].matchAll(/(\d+)([MIDNSHP=X])/g)) {
        const len = parseInt(lenStr, 10);
        switch (op) {
            case 'M': case '=': case 'X':
                if (seq) {
                    for (let i = 0; i < len; i++) {
                        const rb = ref.charCodeAt(pos + i), qb = seq.charCodeAt(qpos + i) & 0xDF; // uppercase
                        if (qb !== rb && rb !== 78 && qb !== 78) { mmPos.push(pos + i); mmBase += String.fromCharCode(qb); }
                    }
                }
                pos += len; qpos += len;
                break;
            case 'I':
                ins.push(pos, len);
                qpos += len;
                break;
            case 'D': case 'N':
                if (op === 'D') dels.push(pos, len);
                if (len >= MIN_GAP_BP || op === 'N') { blocks.push(blockStart, pos); blockStart = pos + len; }
                pos += len;
                break;
            case 'S':
                // leading clip if nothing has been aligned yet, otherwise trailing
                if (pos === start && qpos === 0) clipL = len; else clipR = len;
                qpos += len;
                break;
        }
    }
    blocks.push(blockStart, pos);

    return {
        name: f[0], start, end: pos, mapq: parseInt(f[4], 10) || 0,
        reverse: (flag & 0x10) !== 0, supplementary: (flag & 0x800) !== 0,
        blocks, dels: Int32Array.from(dels), ins: Int32Array.from(ins),
        mmPos: Int32Array.from(mmPos), mmBase, clipL, clipR
    };
}

async function fetchContigReads(contig) {
    const countOut = await CLI.exec("samtools", ["view", "-c", "-F", EXCLUDE_FLAGS, "output.sorted.bam", contig.name]);
    const total = parseInt(countOut.trim().split('\n').pop(), 10) || 0;

    const args = ["view", "-F", EXCLUDE_FLAGS];
    let fraction = 1;
    if (total > MAX_VIEW_READS) {
        fraction = MAX_VIEW_READS / total;
        // samtools -s INT.FRAC: integer part is the seed, fractional part the fraction kept
        args.push("-s", `42${fraction.toFixed(6).slice(1)}`);
    }
    args.push("output.sorted.bam", contig.name);
    const sam = await CLI.exec("samtools", args);

    // Reference bases for mismatch calling (skipped for very long contigs)
    let ref = '';
    if (referenceFile && contig.length <= MAX_SEQ_BP) {
        ref = await readContigSeq(referenceFile, contig.name);
        if (ref.length !== contig.length) ref = ''; // name/length mismatch: don't guess
    }

    const reads = [];
    for (const line of sam.split('\n')) {
        const r = parseSamRead(line, ref || null);
        if (r) reads.push(r);
    }
    return { total, fraction, reads, ref };
}

// One pending payload per opened tab, delivered when that tab reports it's ready
const readsViews = new Map(); // window -> { ready, payload }

function deliverReads(win) {
    const v = readsViews.get(win);
    if (v && v.ready && v.payload) {
        win.postMessage({ type: 'reads-data', ...v.payload }, window.location.origin);
        readsViews.delete(win);
    }
}

window.addEventListener('message', (e) => {
    if (e.origin !== window.location.origin || e.data?.type !== 'reads-ready') return;
    const v = readsViews.get(e.source);
    if (v) { v.ready = true; deliverReads(e.source); return; }
    // Unknown tab (e.g. the reads tab was reloaded): fetch again for the contig it names
    const contig = contigResults.find(c => c.name === e.data.contig);
    if (contig) sendReadsTo(e.source, contig, true);
});

async function openReadsView(contig) {
    // Open the tab synchronously in the click handler so popup blockers allow it
    const win = window.open(`reads.html#${encodeURIComponent(contig.name)}`, '_blank');
    if (!win) {
        alert('Please allow pop-ups for this page to open the reads view.');
        return;
    }
    sendReadsTo(win, contig, false);
}

async function sendReadsTo(win, contig, ready) {
    readsViews.set(win, { ready, payload: null });
    let payload;
    try {
        const { total, fraction, reads, ref } = await fetchContigReads(contig);
        payload = {
            contig: { name: contig.name, length: contig.length, meandepth: contig.meandepth, coverage: contig.coverage, binSize: contig.binSize },
            profile: contig.profile, total, fraction, reads, ref, limit: MAX_VIEW_READS
        };
    } catch (err) {
        console.error(err);
        payload = { error: String(err?.message || err) };
    }
    const v = readsViews.get(win);
    if (v) { v.payload = payload; deliverReads(win); }
}

document.getElementById("depthList").addEventListener("click", (e) => {
    const item = e.target.closest('.depth-item[data-index]');
    if (!item) return;
    const contig = contigResults[parseInt(item.dataset.index, 10)];
    if (contig) openReadsView(contig);
});

// Attach event listeners
document.getElementById("btn").addEventListener("click", run);
document.getElementById("reset-btn").addEventListener("click", resetPage);