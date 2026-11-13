// Results tables start with 5 rows; a "Rows per page" menu beside the pager shows
// more (or all). The rows scroll in a box beneath a frozen header row (see
// .tableScroll and .tableFooter in styles.css).
// The not-evaluated tables set their own dom and keep their single page.
$.extend(true, $.fn.dataTable.defaults, {
    dom: '<"tableScroll"t><"tableFooter"ilp>',
    pageLength: 5,
    lengthMenu: [[5, 10, 25, 50, -1], [5, 10, 25, 50, 'All']],
    language: {lengthMenu: 'Rows per page: _MENU_'}
});


// Clicking an example adds it to the input box and selects its organism.
// Variants of the other organism can't be run together, so they are replaced.
document.querySelectorAll('.exampleVariant').forEach(function (example) {
    example.addEventListener('click', function () {
        const input = document.getElementById('variantInput');
        const organism = document.getElementById('queryOrganism');
        const current = input.value.trim();

        if (!current || organism.value !== example.dataset.organism) {
            input.value = example.dataset.variant;
        } else if (!current.split('\n').map(s => s.trim()).includes(example.dataset.variant)) {
            input.value = current + '\n' + example.dataset.variant;
        }
        organism.value = example.dataset.organism;
        input.focus();
    });
});

// Run variant-model mapping for submitted variants on click.
// Clear empties the variant box; results already shown stay until the next submission
document.getElementById('clearButton').addEventListener('click', function () {
    const input = document.getElementById('variantInput');
    input.value = '';
    showError('');
    input.focus();
});


document.getElementById('runButton').addEventListener('click', function () {
    const vars = document.getElementById('variantInput').value;
    const org = document.getElementById('queryOrganism').value;

    showError('');
    const stopStatus = showStatus(vars.split('\n').filter(line => line.trim()).length);

    fetch('/run_variants', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({ variants: vars, organism: org})
    })
    .then(async res => {
        // The server sends {error: message} when the analysis fails
        const data = await res.json().catch(() => ({}));
        if (!res.ok || data.error) {
            throw new Error(data.error || 'The analysis could not be completed. Please confirm that the Organism setting matches your coordinates.');
        }
        return data;
    })
    .then(data => {
        // store full results, with what was run so a downloaded copy can be reopened
        data.organism = org;
        data.input = vars;
        data.run_at = new Date().toISOString();
        window.fullJSON = data;

        // build sidebar; a failure here is in the page, not the analysis
        try {
            buildSidebar(data);
        } catch (err) {
            throw new Error(`The analysis finished, but the results could not be displayed (${err.message}). Try reloading the page.`);
        }
    })
    .catch(err => {
        console.error(err);

        // Hide results from any earlier search so they aren't mistaken for this one
        window.fullJSON = undefined;
        document.getElementById('layoutOutput').style.display = 'none';
        document.getElementById('resultsPanel').hidden = true;
        showError(err.message);
    })
    .finally(stopStatus);
});


// While a query runs: hide earlier results, show the status panel with a running
// clock, and disable Submit so the query isn't sent twice. Returns the function
// that puts everything back.
function showStatus(variantCount) {
    const panel = document.getElementById('statusPanel');
    const elapsed = document.getElementById('statusElapsed');
    const button = document.getElementById('runButton');
    const start = Date.now();

    document.getElementById('resultsPanel').hidden = true;
    document.getElementById('layoutOutput').style.display = 'none';

    document.getElementById('statusText').textContent =
        `Mapping ${variantCount} ${variantCount === 1 ? 'variant' : 'variants'} with Ensembl VEP, ClinVar and MGI. ` +
        'This can take a minute or more for genes with many alleles or ClinVar variants.';
    elapsed.textContent = '0:00';
    panel.hidden = false;

    button.disabled = true;
    button.textContent = 'Running…';

    const timer = setInterval(() => {
        const seconds = Math.floor((Date.now() - start) / 1000);
        elapsed.textContent = `${Math.floor(seconds / 60)}:${String(seconds % 60).padStart(2, '0')}`;
    }, 1000);

    return () => {
        clearInterval(timer);
        panel.hidden = true;
        button.disabled = false;
        button.textContent = 'Submit';
    };
}


// Show a message under the Submit button; an empty message hides it
function showError(message) {
    const box = document.getElementById('errorMessage');
    box.textContent = message;
    box.hidden = !message;
}



// Sort score rows so each allele's transcripts stay together. Alleles (ClinVar variants
// on mouse queries) are ordered by their Best Tier, then by Disease Match; within an
// allele, transcripts by their own tier, protein-coding transcripts first.
function groupedScoreSort(scores) {
    const key = row => String(row['AlleleID'] ?? row['Name'] ?? '');
    const tier = value => (value === null || value === undefined) ? 99 : Number(value);
    const disease = {};
    scores.forEach(row => { if (row['Disease Match'] === true) disease[key(row)] = true; });

    return [...scores].sort((a, b) =>
        tier(a['Best Tier']) - tier(b['Best Tier'])
        || (disease[key(b)] ? 1 : 0) - (disease[key(a)] ? 1 : 0)
        || key(a).localeCompare(key(b))
        || tier(a['Equivalency Tier']) - tier(b['Equivalency Tier'])
        || (b['Counts Toward Tier'] ? 1 : 0) - (a['Counts Toward Tier'] ? 1 : 0));
}


// Build the sidebar with one button per input variant
function buildSidebar(data) {
    const sidebar = document.getElementById('sidebar');

    // Show the results panel beside the input box; the tables appear when a variant is clicked
    document.getElementById('resultsPanel').hidden = false;
    document.getElementById('layoutOutput').style.display = 'none';

    const queryVariants = data.query_variants || [];
    const unevaluated = data.unevaluated_variants || [];

    // Create sidebar header
    sidebar.innerHTML = queryVariants.length > 0 || unevaluated.length === 0 ? '<h3>Variant Inputs:</h3>' : '';

    // Handle no results
    if (queryVariants.length === 0 && unevaluated.length === 0) {
        sidebar.innerHTML += '<h2>No inputs found.</h2>';
        return;
    }

    // Labels entered after a variant are shown beside it
    const labels = data.labels || {};
    function addButton(name, render) {
        const btn = document.createElement('button');
        btn.className = 'sidebar-button';
        btn.textContent = name;
        if (labels[name]) {
            const tag = document.createElement('span');
            tag.className = 'variantLabel';
            tag.textContent = labels[name];
            btn.append(' ', tag);
        }
        btn.onclick = () => {
            // deactivate all other buttons
            document.querySelectorAll('#sidebar .sidebar-button')
            .forEach(b => b.classList.remove('active'));

            // activate clicked button
            btn.classList.add('active')

            // render output tables and jump down to them
            render();
            const output = document.getElementById('layoutOutput');
            output.style.display = 'block';
            output.scrollIntoView({behavior: 'smooth', block: 'start'});
        }
        sidebar.appendChild(btn);
    }

    // Create sidebar buttons
    queryVariants.forEach(row => addButton(row['Name'], () => renderTables(row['Name'], row['Gene Symbol'])));

    // Variants without an amino-acid change get their own group
    const unevaluatedNames = [...new Set(unevaluated.map(row => row['Input Name']))];
    if (unevaluatedNames.length > 0) {
        sidebar.insertAdjacentHTML('beforeend', '<h3>Not Evaluated:</h3>');
        unevaluatedNames.forEach(name => addButton(name, () => renderUnevaluated(name)));
    }
}


// Back to the variant list from the results tables
document.querySelector('.backToVariants').addEventListener('click', function (event) {
    event.preventDefault();
    document.querySelector('.topRow').scrollIntoView({behavior: 'smooth', block: 'start'});
});


// Render all tables for the selected input
function renderTables(name, gene) {
    const org = document.getElementById('queryOrganism').value;

    // Define organism header identifiers
    if (org === 'human') {
        queryOrganism = 'Human';
        queryShort = 'Hum';
        targetOrganism = 'Mouse';
        targetShort = 'Mus';
    } else if (org === 'mouse') {
        queryOrganism = 'Mouse';
        queryShort = 'Mus';
        targetOrganism = 'Human';
        targetShort = 'Hum';
    }


    // Extract homolog
    homolog = window.fullJSON.gene_mapping.find(row => row[`${queryShort} Gene`] === gene)?.[`${targetShort} Gene`]


    // Extract relevant tables
    if (org === 'human') {
        queryGenes = window.fullJSON.query_genes.filter(row => row['Gene Symbol'] === gene);
        queryProteins = window.fullJSON.query_variants.filter(row => row['Name'] === name);
        outputGenes = window.fullJSON.output_genes.filter(row => row['Gene Symbol'] === homolog);
        outputProteins = window.fullJSON.output_variants.filter(row => row['Gene Symbol'] === homolog);
        scores = window.fullJSON.scores.filter(row => row['Input Name'] === name);
        // The coverage table lists every allele of the homolog with its phenotypes,
        // shown in dropdowns, so human queries need no separate phenotype table
        alleleCoverage = (window.fullJSON.allele_coverage || []).filter(row => row['Gene Symbol'] === homolog);
        phenotypesByAllele = new Map(alleleCoverage.map(row => [row['AlleleID'], row]));

        // Sort output proteins and score dataframes
        sortedScores = groupedScoreSort(scores);
        const rank = new Map(sortedScores.map((row, index) => [`${row['AlleleID']}||${row['Transcript ID']}`, index]));

        console.log(rank)

        sortedOutputProteins = [...outputProteins].sort((a, b) =>
            (rank.get(`${a['AlleleID']}||${a['Transcript ID']}`) ?? Infinity) -
            (rank.get(`${b['AlleleID']}||${b['Transcript ID']}`) ?? Infinity)
        );

        console.log(outputProteins)
        console.log(sortedOutputProteins)


    }
    else if (org === 'mouse') {
        queryGenes = window.fullJSON.query_genes.filter(row => row['Gene Symbol'] === gene);
        queryProteins = window.fullJSON.query_variants.filter(row => row['Name'] === name);
        outputGenes = window.fullJSON.output_genes.filter(row => row['Gene Symbol'] === homolog);
        // Coding variants match by amino-acid change, splice-site and intronic variants by intron offset
        outputProteins = window.fullJSON.output_variants.filter(row => row['Gene Symbol'] === homolog
                                                                        && queryProteins.some(r => r['Intron Offset']
                                                                            ? r['Intron Offset'] === row['Intron Offset']
                                                                            : r.refAA === row.refAA && r.varAA === row.varAA));
        scores = window.fullJSON.scores.filter(row => row['Input Name'] === name);
        phenotypes = window.fullJSON.phenotypes.filter(row => outputProteins.some(r => r.Name === row.Name));

        // Mouse queries list the queried gene's alleles, marking the ones this variant matches
        alleleCoverage = markEvaluated(window.fullJSON.allele_coverage, gene, queryProteins);
        phenotypesByAllele = new Map(alleleCoverage.map(row => [row['AlleleID'], row]));

        // Sort output proteins and score dataframes
        sortedScores = groupedScoreSort(scores);
        const rank = new Map(sortedScores.map((row, index) => [`${row['Name']}||${row['Transcript ID']}`, index]));

        console.log(rank)

        sortedOutputProteins = [...outputProteins].sort((a, b) =>
            (rank.get(`${a['Name']}||${a['Transcript ID']}`) ?? Infinity) -
            (rank.get(`${b['Name']}||${b['Transcript ID']}`) ?? Infinity)
        );
    }





    // Query variant fields that are not shown (MONDO_set only on human queries)
    const queryVariantHidden = ['Name', 'Submission', 'Gene Symbol', 'Biotype', 'refAA', 'protein_start', 'varAA',
                                ...(org === 'human' ? ['MONDO_set'] : [])];

    // Mouse Allele Coverage: for human queries it comes before the mouse variants, as the
    // overview of how many alleles could be evaluated ahead of the per-allele detail
    const coverage = coverageSection(alleleCoverage, org === 'human' ? homolog : gene, org);

    // Create HTML tables
    results.innerHTML = `
        <div class='tableWrapper'>
            ${infoHeading('submittedInfo', `Submitted ${queryOrganism} Variant`, INFO.submitted(org))}
            ${jsonToFieldList(queryProteins.map(row => withLabel(name, variantFields(row))), "queryVariantTable", queryVariantHidden)}
        </div>

        <div class='tableWrapper'>
            ${infoHeading('genesInfo', 'Genes', INFO.genes(org))}
            ${jsonToFieldList([...queryGenes, ...outputGenes].map(geneFields), "genesTable", [],
                              [...queryGenes.map(() => `${queryOrganism} Gene`), ...outputGenes.map(() => `${targetOrganism} Homolog`)])}
        </div>

        ${org === 'human' ? coverage : ''}

        <div class='tableWrapper'>
            ${infoHeading('variantsInfo', `${targetOrganism} Variants Evaluated`, INFO.variantsEvaluated(org))}
            ${jsonToHTMLTable(sortedOutputProteins.map(row => withAlignButton(displayFields(row), row, org)), "outputVariantTable")}
        </div>

        <div class='tableWrapper'>
            ${infoHeading('scoringInfo', 'Variant Equivalency Scoring', INFO.scoring(org))}
            ${jsonToHTMLTable(sortedScores, "scoresTable")}
        </div>

        ${org === 'human' ? phenotypeSimilaritySection(name, gene, sortedScores) : ''}


        ${org === 'human' ? '' : coverage}

        ${org === 'human' ? '' : `
        <div class='tableWrapper'>
            ${infoHeading('phenotypesInfo', `${targetOrganism} Phenotypes`, INFO.humanPhenotypes())}
            ${jsonToHTMLTable(phenotypes, "phenotypeTable")}
        </div>`}
    `;


    // Function to hide columns in DataTables by name
    function hideColumnsByName(table, names) {
        table.columns().every(function () {
            const header = $(this.header()).text().trim();
            if (names.includes(header)) {
                this.visible(false);
            }
        });
    }


    // CREATE OUTPUT VARIANT DATATABLE ELEMENT
    const outputVariantTable = $('#outputVariantTable').DataTable({
        ordering: false
    });

    // Columns that are be hidden by default
    const outputHiddenCols = [
                    'Submission',
                    'Input Name',
                    'Gene Symbol',
                    'Biotype',
                    'refAA',
                    'varAA'
                ];

    // Extra human-only columns to hide
    const outputHumanHidden = [
                    'protein_start'
                ];

    // Extra mouse-only columns to hide
    const outputMouseHidden = [
                    'MONDO_set'
                ];

    // Hide columns
    hideColumnsByName(outputVariantTable, outputHiddenCols);

    if (targetOrganism === 'Human') {hideColumnsByName(outputVariantTable, outputHumanHidden);}

    if (targetOrganism === 'Mouse') {hideColumnsByName(outputVariantTable, outputMouseHidden);}

    // Align buttons open the transcript alignment beneath the row
    initAlignmentToggles(outputVariantTable, org, queryProteins, scores);
    

    // CREATE SCORE DATATABLE ELEMENT
    const scoreTable = $('#scoresTable').DataTable({
        ordering: false
    });

    // Columns that are be hidden by default
    // AlleleID is linked from the Allele Symbol instead. Phenotype similarity has its own
    // table below (one row per allele); its columns stay in the JSON download
    // The tier's basis is shown in the tier cell. Total Score is superseded by the tier and
    // stays in the download only
    const scoreHiddenCols = ['Input Name', 'AlleleID', 'Matched Disease(s)', 'Phenotype Score (variant)',
                             'Phenotype Score (gene)', 'MP Terms Mapped', 'MP Terms Total',
                             'Tier Basis', 'Counts Toward Tier', 'Best Tier', 'Total Score'];

    // Criteria that don't apply to any row, such as the coding criteria for a splice-site variant
    const notApplicable = Object.keys(sortedScores[0] || {})
        .filter(col => col.endsWith('Match') && sortedScores.every(row => row[col] === null));

    // Hide columns
    // The query transcript matters only when the variant falls in several of them
    const queryTranscripts = new Set(sortedScores.map(row => row['Query Transcript ID']));
    const singleQueryTranscript = queryTranscripts.size <= 1 ? ['Query Transcript ID'] : [];
    hideColumnsByName(scoreTable, [...scoreHiddenCols, ...notApplicable, ...singleQueryTranscript]);

    // Aligned Position buttons open the transcript overview beneath the row
    initTranscriptOverview(scoreTable, org, queryProteins, sortedOutputProteins, sortedScores,
                           org === 'human' ? homolog : gene);


    // CREATE PHENOTYPE DATATABLE ELEMENT
    $('#phenotypeTable').DataTable({
        ordering: false
    });


    initCoverageTable(phenotypesByAllele);
    initPhenotypeSimilarityTable(name);

}


// Allele rows for a mouse gene, with Evaluated marking the alleles the
// given variant rows match (the server marks every allele in the batch)
function markEvaluated(alleleRows, gene, variantRows) {
    const ids = new Set(variantRows.flatMap(row => String(row['AlleleID'] || '').split(', ')));
    return (alleleRows || [])
        .filter(row => row['Gene Symbol'] === gene)
        .map(row => ({...row, Evaluated: ids.has(row['AlleleID'])}));
}


// Open or close an allele's phenotype list beneath its row
function initPhenotypeToggles(tableId, dataTable, phenotypesByAllele) {
    $(`#${tableId} tbody`).on('click', 'button.phenotypeToggle', function () {
        const row = dataTable.row($(this).closest('tr'));
        if (row.child.isShown()) {
            row.child.hide();
            this.setAttribute('aria-expanded', 'false');
        } else {
            const allele = phenotypesByAllele.get(this.dataset.allele);
            row.child(`<div class="phenotypeList">${phenotypeLists(allele.Phenotypes, allele['Normal Phenotypes'])}</div>`).show();
            this.setAttribute('aria-expanded', 'true');
        }
    });
}


// PROTEIN DOMAINS on the alignment views: one colour per Pfam domain, shared by both views
const DOMAIN_COLORS = ['#2b9a9a', '#c9a227', '#c2577a', '#3fa7d6', '#8d6e63', '#7a8b2e'];

function domainColors(domainLists) {
    const colors = new Map();
    domainLists.flat().forEach(d => {
        if (!colors.has(d.id)) colors.set(d.id, DOMAIN_COLORS[colors.size % DOMAIN_COLORS.length]);
    });
    return colors;
}

function domainLegend(domainLists, colors) {
    const names = new Map();
    domainLists.flat().forEach(d => { if (!names.has(d.id)) names.set(d.id, d.name); });
    if (names.size === 0) return '<p class="alignNote">No Pfam domains are annotated on these proteins.</p>';
    return `<p class="domainLegend">${[...names].map(([id, name]) =>
        `<span class="domainKey"><span class="domainSwatch" style="background:${colors.get(id)}"></span>`
        + `<a href="https://www.ebi.ac.uk/interpro/entry/pfam/${escapeText(id)}" target="_blank" rel="noopener noreferrer">${escapeText(id)}</a><span>${escapeText(name || '')}</span></span>`).join('')}</p>`;
}

// Domains split into lanes so overlapping ones (e.g. two kinase-domain models) sit one above the other
function domainLanes(domains) {
    const lanes = [];
    [...domains].sort((a, b) => a.start - b.start).forEach(d => {
        const lane = lanes.find(l => l[l.length - 1].end < d.start);
        if (lane) lane.push(d); else lanes.push([d]);
    });
    return lanes;
}


// TRANSCRIPT ALIGNMENT
// Each Variants Evaluated row has an Align button that opens, beneath the row, a
// codon-aware alignment of the query transcript and the row's transcript, centred on the
// variant. The server aligns the proteins as the scoring does and lines the codons up under them.

const ALIGN_FLANK = 30;       // residues either side of the variant shown at first
const ALIGN_BLOCK = 20;       // alignment columns per printed line

function residueOf(aminoAcids) {
    const match = /^[A-Z*-]+(\d+)/.exec(aminoAcids || '');
    return match ? Number(match[1]) : null;
}

// Variant row for display with an Align button after its Transcript ID
function withAlignButton(displayed, row, org) {
    const key = org === 'human' ? row['AlleleID'] : row['Name'];
    const button = `<button type="button" class="alignToggle" aria-expanded="false"
        data-key="${escapeText(String(key ?? ''))}" data-tx="${escapeText(row['Transcript ID'] || '')}"
        data-residue="${residueOf(row['Amino Acids']) ?? ''}">Align</button>`;
    const out = {};
    Object.entries(displayed).forEach(([col, value]) => {
        out[col] = value;
        if (col === 'Transcript ID') out['Alignment'] = button;
    });
    if (!('Alignment' in out)) out['Alignment'] = button;
    return out;
}

function initAlignmentToggles(dataTable, org, queryProteins, scores) {
    const querySpecies = org;
    const targetSpecies = org === 'human' ? 'mouse' : 'human';
    const keyCol = org === 'human' ? 'AlleleID' : 'Name';

    $(dataTable.table().node()).find('tbody').on('click', 'button.alignToggle', function () {
        const row = dataTable.row($(this).closest('tr'));
        if (row.child.isShown()) {
            row.child.hide();
            this.setAttribute('aria-expanded', 'false');
            return;
        }
        this.setAttribute('aria-expanded', 'true');

        // Query transcripts this row was scored against, best tier first
        const { key, tx, residue } = this.dataset;
        const pairs = scores
            .filter(s => String(s[keyCol]) === key && s['Transcript ID'] === tx && s['Query Transcript ID'])
            .sort((a, b) => (a['Equivalency Tier'] ?? 99) - (b['Equivalency Tier'] ?? 99));
        const queryTxs = [...new Set(pairs.map(s => s['Query Transcript ID']))];
        if (queryTxs.length === 0) {
            row.child(`<div class="alignmentView"><p class="alignNote">No scored query transcript to align with.</p>${alignActions()}</div>`).show();
            initCloseButton($(row.child()).find('.alignmentView'), row, this);
            return;
        }

        const picker = queryTxs.length > 1
            ? `<label class="alignPicker">Query transcript <select>${queryTxs.map(q => `<option>${escapeText(q)}</option>`).join('')}</select></label>`
            : '';
        row.child(`<div class="alignmentView">${picker}<div class="alignBody"><p class="alignNote">Aligning…</p></div></div>`).show();
        const view = $(row.child()).find('.alignmentView');
        initCloseButton(view, row, this);

        const load = (queryTx, full) => {
            const queryRow = queryProteins.find(r => r['Transcript ID'] === queryTx);
            const params = new URLSearchParams({
                query_tx: queryTx, query_species: querySpecies, target_tx: tx, target_species: targetSpecies,
            });
            const queryResidue = queryRow ? residueOf(queryRow['Amino Acids']) : null;
            if (queryResidue) params.set('query_residue', queryResidue);
            if (residue) params.set('target_residue', residue);
            if (!full && queryResidue) params.set('flank', ALIGN_FLANK);

            view.find('.alignBody').html('<p class="alignNote">Aligning…</p>');
            fetch(`/alignment?${params}`)
                .then(response => response.json().then(body => ({ ok: response.ok, body })))
                .then(({ ok, body }) => {
                    if (!ok) throw new Error(body.error || 'The alignment could not be made.');
                    view.find('.alignBody').html(alignmentView(body, queryRow));
                    view.find('button.alignFull').on('click', () => load(queryTx, true));
                    view.find('button.alignWindow').on('click', () => load(queryTx, false));
                })
                .catch(error => view.find('.alignBody').html(`<p class="alignNote">${escapeText(error.message)}</p>${alignActions()}`));
        };

        view.find('select').on('change', function () { load(this.value, false); });
        load(queryTxs[0], false);
    });
}

const SPECIES_LABEL = { human: 'Human', mouse: 'Mouse' };

// Bottom row of an alignment or overview: optional buttons, then Close beside them
function alignActions(left = '') {
    return `<div class="alignActions">${left}<button type="button" class="alignClose">Close</button></div>`;
}

// Close closes the panel beneath the row, as clicking the row's own button again does
function initCloseButton(view, row, button) {
    view.on('click', 'button.alignClose', () => {
        row.child.hide();
        button.setAttribute('aria-expanded', 'false');
        button.focus();
    });
}

// Summary, gap list and the printed alignment
function alignmentView(data, queryRow) {
    const q = SPECIES_LABEL[data.query.species], t = SPECIES_LABEL[data.target.species];
    const st = data.stats;
    const cols = data.columns;

    const variantCol = cols.find(c => c.qp === data.query.residue);
    let lead = '';
    if (variantCol) {
        const aligned = variantCol.t === '-'
            ? `falls in a gap in the ${t.toLowerCase()} protein`
            : `aligns to ${t.toLowerCase()} ${variantCol.t}${variantCol.tp}`;
        const codons = variantCol.qc && variantCol.tc ? ` (codons ${variantCol.qc} and ${variantCol.tc})` : '';
        lead = `${q} ${variantCol.q}${variantCol.qp} ${aligned}${codons}. `;
        if (data.target.residue && variantCol.tp !== data.target.residue) {
            lead += `The ${t.toLowerCase()} variant is at residue ${data.target.residue}, which does not align to it. `;
        }
    }
    const change = queryRow && queryRow['Sequence Change'] ? ` The submitted change is ${escapeText(queryRow['Sequence Change'])}.` : '';

    const identity = `Proteins ${data.query.length} and ${data.target.length} residues long; `
        + `${st.protein_identity}% identical over ${st.aligned_residues} aligned residues`
        + (st.cds_identity !== null ? `, coding sequences ${st.cds_identity}% identical` : '')
        + (st.gap_columns ? `; ${st.gap_columns} gap ${st.gap_columns === 1 ? 'column' : 'columns'}.` : '; no gaps.');

    const gaps = data.gaps.map(g => {
        const lacking = g.missing_in === 'target' ? t : q;
        const having = g.missing_in === 'target' ? q : t;
        const range = g.from === g.to ? `residue ${g.from}` : `residues ${g.from}–${g.to}`;
        return `<li>${lacking} lacks ${g.length} ${g.length === 1 ? 'residue' : 'residues'} (${having.toLowerCase()} ${range})</li>`;
    });

    const notes = [];
    if (!data.query.codons || !data.target.codons) notes.push('Codons are not shown where a coding sequence is unavailable or does not match its protein.');
    if (data.windowed) notes.push('These proteins are long, so only the region around the variant was aligned.');

    const colors = domainColors([data.query.domains || [], data.target.domains || []]);

    const toggle = data.windowed ? ''
        : data.complete
            ? (data.query.residue ? '<button type="button" class="alignWindow">Show the region around the variant</button>' : '')
            : '<button type="button" class="alignFull">Show the full alignment</button>';

    return `
        <p class="alignSummary">${escapeText(data.query.transcript)} (${q}) vs ${escapeText(data.target.transcript)} (${t}). ${lead}${change}</p>
        <p class="alignStats">${identity}</p>
        ${gaps.length ? `<ul class="alignGaps">${gaps.join('')}</ul>` : ''}
        ${notes.map(n => `<p class="alignNote">${n}</p>`).join('')}
        ${domainLegend([data.query.domains || [], data.target.domains || []], colors)}
        <div class="alignScroll"><pre class="alignText">${alignmentText(cols, data, q, t, colors)}</pre></div>
        <p class="alignKey"><span class="alnVariant">Highlighted</span>: variant residue. <span class="alnDiff">Bold red</span>: base differs.
           | identical residue, . different residue.</p>
        ${alignActions(toggle)}`;
}

// Printed alignment, ALIGN_BLOCK columns per line group: query codons, query residues,
// match line, target residues, target codons
function alignmentText(cols, data, q, t, colors) {
    const showQ = data.query.codons, showT = data.target.codons;
    const width = Math.max(q.length, t.length) + 7;
    const label = (text) => text.padEnd(width);
    const codon = (mine, other) => !mine
        ? '---'
        : [...mine].map((base, i) => other && other[i] !== base ? `<b class="alnDiff">${base}</b>` : base).join('');
    const mark = (c, text) => {
        if (c.qp && c.qp === data.query.residue) return `<span class="alnVariant">${text}</span>`;
        if (c.tp && c.tp === data.target.residue) return `<span class="alnVariant alnTarget">${text}</span>`;
        return text;
    };
    // A line of bars under (or over) a sequence where its residues are in a domain, labelled where each domain starts
    const names = new Map([...(data.query.domains || []), ...(data.target.domains || [])].map(d => [d.id, d.name]));
    const domainLine = (block, key) => {
        if (!block.some(c => c[key])) return null;
        let out = '';
        for (let i = 0; i < block.length;) {
            const id = block[i][key];
            let j = i;
            while (j < block.length && block[j][key] === id) j++;
            const width = (j - i) * 4 - 1;
            if (!id) {
                out += ' '.repeat(width);
            } else {
                const name = `${id} ${names.get(id) || ''}`.trim();
                const fits = [name, id].find(l => l.length + 1 <= width);
                const text = (fits ? `${fits} ` : '').padEnd(width, '━').slice(0, width);
                out += `<span class="alnDomain" style="color:${colors.get(id)}" title="${escapeText(name)}">${escapeText(text)}</span>`;
            }
            if (j < block.length) out += ' ';
            i = j;
        }
        return out;
    };

    const first = (block, key) => block.find(c => c[key])?.[key] ?? '';
    const last = (block, key) => [...block].reverse().find(c => c[key])?.[key] ?? '';

    const groups = [];
    for (let i = 0; i < cols.length; i += ALIGN_BLOCK) {
        const block = cols.slice(i, i + ALIGN_BLOCK);
        const lines = [];
        const qDomains = domainLine(block, 'qd'), tDomains = domainLine(block, 'td');
        if (qDomains) lines.push(label('') + qDomains);
        if (showQ) lines.push(label('') + block.map(c => mark(c, codon(c.qc, c.tc))).join(' '));
        lines.push(label(`${q} ${first(block, 'qp')}`) + block.map(c => mark(c, ` ${c.q} `)).join(' ') + `  ${last(block, 'qp')}`);
        lines.push(label('') + block.map(c => c.q === '-' || c.t === '-' ? '   ' : c.q === c.t ? ' | ' : ' . ').join(' '));
        lines.push(label(`${t} ${first(block, 'tp')}`) + block.map(c => mark(c, ` ${c.t} `)).join(' ') + `  ${last(block, 'tp')}`);
        if (showT) lines.push(label('') + block.map(c => mark(c, codon(c.tc, c.qc))).join(' '));
        if (tDomains) lines.push(label('') + tDomains);
        groups.push(lines.join('\n'));
    }
    return groups.join('\n\n');
}


// TRANSCRIPT OVERVIEW
// Every transcript scored for an allele (or, for mouse queries, a human variant), drawn
// against the query transcript's residue numbering: solid bars where a transcript's
// residues align to the query protein, a thin line where it lacks the query's residues,
// ticks where it has residues the query lacks, and a mark at its own variant residue.

function initTranscriptOverview(dataTable, org, queryProteins, outputProteins, scores, mouseGene) {
    const keyCol = org === 'human' ? 'AlleleID' : 'Name';

    $(dataTable.table().node()).find('tbody').on('click', 'button.overviewToggle', function () {
        const row = dataTable.row($(this).closest('tr'));
        if (row.child.isShown()) {
            row.child.hide();
            this.setAttribute('aria-expanded', 'false');
            return;
        }
        this.setAttribute('aria-expanded', 'true');
        const { key, tx, qtx } = this.dataset;

        // The human transcript is the reference and the mouse transcripts compared with it
        // are drawn beneath. They are transcripts of the mouse gene, each marked with where the
        // allele's (human query) or submitted variant's (mouse query) change falls in it
        const human = org === 'human'
            ? { tx: qtx, variant: queryProteins.find(r => r['Transcript ID'] === qtx) }
            : { tx: tx, variant: outputProteins.find(r => String(r[keyCol]) === key && r['Transcript ID'] === tx) };
        const group = scores.filter(s => String(s[keyCol]) === key
            && (org === 'human' ? s['Query Transcript ID'] === qtx : s['Transcript ID'] === tx));
        const mouseTxCol = org === 'human' ? 'Transcript ID' : 'Query Transcript ID';
        const seen = new Set();
        const tracks = group.filter(s => !seen.has(s[mouseTxCol]) && seen.add(s[mouseTxCol])).map(s => {
            const mouseTx = s[mouseTxCol];
            const variant = org === 'human'
                ? outputProteins.find(v => String(v[keyCol]) === key && v['Transcript ID'] === mouseTx)
                : queryProteins.find(v => v['Transcript ID'] === mouseTx);
            return { tx: mouseTx, aminoAcids: variant ? variant['Amino Acids'] : null, score: s };
        });
        const current = org === 'human' ? tx : qtx;
        const humanChange = human.variant ? human.variant['Amino Acids'] : null;
        const label = org === 'human' ? `the ${(group[0] && group[0]['Allele Symbol']) || key} variant` : 'the submitted mouse variant';

        row.child('<div class="txOverview"><p class="alignNote">Aligning transcripts…</p></div>').show();
        const view = $(row.child()).find('.txOverview');
        initCloseButton(view, row, this);

        fetch('/transcript_overview', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                query_tx: human.tx, query_species: 'human', query_residue: residueOf(humanChange),
                target_species: 'mouse', targets: tracks.map(t => [t.tx, residueOf(t.aminoAcids)]),
            }),
        })
            .then(response => response.json().then(body => ({ ok: response.ok, body })))
            .then(({ ok, body }) => {
                if (!ok) throw new Error(body.error || 'The transcripts could not be aligned.');
                view.html(overviewView(body, tracks, humanChange, label, current, mouseGene) + alignActions());
            })
            .catch(error => view.html(`<p class="alignNote">${escapeText(error.message)}</p>${alignActions()}`));
    });
}

function axisStep(length) {
    return length <= 200 ? 50 : length <= 600 ? 100 : length <= 1500 ? 250 : length <= 4000 ? 500 : 5000;
}

function overviewView(data, targets, qChange, variantLabel, clickedTx, gene) {
    const q = SPECIES_LABEL[data.query.species], t = SPECIES_LABEL[data.target_species];
    const L = data.query.length;
    const pct = residue => `${(100 * (residue - 0.5) / L).toFixed(3)}%`;
    const span = (from, to) => `left:${(100 * (from - 1) / L).toFixed(3)}%;width:${(100 * (to - from + 1) / L).toFixed(3)}%`;
    const qResidue = data.query.residue;
    const qLabel = qChange ? qChange.replace(/[A-Z*]+$/, '') : (qResidue ? `residue ${qResidue}` : '');
    const guide = qResidue ? `<div class="ovGuide" style="left:${pct(qResidue)}"></div>` : '';

    const ticks = [];
    const step = axisStep(L);
    // Leave out a tick that would crowd the end label
    for (let r = step; r <= L - 0.3 * step; r += step) ticks.push(r);
    const axis = `<div class="ovRow ovAxisRow"><div class="ovLabel"></div><div class="ovTrack ovAxis">
        ${[1, ...ticks, L].map(r => `<span class="ovTick" style="left:${pct(r)}">${r}</span>`).join('')}</div></div>`;

    const colors = domainColors([data.query.domains || [], ...data.targets.map(d => d.domain_segments || [])]);
    const strips = segments => domainLanes(segments || []).map((lane, n) => lane.map(d =>
        `<div class="ovDomain" style="${span(d.start, d.end)};top:${23 + 6 * n}px;background:${colors.get(d.id)}"
            title="${escapeText(d.id)} ${escapeText(d.name || '')} (${q.toLowerCase()} residues ${d.start}–${d.end})"></div>`).join('')).join('');
    const trackHeight = segments => `height:${26 + 6 * Math.max(1, domainLanes(segments || []).length)}px`;

    const queryTrack = `<div class="ovRow ovQuery">
        <div class="ovLabel"><strong>${q} ${escapeText(data.query.transcript)}</strong><span>${L} aa${qChange ? ` · ${escapeText(qChange)}` : ''}</span></div>
        <div class="ovTrack" style="${trackHeight(data.query.domains)}">${guide}<div class="ovBlock" style="${span(1, L)}"></div>${strips(data.query.domains)}
            ${qResidue ? `<div class="ovMark ovOn" style="left:${pct(qResidue)}"><span>${escapeText(qLabel)}</span></div>` : ''}</div></div>`;

    const byTx = new Map(data.targets.map(d => [d.transcript, d]));
    const rows = targets.map(({ tx, score, aminoAcids }) => {
        const d = byTx.get(tx);
        const txLabel = escapeText(tx);
        const current = tx === clickedTx ? ' ovCurrent' : '';
        if (!d || d.error) {
            return `<div class="ovRow${current}"><div class="ovLabel"><strong>${txLabel}</strong><span>${escapeText(d?.error || 'Not aligned')}</span></div><div class="ovTrack"></div></div>`;
        }
        const first = d.blocks[0], last = d.blocks[d.blocks.length - 1];
        const baseline = first ? `<div class="ovLine" style="${span(first[0], last[1])}"></div>` : '';
        const blocks = d.blocks.map(b => `<div class="ovBlock" style="${span(b[0], b[1])}" title="${t} ${b[2]}–${b[3]} align to ${q.toLowerCase()} ${b[0]}–${b[1]}"></div>`).join('');
        const inserts = d.insertions.map(i => `<div class="ovInsert" style="left:${pct(Math.max(i.after, 0.5) + 0.5)}"
            title="${i.count} ${t.toLowerCase()} ${i.count === 1 ? 'residue' : 'residues'} (${i.from}–${i.to}) with no ${q.toLowerCase()} counterpart, after ${q.toLowerCase()} residue ${i.after}">${i.count >= 3 ? `+${i.count}` : ''}</div>`).join('');

        const change = aminoAcids ? aminoAcids.replace(/[A-Z*]+$/, '') : (d.residue ? `residue ${d.residue}` : '');
        let mark = '', note = '';
        if (d.residue && d.residue_at_query) {
            const on = d.residue_at_query === qResidue;
            mark = `<div class="ovMark ${on ? 'ovOn' : 'ovOff'}" style="left:${pct(d.residue_at_query)}"><span>${escapeText(change)}</span></div>`;
            if (!on) note = ` · aligns to ${q.toLowerCase()} ${d.residue_at_query}`;
        } else if (d.residue) {
            note = ' · variant residue has no counterpart';
        }
        const total = score['Equivalency Tier'] !== null && score['Equivalency Tier'] !== undefined ? ` · Tier ${score['Equivalency Tier']}` : '';
        const details = `${d.length} aa${aminoAcids ? ` · ${escapeText(aminoAcids)}` : ''}${note}${total}`;
        return `<div class="ovRow${current}">
            <div class="ovLabel"><strong>${txLabel}</strong><span>${details}</span></div>
            <div class="ovTrack" style="${trackHeight(d.domain_segments)}">${guide}${baseline}${blocks}${strips(d.domain_segments)}${inserts}${mark}</div></div>`;
    }).join('');

    const offCount = data.targets.filter(d => !d.error && d.residue && d.residue_at_query !== qResidue).length;
    // The tracks are the gene's transcripts (not the allele's); each shows where the one variant falls in it
    const n = targets.length;
    const geneName = gene ? `${t.toLowerCase()} <em>${escapeText(gene)}</em> gene` : `${t.toLowerCase()} gene`;
    const one = n === 1;
    const summary = `Alignment of ${n} ${one ? 'transcript' : 'transcripts'} of the ${geneName}, `
        + `${one ? '' : 'each '}aligned to ${q.toLowerCase()} ${escapeText(data.query.transcript)} and drawn on its residue numbering (1–${L}). `
        + `The mark on ${one ? 'the' : 'each'} transcript shows where ${variantLabel} falls in that transcript. `
        + (targets.some(x => x.tx === clickedTx) ? 'The highlighted row is the transcript of the table row selected. ' : '')
        + (qResidue ? `The dashed line marks ${q.toLowerCase()} ${escapeText(qLabel)}. ` : '')
        + (offCount ? `In ${offCount} ${offCount === 1 ? 'transcript' : 'transcripts'} the variant residue aligns somewhere else.`
            : one ? 'In the transcript the variant residue aligns with it, whatever its residue number in that transcript.'
                  : 'In every transcript the variant residue aligns with it, whatever its residue number in that transcript.');

    return `
        <p class="alignSummary">${summary}</p>
        ${domainLegend([data.query.domains || [], ...data.targets.map(d => d.domain_segments || [])], colors)}
        <div class="ovChart">${axis}${queryTrack}${rows}</div>
        <p class="alignKey"><span class="ovKeyBlock"></span> aligned residues
            <span class="ovKeyLine"></span> ${q.toLowerCase()} residues this transcript lacks
            <span class="ovKeyInsert"></span> residues only in this transcript
            <span class="ovKeyMark ovOn"></span> variant residue aligned with the ${q.toLowerCase()} variant
            <span class="ovKeyMark ovOff"></span> variant residue elsewhere.
            Coloured strips beneath each bar: Pfam domains, placed by the alignment on ${q.toLowerCase()} numbering.
            Hover over a bar or strip for its residue numbers.</p>`;
}


// Mouse Allele Coverage table and its summary line
// PHENOTYPE SIMILARITY (human queries)
// One row per evaluated allele, in the order of the scoring table, with each phenotype
// score and a button that opens the term mapping behind it in a row beneath the allele.

function phenotypeSimilarityRows(name) {
    return (window.fullJSON.phenotype_similarity || []).filter(row => row['Input Name'] === name);
}

function phenotypeSimilaritySection(name, gene, sortedScores) {
    const rows = phenotypeSimilarityRows(name);
    if (rows.length === 0) return '';

    // Same allele order as the scoring table; symbols come from the score rows
    const order = new Map();
    const symbols = new Map();
    sortedScores.forEach(r => {
        if (!order.has(r['AlleleID'])) order.set(r['AlleleID'], order.size);
        symbols.set(r['AlleleID'], r['Allele Symbol']);
    });
    rows.sort((a, b) => (order.get(a['AlleleID']) ?? Infinity) - (order.get(b['AlleleID']) ?? Infinity));

    const score = v => (v === null || v === undefined) ? '<span class="naScore">—</span>' : Number(v).toFixed(2);
    const button = (row, kind) =>
        `<button type="button" class="phenotypeToggle mappingToggle" data-allele="${row['AlleleID']}" data-kind="${kind}" aria-expanded="false">Phenotype term mapping</button>`;

    const body = rows.map(row => {
        const aid = row['AlleleID'];
        const symbol = symbols.get(aid) || aid;
        return `<tr>
            <td><a href="https://www.informatics.jax.org/allele/${aid}" target="_blank" rel="noopener noreferrer">${symbol}</a><br><span class="alleleId">${aid}</span></td>
            <td>${escapeText(row['Matched Disease(s)'] || '')}</td>
            <td>${row['MP Terms Mapped']} of ${row['MP Terms Total']}</td>
            <td class="scoreCell"><span class="scoreValue">${score(row['Phenotype Score (variant)'])}</span>${button(row, 'variant')}</td>
            <td class="scoreCell"><span class="scoreValue">${score(row['Phenotype Score (gene)'])}</span>${button(row, 'gene')}</td>
        </tr>`;
    }).join('');

    return `
        <div class='tableWrapper'>
            ${infoHeading('similarityInfo', 'Phenotype Similarity', INFO.phenotypeSimilarity())}
            <p class="coverageSummary">How closely each evaluated allele's MGI phenotypes match the human phenotypes.
            Mouse (MP) terms are translated to human (HPO) terms, and each translated term is matched to its closest
            human term by semantic similarity (1 for the identical term, less for a broader, narrower or related term).
            The score, from 0 to 1, is how much of the allele's phenotype the human terms account for.
            <strong>Variant</strong> compares with each disease of the human variant and reports the best;
            <strong>Gene</strong> compares with all HPO terms of ${escapeText(gene)}.
            Open <em>Phenotype term mapping</em> to see every term behind a score.</p>
            <table id="phenotypeSimilarityTable" class="table.dataTable">
                <thead><tr>
                    <th>Mouse Allele</th>
                    <th>Matched Disease(s)</th>
                    <th>MP Terms Mapped to HPO</th>
                    <th>Variant Phenotype Score</th>
                    <th>Gene Phenotype Score</th>
                </tr></thead>
                <tbody>${body}</tbody>
            </table>
        </div>`;
}

function initPhenotypeSimilarityTable(name) {
    const table = $('#phenotypeSimilarityTable');
    if (table.length === 0) return;
    const dataTable = table.DataTable({ ordering: false });
    const byAllele = new Map(phenotypeSimilarityRows(name).map(row => [row['AlleleID'], row]));

    // One mapping open per allele: the other score's button switches it, the same button closes it
    $('#phenotypeSimilarityTable tbody').on('click', 'button.mappingToggle', function () {
        const tr = $(this).closest('tr');
        const row = dataTable.row(tr);
        const kind = this.dataset.kind;
        const open = row.child.isShown() && tr.attr('data-open') === kind;

        tr.find('button.mappingToggle').attr('aria-expanded', 'false');
        if (open) {
            row.child.hide();
            tr.removeAttr('data-open');
        } else {
            row.child(mappingView(byAllele.get(this.dataset.allele), kind)).show();
            tr.attr('data-open', kind);
            this.setAttribute('aria-expanded', 'true');
        }
    });
}

// The terms behind one phenotype score. Each mouse term is translated to HPO, then compared
// with its closest human term by semantic similarity (1 = the identical term)
function mappingView(row, kind) {
    const m = row[kind === 'variant' ? 'Variant Mapping' : 'Gene Mapping'];
    const link = (url, label, id) => `<a href="${url}" target="_blank" rel="noopener noreferrer">${escapeText(label || id)}</a> <span class="termId">${id}</span>`;
    const mpTerm = (id, label) => link(`https://www.informatics.jax.org/vocab/mp_ontology/${id}`, label, id);
    const hpTerm = (id, label) => link(`https://hpo.jax.org/browse/term/${id}`, label, id);
    const byName = (a, b) => String(a.mp_name || a.mp).localeCompare(String(b.mp_name || b.mp), undefined, {sensitivity: 'base'});

    // Terms outside the score: no HPO translation, or general terms that are not scored
    const listOf = (items, summary) => items && items.length
        ? `<details class="unmappedTerms"><summary>${summary(items.length)}</summary>
            ${termList(items.map(u => ({id: u.mp, name: u.mp_name || ''})))}</details>`
        : '';
    const notScored = listOf(m && m.not_scored, n => `${n} general term${n === 1 ? '' : 's'} not scored (top-level system terms and "no abnormal phenotype detected")`);
    const unmapped = listOf(m && m.unmapped, n => `${n} mouse phenotype${n === 1 ? '' : 's'} with no HPO translation (not part of the score)`);
    const normal = listOf(m && m.normal, n => `${n} system${n === 1 ? '' : 's'} or phenotype${n === 1 ? '' : 's'} MGI recorded as <strong>normal</strong> (examined, no abnormality found; not part of the score)`);
    const outside = `${normal}${unmapped}${notScored}`;

    if (!m || m.total_mp_count === 0) {
        const why = m && m.not_scored && m.not_scored.length
            ? 'This allele has only general MGI phenotype terms, which are not scored, so there is nothing to compare.'
            : m && m.normal && m.normal.length
                ? 'MGI records only normal findings for this allele, so there is no abnormal phenotype to compare.'
                : 'This allele has no MGI phenotype annotations, so there is nothing to compare.';
        return `<div class="mappingView"><p class="mappingNote">${why}</p>${outside}</div>`;
    }
    if (!m.sections || m.sections.length === 0) {
        const why = kind === 'variant'
            ? 'none of the variant\'s diseases has HPO phenotype annotations, or none of the allele\'s phenotypes translates to HPO.'
            : 'the human gene has no HPO annotations, or none of the allele\'s phenotypes translates to HPO.';
        return `<div class="mappingView"><p class="mappingNote">No ${kind} phenotype score: ${why}</p>${outside}</div>`;
    }

    const key = `<p class="mappingKey"><strong>Translation</strong>: how closely the mouse term corresponds to an HPO term
        (exact 1, narrow or broad 0.75, close 0.5, related 0.25), with the weight it carries.
        <strong>Similarity</strong>: how close that HPO term is to the nearest ${kind === 'variant' ? "disease" : "gene"} HPO term,
        from the HPO hierarchy and how specific the terms are: 1 for the identical term, high for a parent or child,
        0 when they share no more than an organ system.</p>`;

    const linkText = p => ({
        'identical': 'identical',
        'mouse term more specific': 'mouse term is more specific',
        'mouse term more general': 'mouse term is more general',
        'shared ancestor': p.mica ? `related through ${escapeText(p.mica_name || p.mica)}` : 'related',
    }[p.link] || '');
    const level = s => s >= 0.999 ? 'simFull' : s >= 0.5 ? 'simHigh' : s > 0 ? 'simLow' : 'simNone';

    const section = sec => {
        const pairs = [...sec.pairs].sort((a, b) => (b.similarity - a.similarity) || byName(a, b));
        const heading = kind === 'variant'
            ? `Compared with the HPO terms of <strong>${escapeText(sec.reference)}</strong>`
            : 'Compared with all HPO terms of the human gene';
        const rows = pairs.map(p => `
            <tr class="${level(p.similarity)}">
                <td>${mpTerm(p.mp, p.mp_name)}</td>
                <td class="relation">${(p.relation || '').replace('Match', '') || 'mapped'} <span class="termId">${p.weight}</span></td>
                <td>${hpTerm(p.hp, p.hp_name)}${p.counted ? '' : '<span class="countedNote">counted under another mouse term</span>'}</td>
                <td>${p.human_hp ? hpTerm(p.human_hp, p.human_hp_name) : '<span class="naScore">none related</span>'}</td>
                <td class="simCell"><span class="simBar"><span style="width:${Math.round(p.similarity * 100)}%"></span></span>
                    <span class="simValue">${Number(p.similarity).toFixed(2)}</span>
                    ${p.human_hp ? `<span class="linkNote">${linkText(p)}</span>` : ''}</td>
            </tr>`).join('');
        return `<div class="mappingSection">
            <p class="mappingNote">${heading}: <strong>score ${Number(sec.score).toFixed(2)}</strong>
                = weight × similarity summed (${sec.matched_weight}) ÷ weight summed (${sec.total_weight}).
                Each HPO term counts once, at its strongest translation.</p>
            <div class="mappingScroll"><table class="mappingTable">
                <thead><tr><th>Mouse phenotype (MP)</th><th>Translation, weight</th><th>Translated to (HPO)</th>
                    <th>Closest ${kind === 'variant' ? 'disease' : 'gene'} term (HPO)</th><th>Similarity</th></tr></thead>
                <tbody>${rows}</tbody>
            </table></div></div>`;
    };

    const tie = m.sections.length > 1
        ? `<p class="mappingNote">${m.sections.length} diseases tie for the best score; each is shown with its own comparison.</p>`
        : '';
    return `<div class="mappingView">${tie}${key}${m.sections.map(section).join('')}${outside}</div>`;
}


function coverageSection(alleleCoverage, gene, org) {
    const total = alleleCoverage.length;
    let summary;

    if (total === 0) {
        summary = 'No mouse alleles found in MGI.';
    } else if (org === 'human') {
        const evaluated = alleleCoverage.filter(r => r.Evaluated).length;
        const withDisease = alleleCoverage.filter(r => r.Diseases && r.Diseases.length > 0).length;
        summary = `VarLift evaluated <strong>${evaluated} of ${total}</strong> ${gene} alleles in MGI (${withDisease} with a disease association). ${total - evaluated} ${total - evaluated === 1 ? 'allele exists' : 'alleles exist'} but could not be evaluated for equivalency (see table for reasons).`;
    } else {
        const submitted = alleleCoverage.filter(r => r.Evaluated).map(r => r.AlleleSymbol);
        const withDisease = alleleCoverage.filter(r => r.Diseases && r.Diseases.length > 0).length;
        // Alleles of the gene's cluster (e.g. Pcdhg) and of overlapping genes are listed too
        const own = alleleCoverage.filter(r => !r.Marker || r.Marker === gene).length;
        const otherMarkers = [...new Set(alleleCoverage.filter(r => r.Marker && r.Marker !== gene).map(r => r.Marker))];
        const noCoordinates = alleleCoverage.filter(r => !r.Consequence).length;
        const plural = n => n === 1 ? 'allele' : 'alleles';
        const perMarker = otherMarkers.map(m => `${m} (${alleleCoverage.filter(r => r.Marker === m).length})`);
        const counts = otherMarkers.length === 0
            ? `MGI has <strong>${total}</strong> ${gene} ${plural(total)} (${withDisease} with a disease association).`
            : `MGI has ${own === 0 ? 'no' : `<strong>${own}</strong>`} ${gene} ${plural(own)}, and <strong>${total - own}</strong> ${plural(total - own)} of other markers at this position: ${perMarker.join(', ')}. ${withDisease === 0 ? 'None has' : `${withDisease} of the ${total} ${withDisease === 1 ? 'has' : 'have'}`} a disease association.`;
        const coordinates = noCoordinates === 0 ? ''
            : noCoordinates === total
                ? ` ${total === 1 ? 'It has' : 'None of them have'} variant coordinates in MGI (e.g. large deletions or knockouts).`
                : ` ${noCoordinates} of them have no variant coordinates in MGI (e.g. large deletions or knockouts).`;
        summary = `${counts} ${submitted.length > 0
            ? `The submitted variant corresponds to ${submitted.join(', ')}.`
            : `The submitted variant's coordinates do not match any of them.${coordinates}`}`;
    }

    return `
        <div class='tableWrapper'>
            ${infoHeading('coverageInfo', 'Mouse Allele Coverage', INFO.coverage(org))}
            <p class="coverageSummary">${summary}</p>
            ${jsonToHTMLTable(alleleCoverage, "alleleCoverageTable")}
        </div>`;
}


// Make the coverage table sortable and open each allele's phenotypes beneath its row
function initCoverageTable(phenotypesByAllele) {
    const table = $('#alleleCoverageTable');
    if (table.length === 0) return;

    const headers = table.find('thead th').map((i, th) => $(th).text().trim()).get();

    // Marker matters only when alleles of a cluster or overlapping gene are listed with the gene's own
    const markerCol = headers.indexOf('Marker');
    const showMarker = markerCol !== -1 &&
        [...phenotypesByAllele.values()].some(row => row['Marker'] && row['Marker'] !== row['Gene Symbol']);

    const alleleCoverageTable = table.DataTable({
        ordering: true,
        order: [[headers.indexOf('Evaluated'), 'desc'],
                ...(showMarker ? [[markerCol, 'asc']] : []),
                [headers.indexOf('AlleleSymbol'), 'asc']],
        columnDefs: [{ targets: headers.indexOf('Phenotypes'), orderable: false }]
    });
    if (markerCol !== -1 && !showMarker) alleleCoverageTable.column(markerCol).visible(false);

    // Hide Gene Symbol (same for all rows of a single-gene variant query), and
    // AlleleID, which is shown beneath the allele symbol
    alleleCoverageTable.column(headers.indexOf('Gene Symbol')).visible(false);
    alleleCoverageTable.column(headers.indexOf('AlleleID')).visible(false);
    // Reason opens from the Evaluated cell of alleles that weren't evaluated
    if (headers.includes('Reason')) alleleCoverageTable.column(headers.indexOf('Reason')).visible(false);
    // Normal findings open with the abnormal ones in the Phenotypes dropdown
    if (headers.includes('Normal Phenotypes')) alleleCoverageTable.column(headers.indexOf('Normal Phenotypes')).visible(false);

    initPhenotypeToggles('alleleCoverageTable', alleleCoverageTable, phenotypesByAllele);
}


// Show what is known about a variant VarLift could not evaluate
function renderUnevaluated(name) {
    const org = document.getElementById('queryOrganism').value;
    const queryOrganism = org === 'human' ? 'Human' : 'Mouse';
    const isEmpty = v => v === null || v === undefined || v === '' || (Array.isArray(v) && v.length === 0);
    const matches = window.fullJSON.unevaluated_variants.filter(row => row['Input Name'] === name);
    const hasAlleles = matches.some(row => row['AlleleID']);

    // Mouse variants inside a gene also get that gene's coverage table
    const gene = matches[0]['Gene Symbol'];
    const alleleCoverage = org === 'mouse' && gene ? markEvaluated(window.fullJSON.allele_coverage, gene, matches) : [];

    // Allele phenotypes go in their own table so long lists don't stretch the rows
    // (alleles already listed in the coverage table are left out)
    const listed = new Set(alleleCoverage.map(row => row['AlleleID']));
    const phenotypes = matches
        .filter(row => (!isEmpty(row['Phenotypes']) || !isEmpty(row['Normal Phenotypes'])) && !listed.has(row['AlleleID']))
        .map(row => ({AlleleID: row['AlleleID'], AlleleSymbol: row['AlleleSymbol'], Phenotypes: row['Phenotypes'] || [],
                      'Normal Phenotypes': row['Normal Phenotypes'] || []}));

    // Biotype is a gene property, shown only in the Genes table
    const rows = matches.map(({Phenotypes, 'Normal Phenotypes': _normal, Biotype, ...row}) => withLabel(name, displayFields({
        ...row,
        'Associated Diseases': isEmpty(row['Associated Diseases']) ? 'None reported' : row['Associated Diseases']
    })));

    results.innerHTML = `
        <div class='tableWrapper'>
            ${infoHeading('unevaluatedInfo', `Submitted ${queryOrganism} Variant (not evaluated)`, INFO.unevaluated(org))}
            <p class="coverageSummary">${rows[0]['Reason']}. No homologous variants or similarity scores are reported. ${org === 'mouse'
                ? (hasAlleles ? 'The MGI allele(s) at this position are listed below.' : 'No MGI allele matches this position and allele.')
                : 'ClinVar disease associations for this variant are listed below.'}</p>
            ${jsonToHTMLTable(rows, 'unevaluatedTable')}
        </div>

        ${phenotypes.length === 0 ? '' : `
        <div class='tableWrapper'>
            ${infoHeading('unevaluatedPhenotypesInfo', `${queryOrganism} Phenotypes`, INFO.unevaluatedPhenotypes())}
            ${jsonToHTMLTable(phenotypes, 'unevaluatedPhenotypeTable')}
        </div>`}

        ${alleleCoverage.length === 0 ? '' : coverageSection(alleleCoverage, gene, org)}
    `;

    const table = $('#unevaluatedTable').DataTable({
        dom: 't',
        pageLength: 20,
        ordering: false
    });

    // Input Name and Reason are already shown above the table, and columns
    // with no values (e.g. Exon for an intronic variant) are left out
    const hidden = ['Input Name', 'Reason', 'MONDO'];

    table.columns().every(function () {
        const header = $(this.header()).text().trim();
        const unused = rows.every(row => isEmpty(row[header]));
        if (hidden.includes(header) || unused) this.visible(false);
    });

    $('#unevaluatedPhenotypeTable').DataTable({
        dom: 't',
        pageLength: 20,
        ordering: false
    });

    initCoverageTable(new Map(alleleCoverage.map(row => [row['AlleleID'], row])));
}


// An allele's abnormal phenotypes, then the systems or phenotypes MGI recorded as normal
function phenotypeLists(phenotypes, normal) {
    const abnormal = phenotypes && phenotypes.length ? phenotypeLinks(phenotypes) : '';
    const normalList = normal && normal.length
        ? `<div class="normalPhenotypes"><span class="normalHeading">Normal: examined, no abnormality found</span>${phenotypeLinks(normal)}</div>`
        : '';
    return abnormal + normalList;
}


// "ID,label" phenotype strings as a term list (see termList)
function phenotypeLinks(phenotypes) {
    return termList(phenotypes.map(v => {
        // Split on the first comma only; labels can contain commas
        const text = String(v);
        const comma = text.indexOf(',');
        return comma < 0 ? {id: text.trim(), name: ''}
            : {id: text.slice(0, comma).trim(), name: text.slice(comma + 1).trim()};
    }));
}


// Turn an array of objects into an HTML table
function jsonToHTMLTable(data, tableId) {
    if (!data || data.length === 0) return `<p>No data</p>`;

    const org = document.getElementById('queryOrganism').value;
    const columns = Object.keys(data[0]);

    let html = `<table id="${tableId}" class="table.dataTable">`;
    html += `<thead><tr>`;

    // Coverage table: the allele's disease models are human diseases. Human Variants
    // Evaluated: ClinVar's submitted conditions
    const label = col => (tableId === 'alleleCoverageTable' && col === 'Diseases') ? 'Human Diseases Modeled'
        : (tableId === 'outputVariantTable' && org === 'mouse' && col === 'Associated Diseases') ? 'Associated Conditions'
        : (tableId === 'scoresTable' && col === 'Equivalency Tier') ? 'Equivalency Tier (1–5)'
        : col;

    columns.forEach(col => {
        html += `<th>${label(col)}</th>`;
    });

    html += `</tr></thead><tbody>`;

    data.forEach(row => {
        html += `<tr>`;
        columns.forEach(col => {
            html += `<td>${formatCell(row, col, tableId, org)}</td>`;
        });
        html += `</tr>`;
    });

    html += `</tbody></table>`;

    return html;
}


// Link for a MONDO disease ID, to its Monarch Initiative page
function mondoLink(id) {
    return `<a href="https://monarchinitiative.org/${id}" target="_blank" rel="noopener noreferrer">${id}</a>`;
}


// Link for an HPO term ID
// HP terms as a list sorted by name: the term name links to HPO, its ID follows in small grey text.
// Names come from the server with the results; a term without one shows its ID only.
function hpTermList(ids) {
    const names = (window.fullJSON && window.fullJSON.hp_names) || {};
    return termList([...new Set(ids)].map(id => ({id, name: names[id] || ''})));
}


// A list of ontology terms (HP or MP), one per row: the name, linked, with its ID in grey
// after it, sorted by name (terms without a name last)
function termList(terms) {
    const seen = new Set();
    const items = terms
        .filter(({id}) => !seen.has(id) && seen.add(id))
        .sort((a, b) => (a.name === '') - (b.name === '')
            || a.name.localeCompare(b.name, undefined, {sensitivity: 'base'})
            || a.id.localeCompare(b.id))
        .map(({id, name}) => {
            const url = id.startsWith('MP:')
                ? `https://www.informatics.jax.org/vocab/mp_ontology/${id}`
                : `https://hpo.jax.org/browse/term/${id}`;
            const label = name ? escapeText(name) : id;
            const idTag = name ? ` <span class="hpId">${id}</span>` : '';
            return `<li><a href="${url}" target="_blank" rel="noopener noreferrer">${label}</a>${idTag}</li>`;
        });
    return `<ul class="hpTermList">${items.join('')}</ul>`;
}


// Escape text from data files before it goes into HTML
function escapeText(text) {
    return String(text).replace(/[&<>"']/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}[c]));
}


// Link for a gene accession: MGI IDs to MGI, HGNC IDs to HGNC
function accessionLink(id) {
    if (!id) return '';
    const url = id.startsWith('MGI:') ? `https://www.informatics.jax.org/marker/${id}`
              : id.startsWith('HGNC:') ? `https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/${id}`
              : null;
    return url ? `<a href="${url}" target="_blank" rel="noopener noreferrer">${id}</a>` : id;
}


// Gene row for display: the accession goes in parentheses after the Gene Symbol,
// Description becomes Gene Name, and Chromosome, Start, End and Strand become one
// Genome Location field (e.g. "Chr 2:105,499,234-105,530,952 (+)")
function geneFields(row) {
    const out = {};
    Object.entries(row).forEach(([col, value]) => {
        if (col === 'Gene Symbol') {
            out[col] = row['Accession'] ? `${value} (${accessionLink(row['Accession'])})` : value;
        } else if (col === 'Description') {
            out['Gene Name'] = value;
        } else if (col === 'Chromosome') {
            const start = Number(row['Start']).toLocaleString('en-US');
            const end = Number(row['End']).toLocaleString('en-US');
            out['Genome Location'] = `Chr ${value}:${start}-${end} (${row['Strand']})`;
        } else if (col === 'Ensembl_ID') {
            // Ensembl and NCBI Gene IDs share one line, e.g. "ENSG00000169057, Gene ID: 4204"
            const ncbi = row['NCBI_Gene_ID'];
            const ncbiRef = ncbi ? `Gene ID: <a href="https://www.ncbi.nlm.nih.gov/gene/${ncbi}" target="_blank" rel="noopener noreferrer">${ncbi}</a>` : '';
            const ensemblRef = value ? `<a href="https://www.ensembl.org/id/${value}" target="_blank" rel="noopener noreferrer">${value}</a>` : '';
            out['Cross References'] = [ensemblRef, ncbiRef].filter(Boolean).join(', ');
        } else if (!['Start', 'End', 'Strand', 'Accession', 'NCBI_Gene_ID'].includes(col)) {
            out[col] = value;
        }
    });
    return out;
}


// Variant row for display: a Genome Location field first, taken from the input name
// chromosome:start-end:ref/alt (e.g. "Chr 2:105,522,966 G>T"), and on mouse
// variants the linked allele ID in parentheses after the AlleleSymbol
// Information button for a section heading; it opens the window with the same id
function infoButton(id, title) {
    return `<button type="button" class="infoButton" data-dialog="${id}" aria-label="About the ${title} table"
        title="About this table">i</button>`;
}


// Closable information window explaining a results table
function infoDialog(id, title, body) {
    return `
    <dialog id="${id}" class="infoDialog" aria-labelledby="${id}Title">
        <div class="infoDialogHeader">
            <h3 id="${id}Title">About the ${title} table</h3>
            <button type="button" class="infoClose" aria-label="Close">&times;</button>
        </div>
        <div class="infoDialogBody">${body}</div>
    </dialog>`;
}


// Heading with its information button, followed by the window
function infoHeading(id, title, body, heading = title) {
    return `<h2>${heading} ${infoButton(id, title)}</h2>${infoDialog(id, title, body)}`;
}


// Information buttons open their window; the close button, Escape or a click outside closes it.
// Delegated, so it covers windows added each time the results are drawn
document.addEventListener('click', event => {
    const button = event.target.closest('.infoButton');
    if (button) {
        document.getElementById(button.dataset.dialog)?.showModal();
        return;
    }
    if (event.target.closest('.infoClose')) {
        event.target.closest('dialog').close();
        return;
    }
    if (event.target instanceof HTMLDialogElement && event.target.classList.contains('infoDialog')) {
        event.target.close();
    }
});


const INFO = {
    submitted: (org) => {
        const human = org === 'human';
        return `
            <p>The variant you submitted, as annotated by Ensembl VEP on the ${human ? 'MANE Select human transcript (GRCh38)' : 'canonical mouse transcript (GRCm39)'}.
               If the variant falls in several transcripts, each one gets its own column.</p>
            <dl>
                <dt>Label</dt><dd>The label you entered after the variant, if any. It is not used in the analysis.</dd>
                <dt>Genome Location</dt><dd>Chromosome, position and change, as submitted.</dd>
                <dt>Transcript ID, Exon / Intron</dt><dd>The transcript VEP annotated and where the variant falls in it.</dd>
                <dt>Pfam Domain Name</dt><dd>The protein domain the variant is in, with its linked Pfam ID.</dd>
                <dt>PolyPhen</dt><dd>The predicted effect of a missense change. Shown for reference; not used in scoring.</dd>
                <dt>Molecular Consequence</dt><dd>VEP's consequence. It decides the scoring class (Coding or Splice/Intronic).</dd>
                <dt>Sequence Change, Amino Acid Change</dt><dd>The codon and amino-acid change on this transcript.</dd>
                ${human ? `
                <dt>Associated Conditions</dt><dd>The variant's conditions in ClinVar, with their MONDO IDs.</dd>
                <dt>HP_IDs_variant</dt><dd>HPO terms for this variant, from ClinVar and the HPO annotations of its diseases. Used for the variant phenotype score.</dd>
                <dt>HP_IDs_gene</dt><dd>Every HPO term annotated to the gene. Used for the gene phenotype score.</dd>` : `
                <dt>AlleleSymbol</dt><dd>MGI alleles at exactly this position and change, linked to MGI.</dd>
                <dt>Associated Diseases</dt><dd>The diseases those alleles model, with their MONDO IDs.</dd>`}
            </dl>`;
    },

    genes: (org) => `
        <p>The gene the submitted variant is in and its ${org === 'human' ? 'mouse' : 'human'} ortholog, side by side.
           The ortholog comes from MGI's human–mouse homology report.</p>
        <dl>
            <dt>Gene Symbol</dt><dd>Linked to HGNC (human) or MGI (mouse).</dd>
            <dt>Gene Name, Biotype</dt><dd>The gene's full name and type (for example, protein_coding).</dd>
            <dt>Genome Location</dt><dd>Chromosome, start–end and strand (GRCh38 for human, GRCm39 for mouse).</dd>
            <dt>Cross References</dt><dd>The Ensembl gene ID and NCBI Gene ID, each linked to its database.</dd>
        </dl>`,

    variantsEvaluated: (org) => org === 'human' ? `
        <p>Every MGI allele of the mouse ortholog that has genomic coordinates and a scorable consequence,
           annotated by VEP on <strong>every</strong> mouse transcript it falls in. An allele therefore has one row per transcript.
           Rows are in the same order as the Variant Equivalency Scoring table.</p>
        <dl>
            <dt>AlleleID, AlleleSymbol</dt><dd>The MGI allele, linked to its MGI page.</dd>
            <dt>Align</dt><dd>Opens a codon-by-codon alignment of this transcript with the human transcript around the variant, with Pfam domains marked.</dd>
            <dt>Transcript ID, Exon / Intron</dt><dd>The mouse transcript and where the allele falls in it.</dd>
            <dt>Amino Acid Change</dt><dd>The change on that transcript. Residue numbers can differ between transcripts and between species.</dd>
            <dt>Associated Diseases</dt><dd>The diseases MGI curates the allele as a model of, with their MONDO IDs.</dd>
        </dl>
        <p>Alleles without coordinates, or whose consequence can't be scored, are not listed here; see Mouse Allele Coverage.</p>` : `
        <p>Candidate human variants from ClinVar in the human ortholog. Coding variants are candidates when they make the
           <strong>same amino-acid change</strong> (for example, any R&rarr;C); splice-site and intronic variants when they are
           at the <strong>same intron offset</strong> (for example, +1). Whether a candidate is at the equivalent position is
           tested in the Variant Equivalency Scoring table. Each candidate has one row per human transcript.</p>
        <dl>
            <dt>Name</dt><dd>The ClinVar variant name (HGVS).</dd>
            <dt>Align</dt><dd>Opens a codon-by-codon alignment of this transcript with the mouse transcript around the variant, with Pfam domains marked.</dd>
            <dt>Transcript ID, Exon / Intron</dt><dd>The human transcript and where the variant falls in it.</dd>
            <dt>Amino Acid Change</dt><dd>The change on that transcript.</dd>
            <dt>Associated Conditions</dt><dd>The variant's conditions in ClinVar, with their MONDO IDs.</dd>
        </dl>`,

    scoring: (org) => {
        const human = org === 'human';
        return `
            <p>Each row compares the submitted ${human ? 'human' : 'mouse'} variant with one ${human ? 'mouse allele' : 'human ClinVar variant'}
               on one ${human ? 'mouse' : 'human'} transcript. The tier measures how closely the two variants make the same
               molecular change; it compares annotations and does not predict pathogenicity.</p>

            <h4>Columns</h4>
            <dl>
                <dt>${human ? 'Allele Symbol' : 'Name'}</dt>
                <dd>The ${human ? 'mouse allele (linked to MGI)' : 'human ClinVar variant'} being compared. Its rows are grouped together.</dd>
                <dt>Transcript ID</dt>
                <dd>The ${human ? 'mouse' : 'human'} transcript the comparison is made on. A variant is compared on every transcript it falls in.</dd>
                <dt>Query Transcript ID</dt>
                <dd>The transcript of the submitted variant (shown only when it falls in more than one).</dd>
                <dt>Scoring Class</dt>
                <dd><b>Coding</b> (missense, stop gained, frameshift and similar) or <b>Splice/Intronic</b>. Each class has its own criteria.</dd>
                <dt>Equivalency Tier (1–5)</dt>
                <dd>Tier 1 is the most equivalent. The line under the tier says what it rests on. Transcripts marked
                    &ldquo;not protein-coding&rdquo; (for example, nonsense-mediated decay) are shown but don't set the ${human ? 'allele' : 'variant'}'s rank.</dd>
                <dt>Match columns</dt>
                <dd><code>true</code> or <code>false</code> for each criterion, the evidence behind the tier. N/A means the criterion
                    does not apply or neither variant has a value. Columns that apply to no row are hidden.
                    <ul>
                        <li><b>Biotype Match</b>: both transcripts have the same biotype.</li>
                        <li><b>Consequence Match</b>: the variants share a main VEP consequence.</li>
                        <li><b>AA Match</b>: the same reference and variant amino acids.</li>
                        <li><b>AA Position Match</b>: the residues align in a protein alignment of the two transcripts.</li>
                        <li><b>Domain Match</b>: both residues are in Pfam domains of the same Pfam clan.</li>
                        <li><b>Intron Match</b>: the introns split the coding sequence at the same aligned residue and codon phase.</li>
                        <li><b>Intron Offset Match</b>: the same distance from the exon boundary (for example, +1).</li>
                    </ul></dd>
                <dt>Aligned Position</dt>
                <dd>The residue that aligns with the submitted variant's residue. Click it to see all of the gene's transcripts aligned.</dd>
                <dt>Local Alignment Identity</dt>
                <dd>The fraction of the 10 residues on each side that align to identical residues. Below 0.3 the aligned position is uncertain.</dd>
                ${human ? `<dt>Disease Match</dt>
                <dd><code>true</code> when the allele models at least one of the human variant's conditions (same MONDO ID). Shown beside the tier, not part of it.</dd>` : ''}
            </dl>

            <h4>Tiers</h4>
            <table class="infoTiers">
                <thead><tr><th>Tier</th><th>Coding</th><th>Splice/Intronic</th></tr></thead>
                <tbody>
                    <tr><td>1</td><td><b>Identical change</b> at the aligned residue, where the proteins align well around it (names the Pfam domain, if any)</td><td><b>Same splice site and offset</b></td></tr>
                    <tr><td>2</td><td><b>Same residue</b>, a different change; an identical change where the alignment is uncertain; or two truncating variants within 10 residues</td><td><b>Same splice site</b>, a different offset on the same side of the intron</td></tr>
                    <tr><td>3</td><td><b>Same domain</b>: same consequence, in Pfam domains of the same clan</td><td><b>Same intron</b>, deeper in the intron</td></tr>
                    <tr><td>4</td><td><b>Same consequence</b>, elsewhere in the protein</td><td><b>Different intron</b></td></tr>
                    <tr><td>5</td><td><b>Different consequence</b></td><td><b>Different consequence</b></td></tr>
                </tbody>
            </table>
            <p>${human ? 'Alleles' : 'Variants'} are ordered by their best tier on a protein-coding transcript${human ? ', then by Disease Match' : ''}.
               See the <a href="/help#4-how-variant-equivalence-is-scored" target="_blank">User Guide</a> for details.</p>`;
    },

    phenotypeSimilarity: () => `
        <p>One row per evaluated mouse allele, in the same order as the Variant Equivalency Scoring table. It shows how
           closely the allele's MGI phenotypes match the human phenotypes.</p>
        <dl>
            <dt>Mouse Allele</dt><dd>The MGI allele, linked to MGI.</dd>
            <dt>Matched Disease(s)</dt><dd>The human disease whose phenotypes best match the allele. Ties are all listed, separated by semicolons.</dd>
            <dt>MP Terms Mapped to HPO</dt><dd>How many of the allele's scored mouse (MP) terms have a translation to human (HPO) terms. A low ratio means much of the mouse phenotype can't be compared.</dd>
            <dt>Variant Phenotype Score</dt><dd>Similarity, from 0 to 1, to the phenotypes of the human variant's own disease(s).</dd>
            <dt>Gene Phenotype Score</dt><dd>Similarity, from 0 to 1, to every HPO term annotated to the human gene.</dd>
            <dt>Phenotype term mapping</dt><dd>Opens the term-by-term breakdown behind a score: how each MP term translates to HPO, its closest human term, their similarity and the arithmetic.</dd>
        </dl>
        <p>The score is how much of the allele's phenotype the human terms account for: each translated term is matched to its
           closest human term (1 for the identical term, less for a broader, narrower or related term), weighted by how exact
           the MP&rarr;HPO translation is. A dash (&mdash;) means the allele has no translatable phenotypes, not a score of zero.</p>`,

    coverage: (org) => `
        <p>Every MGI allele of ${org === 'human' ? 'the mouse ortholog' : 'the submitted gene'}, and whether VarLift could evaluate it.
           This shows how much of the gene's mouse allele collection the results cover.${org === 'mouse'
           ? ' Alleles at the submitted position are marked; alleles of other markers at the position (such as gene clusters) are included too.' : ''}</p>
        <dl>
            <dt>AlleleSymbol, AlleleID</dt><dd>The MGI allele, linked to MGI.</dd>
            <dt>Consequence</dt><dd>The allele's VEP consequence, when it has coordinates.</dd>
            <dt>Human Diseases Modeled</dt><dd>The diseases MGI curates the allele as a model of. &ldquo;No disease association&rdquo; alleles are still scored.</dd>
            <dt>Evaluated</dt><dd><code>true</code> if VarLift scored the allele. Click <code>false</code> to see why not: usually the allele has no genomic coordinates in MGI (many knockouts, gene traps and transgenes), or its consequence can't be scored (for example, a 5&prime; UTR change).</dd>
            <dt>Phenotypes</dt><dd>Click to list the allele's MGI phenotypes. &ldquo;Normal&rdquo; annotations are counted separately.</dd>
        </dl>
        <p>Evaluated alleles are listed first.</p>`,

    humanPhenotypes: () => `
        <p>Human Phenotype Ontology (HPO) terms that ClinVar lists for each candidate human variant in the Human Variants Evaluated table.
           Many ClinVar records name only a disease, so some candidates have no rows here.</p>
        <dl>
            <dt>Name</dt><dd>The ClinVar variant name.</dd>
            <dt>Phenotypes</dt><dd>The HPO terms, each with its ID.</dd>
        </dl>
        <p>These are shown for comparison with the submitted mouse variant's phenotypes; no phenotype similarity score is
           calculated for mouse queries.</p>`,

    unevaluated: (org) => `
        <p>The submitted variant has no coding, splice-site or intronic consequence (for example, it is in a UTR, a non-coding
           transcript or between genes), so it can't be compared with ${org === 'human' ? 'mouse alleles' : 'human variants'} and is not scored.</p>
        <dl>
            <dt>Gene Symbol, ${org === 'human' ? 'Mouse' : 'Human'} Homolog</dt><dd>The gene the variant is in, if any, and its ortholog.</dd>
            <dt>Consequence, Transcript, Exon / Intron</dt><dd>What VEP reports for the variant.</dd>
            <dt>Reason</dt><dd>Why the variant was not evaluated.</dd>
            ${org === 'human'
                ? '<dt>Associated Conditions</dt><dd>The variant\'s conditions in ClinVar, with their MONDO IDs.</dd>'
                : '<dt>AlleleID, AlleleSymbol, Associated Diseases</dt><dd>MGI alleles at exactly this position and change, and the diseases they model.</dd>'}
        </dl>`,

    unevaluatedPhenotypes: () => `
        <p>MGI phenotypes of the alleles at the submitted position, for alleles not already listed in the Mouse Allele Coverage table.</p>
        <dl>
            <dt>Phenotypes</dt><dd>The allele's Mammalian Phenotype (MP) annotations.</dd>
            <dt>Normal Phenotypes</dt><dd>Annotations recording that a phenotype was examined and found normal.</dd>
        </dl>`,
};


// The label entered with a submitted variant, as its first field (display only)
function withLabel(name, fields) {
    const label = (window.fullJSON.labels || {})[name];
    return label ? {'Label': escapeText(label), ...fields} : fields;
}


function variantFields(row) {
    const out = {...row};

    if (row['AlleleID'] && row['AlleleSymbol']) {
        // Several alleles at the position: each symbol is followed by its own linked MGI ID
        // (the server joins both lists in the same order)
        const ids = String(row['AlleleID']).split(', ');
        const symbols = String(row['AlleleSymbol']).split(', ');
        out['AlleleSymbol'] = ids.length === symbols.length
            ? symbols.map((symbol, i) => `${symbol} (${formatCell({AlleleID: ids[i]}, 'AlleleID')})`).join(', ')
            : `${row['AlleleSymbol']} (${formatCell(row, 'AlleleID')})`;
        delete out['AlleleID'];
    }

    // The linked Pfam ID goes in parentheses after the domain name
    if (row['Pfam Domain ID'] && row['Pfam Domain Name']) {
        out['Pfam Domain Name'] = `${row['Pfam Domain Name']} (${formatCell(row, 'Pfam Domain ID')})`;
        delete out['Pfam Domain ID'];
    }

    const match = String(row['Name'] || '').match(/^([^:]+):(\d+)-(\d+):([^/]+)\/(.+)$/);
    if (!match) return displayFields(out);

    const [, chrom, start, end, ref, alt] = match;
    const position = start === end
        ? Number(start).toLocaleString('en-US')
        : `${Number(start).toLocaleString('en-US')}-${Number(end).toLocaleString('en-US')}`;
    return displayFields({'Genome Location': `Chr ${chrom}:${position} ${ref}>${alt}`, ...out});
}


// Display names for variant columns
const DISPLAY_NAMES = {
    'Exon Rank': 'Exon',
    'Intron Rank': 'Intron',
    'Amino Acids': 'Amino Acid Change',
};

// Variant row for display: columns renamed as in DISPLAY_NAMES, and each disease
// followed by its linked MONDO ID in parentheses (the server keeps the two lists in
// step); the separate MONDO column is dropped
function displayFields(row) {
    const diseases = row['Associated Diseases'];
    const mondos = row['MONDO'];
    const paired = Array.isArray(diseases) && Array.isArray(mondos) && diseases.length === mondos.length;

    const out = {};
    Object.entries(row).forEach(([col, value]) => {
        if (col === 'MONDO' && paired) return;
        if (col === 'Associated Diseases' && paired) {
            value = diseases.map((name, i) => `${name} (${mondoLink(mondos[i])})`).join('<br>');
        }
        out[DISPLAY_NAMES[col] || col] = value;
    });
    return out;
}


// Single-record tables (the query variant, and the query gene beside its homolog) read
// down the page as Field / Value rows; several records become several value columns,
// labelled by headers when given. Fields that are empty in every record are left out.
function jsonToFieldList(data, tableId, hidden = [], headers = null) {
    if (!data || data.length === 0) return `<p>No data</p>`;

    const org = document.getElementById('queryOrganism').value;
    const isEmpty = v => v === null || v === undefined || v === '' || (Array.isArray(v) && v.length === 0);
    const fields = Object.keys(data[0])
        .filter(col => !hidden.includes(col))
        .filter(col => data.some(row => !isEmpty(row[col])));

    let html = `<table id="${tableId}" class="fieldTable${headers ? ' compare' : ''}">`;
    if (headers) {
        html += `<thead><tr><th></th>${headers.map(h => `<th scope="col">${h}</th>`).join('')}</tr></thead>`;
    }
    html += `<tbody>`;
    fields.forEach(col => {
        // A submitted human variant's diseases are ClinVar's submitted conditions
        const label = tableId === 'queryVariantTable' && org === 'human' && col === 'Associated Diseases' ? 'Associated Conditions' : col;
        html += `<tr><th scope="row">${label}</th>`;
        data.forEach(row => {
            html += `<td>${formatCell(row, col, tableId, org)}</td>`;
        });
        html += `</tr>`;
    });
    html += `</tbody></table>`;

    return html;
}


// Display value of one table cell, with links for IDs
function formatCell(row, col, tableId, org) {
    let value = row[col] !== undefined ? row[col] : '';

    // Allele table phenotypes open in a dropdown beneath the row
    if (col === 'Phenotypes' && tableId === 'alleleCoverageTable') {
        const count = value.length;
        const normal = (row['Normal Phenotypes'] || []).length;
        const label = [count ? `${count} phenotype${count === 1 ? '' : 's'}` : '', normal ? `${normal} normal` : '']
            .filter(Boolean).join(', ');
        value = !label ? 'None' : `<button type="button" class="phenotypeToggle" data-allele="${row['AlleleID']}" aria-expanded="false">${label}</button>`;
    }

    // Seperate Phenotypes into multiple lines
    else if (col === 'Phenotypes' && Array.isArray(value)) {
        value = phenotypeLinks(value);
    }
    else if (col === 'Normal Phenotypes' && Array.isArray(value)) {
        value = value.length ? `<div class="normalPhenotypes">${phenotypeLinks(value)}</div>` : '';
    }

    // HP term lists open in a dropdown
    else if (['HP_IDs_variant', 'HP_IDs_gene'].includes(col) && Array.isArray(value)) {
        const count = value.length;
        value = count === 0 ? 'None' : `<details class="termDropdown"><summary>${count} HP term${count === 1 ? '' : 's'}</summary>`
            + hpTermList(value) + '</details>';
    }

    // MONDO IDs link to their disease pages
    else if (col === 'MONDO' && Array.isArray(value)) {
        value = value.map(mondoLink).join(', ');
    }

    // Coverage alleles without a disease model
    else if (col === 'Diseases' && tableId === 'alleleCoverageTable' && Array.isArray(value) && value.length === 0) {
        value = 'No disease association';
    }

    // Other lists (HP IDs, MONDO IDs, diseases) get spaces so they can wrap
    else if (Array.isArray(value)) {
        value = value.join(', ');
    }

    // Add linking for alleles
    // A mouse variant can match several alleles, listed as "MGI:1, MGI:2"
    if (col === 'AlleleID' && value) {
        value = String(value).split(', ')
            .map(id => `<a href="https://www.informatics.jax.org/allele/${id}" target="_blank" rel="noopener noreferrer">${id}</a>`)
            .join(', ');
    }

    // An allele that wasn't evaluated says why in a dropdown under "false"
    // (the Reason column itself is hidden)
    if (col === 'Evaluated' && tableId === 'alleleCoverageTable' && value === false && row['Reason']) {
        value = `<details class="termDropdown reasonDropdown"><summary>false</summary>`
            + `<div class="reasonText">${row['Reason']}</div></details>`;
    }

    // Coverage rows stack the linked allele ID beneath the allele symbol, to save a column
    if (col === 'AlleleSymbol' && tableId === 'alleleCoverageTable' && row['AlleleID']) {
        value = `${value}<br><span class="alleleId">${formatCell(row, 'AlleleID', tableId, org)}</span>`;
    }

    // Score rows link the allele symbol to its MGI allele page (the AlleleID column is hidden)
    if (col === 'Allele Symbol' && tableId === 'scoresTable' && row['AlleleID'] && value) {
        value = `<a href="https://www.informatics.jax.org/allele/${row['AlleleID']}" target="_blank" rel="noopener noreferrer">${value}</a>`;
    }

    // Aligned Position opens a drawing of all the allele's transcripts aligned to the query transcript
    if (col === 'Aligned Position' && tableId === 'scoresTable' && value !== null && value !== '' && row['Query Transcript ID']) {
        const key = org === 'human' ? row['AlleleID'] : row['Name'];
        value = `<button type="button" class="overviewToggle" aria-expanded="false" title="Compare all transcripts"
            data-key="${escapeText(String(key ?? ''))}" data-tx="${escapeText(row['Transcript ID'] || '')}"
            data-qtx="${escapeText(row['Query Transcript ID'])}">${value}</button>`;
    }

    // Local alignment identity applies only where the residues align
    if (col === 'Local Alignment Identity' && tableId === 'scoresTable' && (value === null || value === '')) {
        value = '';
    }

    // The tier with what it rests on; NMD and non-coding transcripts don't set the allele's tier
    if (col === 'Equivalency Tier' && tableId === 'scoresTable' && value !== null && value !== '') {
        value = `<div class="tierValue">Tier ${value}</div><div class="tierBasis">${escapeText(row['Tier Basis'] || '')}</div>`
            + (row['Counts Toward Tier'] === false ? '<div class="tierNote">not protein-coding; not used for the allele\'s tier</div>' : '');
    }

    // Criteria that don't apply to the variant's scoring class
    if (tableId === 'scoresTable' && col.endsWith('Match') && value === null) {
        value = 'N/A';
    }

    // Add linking for PFAM domains
    if (col === 'Pfam Domain ID' && value) {
        value = `<a href="https://www.ebi.ac.uk/interpro/entry/pfam/${value}" target="_blank" rel="noopener noreferrer">${value}</a>`;
    }

    if (col === 'Name' && org === 'mouse' && tableId == 'outputVariantTable' && value) {
        value = `<a href="https://www.ncbi.nlm.nih.gov/clinvar?term=${value}" target="_blank" rel="noopener noreferrer">${value}</a>`;
    }

    return value;
}


// Download full JSON for all data
// File name from the time the results were produced, e.g. VarLift_results_2026-11-13_0915.json
function resultsFileName(data) {
    const t = data.run_at ? new Date(data.run_at) : new Date();
    const pad = n => String(n).padStart(2, '0');
    return `VarLift_results_${t.getFullYear()}-${pad(t.getMonth() + 1)}-${pad(t.getDate())}_${pad(t.getHours())}${pad(t.getMinutes())}.json`;
}

document.getElementById('downloadJSON').addEventListener('click', function () {
    if (!window.fullJSON) {
        alert('No results yet.');
        return;
    }

    const blob = new Blob([JSON.stringify(window.fullJSON, null, 2)], {type: 'application/json'});
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = resultsFileName(window.fullJSON);
    a.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
});


// Load Results reopens a file saved with Download All Results in the results viewer
document.getElementById('loadJSON').addEventListener('click', function () {
    document.getElementById('loadJSONFile').click();
});

document.getElementById('loadJSONFile').addEventListener('change', async function () {
    const file = this.files[0];
    this.value = '';   // choosing the same file again still triggers a load
    if (!file) return;

    showError('');
    let data;
    try {
        data = JSON.parse(await file.text());
    } catch (err) {
        showError(`${file.name} could not be read as VarLift results (it is not valid JSON).`);
        return;
    }
    if (!data || !Array.isArray(data.query_variants) || !Array.isArray(data.scores)) {
        showError(`${file.name} is not a VarLift results file. Use a file saved with Download All Results.`);
        return;
    }

    // Files saved before the organism was recorded: human queries score mouse alleles (AlleleID)
    const organism = data.organism
        || (data.scores.some(row => 'AlleleID' in row) || (data.phenotype_similarity || []).length > 0 ? 'human' : 'mouse');
    document.getElementById('queryOrganism').value = organism;
    if (typeof data.input === 'string') document.getElementById('variantInput').value = data.input;

    window.fullJSON = data;
    try {
        buildSidebar(data);
    } catch (err) {
        console.error(err);
        window.fullJSON = undefined;
        document.getElementById('resultsPanel').hidden = true;
        showError(`${file.name} could not be displayed (${err.message}).`);
    }
});
