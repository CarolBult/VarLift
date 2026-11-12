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
        // store full results
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



// Sort score rows so each allele's transcripts stay together
function groupedScoreSort(scores) {
    // Each allele's best transcript score decides where the group lands;
    // within an allele, transcripts sort by their own Total Score.
    const bestByAllele = {};
    scores.forEach(row => {
        const aid = row['AlleleID'];
        const s = row['Total Score'] ?? 0;
        if (bestByAllele[aid] === undefined || s > bestByAllele[aid]) {
            bestByAllele[aid] = s;
        }
    });

    return [...scores].sort((a, b) => {
        const aBest = bestByAllele[a['AlleleID']] ?? 0;
        const bBest = bestByAllele[b['AlleleID']] ?? 0;
        if (bBest !== aBest) return bBest - aBest;              // allele groups by best score
        if (a['AlleleID'] !== b['AlleleID']) {
            return String(a['AlleleID']).localeCompare(b['AlleleID']);  // stable tiebreak between alleles
        }
        return (b['Total Score'] ?? 0) - (a['Total Score'] ?? 0);  // within an allele
    });
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

    function addButton(label, render) {
        const btn = document.createElement('button');
        btn.className = 'sidebar-button';
        btn.textContent = label;
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

    // Mouse Model Coverage: for human queries it comes before the mouse variants, as the
    // overview of how many alleles could be evaluated ahead of the per-allele detail
    const coverage = coverageSection(alleleCoverage, org === 'human' ? homolog : gene, org);

    // Create HTML tables
    results.innerHTML = `
        <div class='tableWrapper'>
            <h2>Submitted ${queryOrganism} Variant</h2>
            ${jsonToFieldList(queryProteins.map(variantFields), "queryVariantTable", queryVariantHidden)}
        </div>

        <div class='tableWrapper'>
            <h2>Genes</h2>
            ${jsonToFieldList([...queryGenes, ...outputGenes].map(geneFields), "genesTable", [],
                              [...queryGenes.map(() => `${queryOrganism} Gene`), ...outputGenes.map(() => `${targetOrganism} Homolog`)])}
        </div>

        ${org === 'human' ? coverage : ''}

        <div class='tableWrapper'>
            <h2>${targetOrganism} Variants Evaluated</h2>
            ${jsonToHTMLTable(sortedOutputProteins.map(displayFields), "outputVariantTable")}
        </div>

        <div class='tableWrapper'>
            <h2>Variant Equivalency Scoring</h2>
            ${jsonToHTMLTable(sortedScores, "scoresTable")}
        </div>

        ${org === 'human' ? phenotypeSimilaritySection(name, gene, sortedScores) : ''}


        ${org === 'human' ? '' : coverage}

        ${org === 'human' ? '' : `
        <div class='tableWrapper'>
            <h2>${targetOrganism} Phenotypes</h2>
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
    

    // CREATE SCORE DATATABLE ELEMENT
    const scoreTable = $('#scoresTable').DataTable({
        ordering: false
    });

    // Columns that are be hidden by default
    // AlleleID is linked from the Allele Symbol instead. Phenotype similarity has its own
    // table below (one row per allele); its columns stay in the JSON download
    const scoreHiddenCols = ['Input Name', 'AlleleID', 'Matched Disease(s)', 'Phenotype Score (variant)',
                             'Phenotype Score (gene)', 'MP Terms Mapped', 'MP Terms Total'];

    // Criteria that don't apply to any row, such as the coding criteria for a splice-site variant
    const notApplicable = Object.keys(sortedScores[0] || {})
        .filter(col => col.endsWith('Match') && sortedScores.every(row => row[col] === null));

    // Hide columns
    hideColumnsByName(scoreTable, [...scoreHiddenCols, ...notApplicable]);


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


// Mouse Model Coverage table and its summary line
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
        `<button type="button" class="phenotypeToggle mappingToggle" data-allele="${row['AlleleID']}" data-kind="${kind}" aria-expanded="false">Term mapping</button>`;

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
            <h2>Phenotype Similarity</h2>
            <p class="coverageSummary">How closely each evaluated allele's MGI phenotypes match the human phenotypes.
            Mouse (MP) terms are translated to human (HPO) terms, and each translated term is matched to its closest
            human term by semantic similarity (1 for the identical term, less for a broader, narrower or related term).
            The score, from 0 to 1, is how much of the allele's phenotype the human terms account for.
            <strong>Variant</strong> compares with each disease of the human variant and reports the best;
            <strong>Gene</strong> compares with all HPO terms of ${escapeText(gene)}.
            Open <em>Term mapping</em> to see every term behind a score.</p>
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
            <h2>Mouse Model Coverage</h2>
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
    const rows = matches.map(({Phenotypes, 'Normal Phenotypes': _normal, Biotype, ...row}) => displayFields({
        ...row,
        'Associated Diseases': isEmpty(row['Associated Diseases']) ? 'None reported' : row['Associated Diseases']
    }));

    results.innerHTML = `
        <div class='tableWrapper'>
            <h2>Submitted ${queryOrganism} Variant (not evaluated)</h2>
            <p class="coverageSummary">${rows[0]['Reason']}. No homologous variants or similarity scores are reported. ${org === 'mouse'
                ? (hasAlleles ? 'The MGI allele(s) at this position are listed below.' : 'No MGI allele matches this position and allele.')
                : 'ClinVar disease associations for this variant are listed below.'}</p>
            ${jsonToHTMLTable(rows, 'unevaluatedTable')}
        </div>

        ${phenotypes.length === 0 ? '' : `
        <div class='tableWrapper'>
            <h2>${queryOrganism} Phenotypes</h2>
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

    // Coverage table: the allele's disease models are human diseases
    const label = col => (tableId === 'alleleCoverageTable' && col === 'Diseases') ? 'Human Diseases Modeled' : col;

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
function variantFields(row) {
    const out = {...row};

    if (row['AlleleID'] && row['AlleleSymbol']) {
        out['AlleleSymbol'] = `${row['AlleleSymbol']} (${formatCell(row, 'AlleleID')})`;
        delete out['AlleleID'];
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
        html += `<tr><th scope="row">${col}</th>`;
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

    // The reason is a long sentence; a wider cell keeps the rows short
    if (col === 'Reason' && tableId === 'alleleCoverageTable' && value) {
        value = `<span class="reasonText">${value}</span>`;
    }

    // Coverage rows stack the linked allele ID beneath the allele symbol, to save a column
    if (col === 'AlleleSymbol' && tableId === 'alleleCoverageTable' && row['AlleleID']) {
        value = `${value}<br><span class="alleleId">${formatCell(row, 'AlleleID', tableId, org)}</span>`;
    }

    // Score rows link the allele symbol to its MGI allele page (the AlleleID column is hidden)
    if (col === 'Allele Symbol' && tableId === 'scoresTable' && row['AlleleID'] && value) {
        value = `<a href="https://www.informatics.jax.org/allele/${row['AlleleID']}" target="_blank" rel="noopener noreferrer">${value}</a>`;
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
document.getElementById('downloadJSON').addEventListener('click', function () {
    if (!window.fullJSON) {
        alert('No results yet.');
        return;
    }

    const dataStr = 'data:text/json;charset=utf-8,' +
        encodeURIComponent(JSON.stringify(window.fullJSON, null, 2));

    const a = document.createElement('a');
    a.setAttribute('href', dataStr);
    a.setAttribute('download', 'variant_mapping_results.json');
    a.click();
});
