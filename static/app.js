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
        showError(err.message);
    });
});


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

            // render output tables
            render();
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

    // Unhide results
    document.getElementById('layoutOutput').style.display = 'flex';
}


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
        otherMouseAlleles = markEvaluated(window.fullJSON.other_alleles, gene, queryProteins);
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

    // Create HTML tables
    results.innerHTML = `
        <div class='tableWrapper'>
            <h2>${queryOrganism} Variant</h2>
            ${jsonToFieldList(queryProteins.map(variantFields), "queryVariantTable", queryVariantHidden)}
        </div>

        <div class='tableWrapper'>
            <h2>Genes</h2>
            ${jsonToFieldList([...queryGenes, ...outputGenes].map(geneFields), "genesTable", [],
                              [...queryGenes.map(() => `${queryOrganism} Gene`), ...outputGenes.map(() => `${targetOrganism} Homolog`)])}
        </div>

        <div class='tableWrapper'>
            <h2>${targetOrganism} Variants</h2>
            ${jsonToHTMLTable(sortedOutputProteins, "outputVariantTable")}
        </div>

        <div class='tableWrapper'>
            <h2>Similarity Scores</h2>
            ${jsonToHTMLTable(sortedScores, "scoresTable")}
        </div>


        ${coverageSection(alleleCoverage, org === 'human' ? homolog : gene, org)}

        ${org === 'mouse' ? otherAllelesSection(otherMouseAlleles, gene) : ''}

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
    // MP Terms Mapped/Total are diagnostics for the phenotype scores; they stay in the JSON download
    const scoreHiddenCols = ['Input Name', 'MP Terms Mapped', 'MP Terms Total'];

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
    initOtherAllelesTable(org === 'mouse' ? otherMouseAlleles : []);

}


// Allele rows for a mouse gene, with Evaluated marking the alleles the
// given variant rows match (the server marks every allele in the batch)
function markEvaluated(alleleRows, gene, variantRows) {
    const ids = new Set(variantRows.flatMap(row => String(row['AlleleID'] || '').split(', ')));
    return (alleleRows || [])
        .filter(row => row['Gene Symbol'] === gene)
        .map(row => ({...row, Evaluated: ids.has(row['AlleleID'])}));
}


// Summary line for the alleles with phenotypes but no disease association (one row per allele)
function otherAllelesSummary(count, gene) {
    return `<p class="coverageSummary">MGI has <strong>${count}</strong> other ${gene} ${count === 1 ? 'allele' : 'alleles'} with phenotype annotations but no disease association, so ${count === 1 ? 'it is' : 'they are'} not part of the coverage count above.</p>`;
}


// Mouse query alleles with phenotypes but no disease model
function otherAllelesSection(otherAlleles, gene) {
    if (otherAlleles.length === 0) return '';
    return `
        <div class='tableWrapper'>
            <h2>Other ${gene} Alleles with Phenotypes</h2>
            ${otherAllelesSummary(otherAlleles.length, gene)}
            ${jsonToHTMLTable(otherAlleles, 'otherAlleleTable')}
        </div>`;
}


function initOtherAllelesTable(otherAlleles) {
    const table = $('#otherAlleleTable');
    if (table.length === 0) return;

    const otherAlleleTable = table.DataTable({
        ordering: false
    });

    // Hide Gene Symbol (same for all rows of a single-gene variant query)
    const headers = table.find('thead th').map((i, th) => $(th).text().trim()).get();
    otherAlleleTable.column(headers.indexOf('Gene Symbol')).visible(false);

    initPhenotypeToggles('otherAlleleTable', otherAlleleTable, new Map(otherAlleles.map(row => [row['AlleleID'], row])));
}


// Open or close an allele's phenotype list beneath its row
function initPhenotypeToggles(tableId, dataTable, phenotypesByAllele) {
    $(`#${tableId} tbody`).on('click', 'button.phenotypeToggle', function () {
        const row = dataTable.row($(this).closest('tr'));
        if (row.child.isShown()) {
            row.child.hide();
            this.setAttribute('aria-expanded', 'false');
        } else {
            const phenos = phenotypesByAllele.get(this.dataset.allele).Phenotypes;
            row.child(`<div class="phenotypeList">${phenotypeLinks(phenos)}</div>`).show();
            this.setAttribute('aria-expanded', 'true');
        }
    });
}


// Mouse Model Coverage table and its summary line
function coverageSection(alleleCoverage, gene, org) {
    const total = alleleCoverage.length;
    let summary;

    if (total === 0) {
        summary = org === 'human' ? 'No mouse alleles found in MGI.' : 'No disease-associated mouse alleles found in MGI.';
    } else if (org === 'human') {
        const evaluated = alleleCoverage.filter(r => r.Evaluated).length;
        const withDisease = alleleCoverage.filter(r => r.Diseases && r.Diseases.length > 0).length;
        summary = `VarLift evaluated <strong>${evaluated} of ${total}</strong> ${gene} alleles in MGI (${withDisease} with a disease association). ${total - evaluated} ${total - evaluated === 1 ? 'allele exists' : 'alleles exist'} but could not be cross-mapped (see table for reasons).`;
    } else {
        const submitted = alleleCoverage.filter(r => r.Evaluated).map(r => r.AlleleSymbol);
        summary = `MGI has <strong>${total}</strong> disease-associated ${gene} ${total === 1 ? 'allele' : 'alleles'}. ${submitted.length > 0
            ? `The submitted variant corresponds to ${submitted.join(', ')}.`
            : 'The submitted variant does not match any of them.'}`;
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

    const alleleCoverageTable = table.DataTable({
        ordering: true,
        order: [[headers.indexOf('Evaluated'), 'desc'], [headers.indexOf('AlleleSymbol'), 'asc']],
        columnDefs: [{ targets: headers.indexOf('Phenotypes'), orderable: false }]
    });

    // Hide Gene Symbol (same for all rows of a single-gene variant query)
    alleleCoverageTable.column(headers.indexOf('Gene Symbol')).visible(false);

    initPhenotypeToggles('alleleCoverageTable', alleleCoverageTable, phenotypesByAllele);
}


// Show what is known about a variant VarLift could not evaluate
function renderUnevaluated(name) {
    const org = document.getElementById('queryOrganism').value;
    const queryOrganism = org === 'human' ? 'Human' : 'Mouse';
    const isEmpty = v => v === null || v === undefined || v === '' || (Array.isArray(v) && v.length === 0);
    const matches = window.fullJSON.unevaluated_variants.filter(row => row['Input Name'] === name);
    const hasAlleles = matches.some(row => row['AlleleID']);

    // Mouse variants inside a gene also get that gene's coverage and other-allele tables
    const gene = matches[0]['Gene Symbol'];
    const alleleCoverage = org === 'mouse' && gene ? markEvaluated(window.fullJSON.allele_coverage, gene, matches) : [];
    const otherAlleles = org === 'mouse' && gene ? markEvaluated(window.fullJSON.other_alleles, gene, matches) : [];

    // Allele phenotypes go in their own table so long lists don't stretch the rows
    // (alleles already listed in the tables above are left out)
    const listed = new Set([...alleleCoverage, ...otherAlleles].map(row => row['AlleleID']));
    const phenotypes = matches
        .filter(row => !isEmpty(row['Phenotypes']) && !listed.has(row['AlleleID']))
        .map(row => ({AlleleID: row['AlleleID'], AlleleSymbol: row['AlleleSymbol'], Phenotypes: row['Phenotypes']}));

    // Biotype is a gene property, shown only in the Genes table
    const rows = matches.map(({Phenotypes, Biotype, ...row}) => ({
        ...row,
        'Associated Diseases': isEmpty(row['Associated Diseases']) ? 'None reported' : row['Associated Diseases']
    }));

    results.innerHTML = `
        <div class='tableWrapper'>
            <h2>${queryOrganism} Variant (not evaluated)</h2>
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

        ${otherAllelesSection(otherAlleles, gene)}
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
    initOtherAllelesTable(otherAlleles);
}


// Turn "ID,label" phenotype strings into links, one per line
function phenotypeLinks(phenotypes) {
    return phenotypes
        .map(v => {
            // Split on the first comma only; labels can contain commas
            const text = String(v);
            const comma = text.indexOf(',');
            const id = text.slice(0, comma).trim();
            const label = text.slice(comma + 1).trim();
            const url = id.startsWith('MP:')
                ? `https://www.informatics.jax.org/vocab/mp_ontology/${id}`
                : `http://purl.obolibrary.org/obo/HP_${id.split(':')[1]}`;
            return `<a href="${url}" target="_blank" rel="noopener noreferrer">${id}: ${label}</a>`;
        })
        .join('<br>');
}


// Turn an array of objects into an HTML table
function jsonToHTMLTable(data, tableId) {
    if (!data || data.length === 0) return `<p>No data</p>`;

    const org = document.getElementById('queryOrganism').value;
    const columns = Object.keys(data[0]);

    let html = `<table id="${tableId}" class="table.dataTable">`;
    html += `<thead><tr>`;

    columns.forEach(col => {
        html += `<th>${col}</th>`;
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
function hpLink(id) {
    return `<a href="https://hpo.jax.org/browse/term/${id}" target="_blank" rel="noopener noreferrer">${id}</a>`;
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
        } else if (!['Start', 'End', 'Strand', 'Accession'].includes(col)) {
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

    // Human variants: each disease with its linked MONDO ID in parentheses. The lists
    // pair up when they come from ClinVar's CLNDISDB; diseases backfilled from
    // variant_summary also include OMIM/Orphanet names, so those keep a MONDO row
    const diseases = row['Associated Diseases'];
    const mondos = row['MONDO'];
    if (Array.isArray(diseases) && Array.isArray(mondos) && diseases.length > 0 && diseases.length === mondos.length) {
        out['Associated Diseases'] = diseases.map((name, i) => `${name} (${mondoLink(mondos[i])})`).join('<br>');
        delete out['MONDO'];
    }

    const match = String(row['Name'] || '').match(/^([^:]+):(\d+)-(\d+):([^/]+)\/(.+)$/);
    if (!match) return out;

    const [, chrom, start, end, ref, alt] = match;
    const position = start === end
        ? Number(start).toLocaleString('en-US')
        : `${Number(start).toLocaleString('en-US')}-${Number(end).toLocaleString('en-US')}`;
    return {'Genome Location': `Chr ${chrom}:${position} ${ref}>${alt}`, ...out};
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
    if (col === 'Phenotypes' && ['alleleCoverageTable', 'otherAlleleTable'].includes(tableId)) {
        const count = value.length;
        value = count === 0 ? 'None' : `<button type="button" class="phenotypeToggle" data-allele="${row['AlleleID']}" aria-expanded="false">${count} phenotype${count === 1 ? '' : 's'}</button>`;
    }

    // Seperate Phenotypes into multiple lines
    else if (col === 'Phenotypes' && Array.isArray(value)) {
        value = phenotypeLinks(value);
    }

    // HP term lists open in a dropdown
    else if (['HP_IDs_variant', 'HP_IDs_gene'].includes(col) && Array.isArray(value)) {
        const count = value.length;
        value = count === 0 ? 'None' : `<details class="termDropdown"><summary>${count} HP term${count === 1 ? '' : 's'}</summary>`
            + `<div class="phenotypeList">${value.map(hpLink).join('<br>')}</div></details>`;
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
