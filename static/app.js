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

        // build sidebar
        buildSidebar(data);
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
        phenotypes = window.fullJSON.phenotypes.filter(row => row['Gene Symbol'] === homolog);
        alleleCoverage = (window.fullJSON.allele_coverage || []).filter(row => row['Gene Symbol'] === homolog);

        // Phenotypes move into the coverage table's dropdowns; the separate
        // table keeps only alleles that are not in the coverage table
        phenotypesByAllele = new Map();
        phenotypes.forEach(row => {
            if (!phenotypesByAllele.has(row['AlleleID'])) phenotypesByAllele.set(row['AlleleID'], row);
        });
        const coverageIds = new Set(alleleCoverage.map(row => row['AlleleID']));
        phenotypes = [...phenotypesByAllele.values()].filter(row => !coverageIds.has(row['AlleleID']));
        alleleCoverage = alleleCoverage.map(row => ({
            ...row,
            Phenotypes: phenotypesByAllele.get(row['AlleleID'])?.Phenotypes || []
        }));

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
        outputProteins = window.fullJSON.output_variants.filter(row => row['Gene Symbol'] === homolog 
                                                                        && queryProteins.some(r => r.refAA === row.refAA)
                                                                        && queryProteins.some(r => r.varAA === row.varAA));
        scores = window.fullJSON.scores.filter(row => row['Input Name'] === name);
        phenotypes = window.fullJSON.phenotypes.filter(row => outputProteins.some(r => r.Name === row.Name));

        // Mouse queries list the queried gene's alleles, marking the ones this variant matches
        alleleCoverage = markSubmitted(window.fullJSON.allele_coverage, gene, queryProteins);
        otherMouseAlleles = markSubmitted(window.fullJSON.other_alleles, gene, queryProteins);
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





    // On human queries, mouse alleles outside the coverage table get their own phenotype table
    const otherAlleles = org === 'human' && alleleCoverage.length > 0;

    // Create HTML tables
    results.innerHTML = `
        <div class='tableWrapper'>
            <h2>${queryOrganism} Gene</h2>
            ${jsonToHTMLTable(queryGenes, "queryGenesTable")}
        </div>

        <div class='tableWrapper'>
            <h2>${queryOrganism} Variant</h2>
            ${jsonToHTMLTable(queryProteins, "queryVariantTable")}
        </div>

        <div class='tableWrapper'>
            <h2>${targetOrganism} Homologs</h2>
            ${jsonToHTMLTable(outputGenes, "outputGenesTable")}
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

        ${otherAlleles && phenotypes.length === 0 ? '' : `
        <div class='tableWrapper'>
            <h2>${otherAlleles ? `Other ${homolog} Alleles with Phenotypes` : `${targetOrganism} Phenotypes`}</h2>
            ${otherAlleles ? `<p class="coverageSummary">These alleles have phenotype annotations in MGI but no disease association, so they are not part of the coverage count above.</p>` : ''}
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


    // CREATE QUERY GENE DATATABLE ELEMENT
    $('#queryGenesTable').DataTable({
        dom: 't',
        pageLength: 20,
        ordering: false
        });


    // CREATE QUERY VARIANT DATATABLE ELEMENT
    const queryVariantTable = $('#queryVariantTable').DataTable({
        dom: 't',
        pageLength: 20,
        ordering: false
    });

    // Columns that are be hidden by default
    const queryHiddenCols = [
                    'Name',
                    'Submission',
                    'Gene Symbol',
                    'refAA',
                    'protein_start',
                    'varAA'
                ];

    // Extra human-only columns to hide
    const queryHumanHidden = ['MONDO_set'];

    // Hide columns
    hideColumnsByName(queryVariantTable, queryHiddenCols);

    if (queryOrganism === 'Human') {hideColumnsByName(queryVariantTable, queryHumanHidden);}


    // CREATE OUTPUT GENE DATATABLE ELEMENT
    $('#outputGenesTable').DataTable({
        dom: 't',
        pageLength: 20,
        ordering: false
    });


    // CREATE OUTPUT VARIANT DATATABLE ELEMENT
    const outputVariantTable = $('#outputVariantTable').DataTable({
        dom: 'tp',
        pageLength: 20,
        ordering: false
    });

    // Columns that are be hidden by default
    const outputHiddenCols = [
                    'Submission',
                    'Input Name',
                    'Gene Symbol',
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
        dom: 'tp',
        pageLength: 20,
        ordering: false
    });

    // Columns that are be hidden by default
    // MP Terms Mapped/Total are diagnostics for the phenotype scores; they stay in the JSON download
    const scoreHiddenCols = ['Input Name', 'MP Terms Mapped', 'MP Terms Total'];

    // Hide columns
    hideColumnsByName(scoreTable, scoreHiddenCols);


    // CREATE PHENOTYPE DATATABLE ELEMENT
    $('#phenotypeTable').DataTable({
        dom: 'tp',
        pageLength: 20,
        ordering: false
    });


    initCoverageTable(phenotypesByAllele);
    initOtherAllelesTable();

}


// Allele rows for a mouse gene, with Submitted marking the alleles the
// given variant rows match (the server marks every allele in the batch)
function markSubmitted(alleleRows, gene, variantRows) {
    const ids = new Set(variantRows.flatMap(row => String(row['AlleleID'] || '').split(', ')));
    return (alleleRows || [])
        .filter(row => row['Gene Symbol'] === gene)
        .map(row => ({...row, Submitted: ids.has(row['AlleleID'])}));
}


// Mouse query alleles with phenotypes but no disease model
function otherAllelesSection(otherAlleles, gene) {
    if (otherAlleles.length === 0) return '';
    return `
        <div class='tableWrapper'>
            <h2>Other ${gene} Alleles with Phenotypes</h2>
            <p class="coverageSummary">These alleles have phenotype annotations in MGI but no disease association, so they are not part of the coverage count above.</p>
            ${jsonToHTMLTable(otherAlleles, 'otherAlleleTable')}
        </div>`;
}


function initOtherAllelesTable() {
    const table = $('#otherAlleleTable');
    if (table.length === 0) return;

    const otherAlleleTable = table.DataTable({
        dom: 'tp',
        pageLength: 20,
        ordering: false
    });

    // Hide Gene Symbol (same for all rows of a single-gene variant query)
    const headers = table.find('thead th').map((i, th) => $(th).text().trim()).get();
    otherAlleleTable.column(headers.indexOf('Gene Symbol')).visible(false);
}


// Mouse Model Coverage table and its summary line
function coverageSection(alleleCoverage, gene, org) {
    const total = alleleCoverage.length;
    let summary;

    if (total === 0) {
        summary = 'No disease-associated mouse alleles found in MGI.';
    } else if (org === 'human') {
        const evaluated = alleleCoverage.filter(r => r.Evaluated).length;
        summary = `VarLift evaluated <strong>${evaluated} of ${total}</strong> disease-associated ${gene} alleles in MGI. ${total - evaluated} ${total - evaluated === 1 ? 'allele exists' : 'alleles exist'} but could not be cross-mapped (see table for reasons).`;
    } else {
        const submitted = alleleCoverage.filter(r => r.Submitted).map(r => r.AlleleSymbol);
        summary = `MGI has <strong>${total}</strong> disease-associated ${gene} ${total === 1 ? 'allele' : 'alleles'}. ${submitted.length > 0
            ? `The submitted variant matches ${submitted.join(', ')}.`
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

    // Human queries flag Evaluated alleles, mouse queries Submitted ones
    const headers = table.find('thead th').map((i, th) => $(th).text().trim()).get();
    const flag = headers.includes('Evaluated') ? 'Evaluated' : 'Submitted';

    const alleleCoverageTable = table.DataTable({
        dom: 'tp',
        pageLength: 20,
        ordering: true,
        order: [[headers.indexOf(flag), 'desc'], [headers.indexOf('AlleleSymbol'), 'asc']],
        columnDefs: [{ targets: headers.indexOf('Phenotypes'), orderable: false }]
    });

    // Hide Gene Symbol (same for all rows of a single-gene variant query)
    alleleCoverageTable.column(headers.indexOf('Gene Symbol')).visible(false);

    // Open or close an allele's phenotype list beneath its row
    $('#alleleCoverageTable tbody').on('click', 'button.phenotypeToggle', function () {
        const row = alleleCoverageTable.row($(this).closest('tr'));
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


// Show what is known about a variant VarLift could not evaluate
function renderUnevaluated(name) {
    const org = document.getElementById('queryOrganism').value;
    const queryOrganism = org === 'human' ? 'Human' : 'Mouse';
    const isEmpty = v => v === null || v === undefined || v === '' || (Array.isArray(v) && v.length === 0);
    const matches = window.fullJSON.unevaluated_variants.filter(row => row['Input Name'] === name);
    const hasAlleles = matches.some(row => row['AlleleID']);

    // Mouse variants inside a gene also get that gene's coverage and other-allele tables
    const gene = matches[0]['Gene Symbol'];
    const alleleCoverage = org === 'mouse' && gene ? markSubmitted(window.fullJSON.allele_coverage, gene, matches) : [];
    const otherAlleles = org === 'mouse' && gene ? markSubmitted(window.fullJSON.other_alleles, gene, matches) : [];

    // Allele phenotypes go in their own table so long lists don't stretch the rows
    // (alleles already listed in the tables above are left out)
    const listed = new Set([...alleleCoverage, ...otherAlleles].map(row => row['AlleleID']));
    const phenotypes = matches
        .filter(row => !isEmpty(row['Phenotypes']) && !listed.has(row['AlleleID']))
        .map(row => ({AlleleID: row['AlleleID'], AlleleSymbol: row['AlleleSymbol'], Phenotypes: row['Phenotypes']}));

    const rows = matches.map(({Phenotypes, ...row}) => ({
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
    initOtherAllelesTable();
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
            let value = row[col] !== undefined ? row[col] : '';

            // Coverage table phenotypes open in a dropdown beneath the row
            if (col === 'Phenotypes' && tableId === 'alleleCoverageTable') {
                const count = value.length;
                value = count === 0 ? 'None' : `<button type="button" class="phenotypeToggle" data-allele="${row['AlleleID']}" aria-expanded="false">${count} phenotype${count === 1 ? '' : 's'}</button>`;
            }

            // Seperate Phenotypes into multiple lines
            else if (col === 'Phenotypes' && Array.isArray(value)) {
                value = phenotypeLinks(value);
            }

            // Add linking for alleles
            // A mouse variant can match several alleles, listed as "MGI:1, MGI:2"
            if (col === 'AlleleID' && value) {
                value = String(value).split(', ')
                    .map(id => `<a href="https://www.informatics.jax.org/allele/${id}" target="_blank" rel="noopener noreferrer">${id}</a>`)
                    .join(', ');
            }

            // Add linking for PFAM domains
            if (col === 'Pfam Domain ID' && value) {
                value = `<a href="https://www.ebi.ac.uk/interpro/entry/pfam/${value}" target="_blank" rel="noopener noreferrer">${value}</a>`;
            }

            if (col === 'Name' && org === 'mouse' && tableId == 'outputVariantTable' && value) {
                value = `<a href="https://www.ncbi.nlm.nih.gov/clinvar?term=${value}" target="_blank" rel="noopener noreferrer">${value}</a>`;
            }

            html += `<td>${value}</td>`;
        });
        html += `</tr>`;
    });

    html += `</tbody></table>`;

    return html;
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
