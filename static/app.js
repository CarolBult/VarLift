// Run variant-model mapping for submitted variants on click.
document.getElementById('runButton').addEventListener('click', function () {
    const vars = document.getElementById('variantInput').value;
    const org = document.getElementById('queryOrganism').value;


    fetch('/run_variants', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({ variants: vars, organism: org})
    })
    .then(res => res.json())
    .then(data => {
        // store full results
        window.fullJSON = data;

        // build sidebar
        buildSidebar(data);
    })
    .catch(err => {console.error(err);});
});



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

    // Create sidebar header
    sidebar.innerHTML = '<h3>Variant Inputs:</h3>';

    // Handle no results
    if (!data.query_variants || data.query_variants.length === 0) {
        sidebar.innerHTML += '<h2>No inputs found.</h2>';
        return;
    }

    // Create sidebar buttons
    data.query_variants.forEach((row, idx) => {
        const btn = document.createElement('button');
        btn.className = 'sidebar-button';
        btn.textContent = row['Name'];
        btn.onclick = () => {
            // deactivate all other buttons
            document.querySelectorAll('#sidebar .sidebar-button')
            .forEach(b => b.classList.remove('active'));

            // activate clicked button
            btn.classList.add('active')

            // render output tables
            renderTables(row['Name'], row['Gene Symbol']); 
        }
        sidebar.appendChild(btn);
    });

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

        // Sort output proteins and score dataframes
   //     sortedScores = [...scores].sort((a, b) => (b['Total Score'] ?? 0) - (a['Total Score'] ?? 0));
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

        // Sort output proteins and score dataframes
 //       sortedScores = [...scores].sort((a, b) => (b['Total Score'] ?? 0) - (a['Total Score'] ?? 0));
        sortedScores = groupedScoreSort(scores);
        const rank = new Map(sortedScores.map((row, index) => [`${row['Name']}||${row['Transcript ID']}`, index]));

        console.log(rank)

        sortedOutputProteins = [...outputProteins].sort((a, b) =>
            (rank.get(`${a['Name']}||${a['Transcript ID']}`) ?? Infinity) -
            (rank.get(`${b['Name']}||${b['Transcript ID']}`) ?? Infinity)
        );
    }





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


        <div class='tableWrapper'>
            <h2>${targetOrganism} Phenotypes</h2>
            ${jsonToHTMLTable(phenotypes, "phenotypeTable")}
        </div>
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
    const scoreHiddenCols = ['Input Name'];

    // Hide columns
    hideColumnsByName(scoreTable, scoreHiddenCols);


    // CREATE PHENOTYPE DATATABLE ELEMENT
    $('#phenotypeTable').DataTable({
        dom: 'tp',
        pageLength: 20,
        ordering: false
    });


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

            // Seperate Phenotypes into multiple lines
            if (col === 'Phenotypes' && Array.isArray(value)) {

                if  (org === 'human') {
                    value = value
                        .map(v => {
                            const [id, label] = String(v).split(',').map(s => s.trim());
                            const url = `https://www.informatics.jax.org/vocab/mp_ontology/${id}`;
                            return `<a href="${url}" target="_blank" rel="noopener noreferrer">${id}: ${label}</a>`;
                        })
                        .filter(v => v.length)
                        .join('<br>');
                } else if (org === 'mouse') {
                        value = value
                        .map(v => {
                            const [fullID, label] = String(v).split(',').map(s => s.trim());
                            const [hp, id] = String(fullID).split(':').map(s => s.trim());
                            const url = `http://purl.obolibrary.org/obo/HP_${id}`;
                            return `<a href="${url}" target="_blank" rel="noopener noreferrer">${fullID}: ${label}</a>`;
                        })
                        .filter(v => v.length)
                        .join('<br>');
                }
            }

            // Add linking for alleles
            if (col === 'AlleleID') {
                value = `<a href="https://www.informatics.jax.org/allele/${value}" target="_blank" rel="noopener noreferrer">${value}</a>`;
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
