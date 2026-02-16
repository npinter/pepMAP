const body = document.body;
const DEFAULT_SESSION_SECONDS = 1800;
const parsedSessionTime = Number.parseInt(body.dataset.sessionSeconds || '', 10);
const sessionTime = Number.isFinite(parsedSessionTime) ? parsedSessionTime : DEFAULT_SESSION_SECONDS;
const staticBase = body.dataset.staticBase || '/static/';
const REPORT_LOGO_MAP = {
    'dia-nn': 'DIA-NN.png',
    'spectronaut': 'Spectronaut.png',
    'maxquant': 'MaxQuant.png',
    'fragpipe': 'FragPipe.png'
};
const SAMPLE_CLEANUP_DEBOUNCE_MS = 150;
const AUTOCOMPLETE_DEBOUNCE_MS = 200;
let sessionId = null;
let timerInterval = null;
let hasCustomFeatures = false;
let customFeaturesLabel = '';
let geneEntries = [];
let rawRuns = [];
let cleanedRuns = [];
let selectedRuns = new Set();
let selectionClearedManually = false;
let advancedFieldsVisible = false;
let sampleCleanupTimer = null;
let sampleCleanupAbortController = null;
let sampleCleanupRequestId = 0;
let autocompleteAbortController = null;
let autocompleteRequestId = 0;
let autocompleteDebounceTimer = null;

const elements = {
    reportFile: document.getElementById('report_file'),
    fastaFile: document.getElementById('fasta_file'),
    customFeaturesFile: document.getElementById('custom_features_file'),
    reportFilename: document.getElementById('report_filename'),
    fastaFilename: document.getElementById('fasta_filename'),
    customFeaturesFilename: document.getElementById('custom_features_filename'),
    customFeaturesPlaceholder: document.getElementById('custom_features_placeholder'),
    organism: document.getElementById('organism'),
    customFeaturesLabel: document.getElementById('custom_features_label'),
    customFeaturesLabelContainer: document.getElementById('custom_features_label_container'),
    uploadButton: document.getElementById('upload_button'),
    newReportButton: document.getElementById('new_report_button'),
    reportTypeBadge: document.getElementById('report_type_badge'),
    reportTypeLogo: document.getElementById('report_type_logo'),
    customOrganismContainer: document.getElementById('custom_organism_container'),
    customOrganismInput: document.getElementById('custom_organism'),
    mappingBlock: document.getElementById('mapping_block'),
    searchContainer: document.getElementById('search_container'),
    geneEntryInput: document.getElementById('gene_entry_input'),
    geneSuggestions: document.getElementById('gene_suggestions'),
    geneSuggestionsStatus: document.getElementById('gene_suggestions_status'),
    addGeneButton: document.getElementById('add_gene_button'),
    geneList: document.getElementById('gene_list'),
    searchInput: document.getElementById('search_input'),
    searchLabels: document.getElementById('search_labels'),
    mapContainer: document.getElementById('map_container'),
    submitButton: document.getElementById('submit_button'),
    proteotypicCheckbox: document.getElementById('proteotypic_checkbox'),
    chargeStateMode: document.getElementById('charge_state_mode'),
    qValueCutoff: document.getElementById('q_value_cutoff'),
    proteotypicContainer: document.getElementById('proteotypic_container'),
    chargeStateContainer: document.getElementById('charge_state_container'),
    qValueContainer: document.getElementById('q_value_container'),
    sampleNameContainer: document.getElementById('sample_name_container'),
    statsContainer: document.getElementById('stats_container'),
    summaryMode: document.getElementById('summary_mode'),
    showAllStatsButton: document.getElementById('show_all_stats_button'),
    clearStatsButton: document.getElementById('clear_stats_button'),
    timer: document.getElementById('timer'),
    peptidesPlot: document.getElementById('peptides_plot'),
    peptidesPlotContent: document.getElementById('peptides_plot_content'),
    featuresPlot: document.getElementById('features_plot'),
    peptidesTable: document.getElementById('peptides_table'),
    copyTableButton: document.getElementById('copy_table_button'),
    submitHint: document.getElementById('submit_hint'),
    sampleNameCleanup: document.getElementById('sample_name_cleanup'),
    sampleNameCustomContainer: document.getElementById('sample_name_custom_container'),
    sampleNameCustomPattern: document.getElementById('sample_name_custom_pattern'),
    sampleNamePreview: document.getElementById('sample_name_preview'),
    samplePickerContainer: document.getElementById('sample_picker_container'),
    samplePickerList: document.getElementById('sample_picker_list'),
    samplePickerCount: document.getElementById('sample_picker_count'),
    sampleSelectAll: document.getElementById('sample_select_all'),
    sampleClearAll: document.getElementById('sample_clear_all'),
    sessionShareContainer: document.getElementById('session_share_container'),
    sessionShareLink: document.getElementById('session_share_link'),
    sessionShareCopy: document.getElementById('session_share_copy'),
    sessionShareHint: document.getElementById('session_share_hint')
};

function setHidden(element, hidden) {
    if (!element) {
        return;
    }
    element.classList.toggle('hidden', hidden);
}

function setFileInputLocked(input, locked) {
    if (!input) {
        return;
    }
    input.disabled = locked;
    if (locked) {
        input.setAttribute('aria-disabled', 'true');
    } else {
        input.removeAttribute('aria-disabled');
    }
}

function buildSessionShareUrl(id) {
    if (!id) {
        return '';
    }
    const url = new URL(window.location.href);
    url.searchParams.set('session_id', id);
    return url.toString();
}

function clearSessionUrlParams() {
    const url = new URL(window.location.href);
    ['session_id', 'genes', 'labels'].forEach(param => url.searchParams.delete(param));
    const query = url.searchParams.toString();
    const nextUrl = url.pathname + (query ? `?${query}` : '') + url.hash;
    window.history.replaceState({}, '', nextUrl);
}

function updateSessionShare(id) {
    if (!elements.sessionShareContainer || !elements.sessionShareLink) {
        return;
    }
    if (!id) {
        elements.sessionShareLink.textContent = '';
        elements.sessionShareLink.dataset.url = '';
        setHidden(elements.sessionShareContainer, true);
        setHidden(elements.sessionShareHint, true);
        return;
    }
    const url = buildSessionShareUrl(id);
    elements.sessionShareLink.textContent = url;
    elements.sessionShareLink.dataset.url = url;
    setHidden(elements.sessionShareContainer, false);
    setHidden(elements.sessionShareHint, false);
}

function renderLoader() {
    return (
        '<div class="flex justify-center py-6">'
        + '<div class="h-6 w-6 animate-spin rounded-full border-2 border-slate-300 border-t-slate-900"></div>'
        + '</div>'
    );
}

function setButtonLoading(button, loading, label) {
    if (!button) {
        return;
    }
    if (loading) {
        button.dataset.originalMinWidth = button.style.minWidth || '';
        button.dataset.originalMinHeight = button.style.minHeight || '';
        button.style.minWidth = `${button.getBoundingClientRect().width}px`;
        button.style.minHeight = `${button.getBoundingClientRect().height}px`;
        button.dataset.original = button.innerHTML;
        const spinner = '<span class="inline-flex h-4 w-4 animate-spin rounded-full border-2 border-slate-400 border-t-slate-700"></span>';
        if (label) {
            button.innerHTML = `<span class="inline-flex items-center gap-2">${spinner}<span>${label}</span></span>`;
        } else {
            button.innerHTML = spinner;
        }
        button.disabled = true;
    } else {
        button.innerHTML = button.dataset.original || button.innerHTML;
        button.style.minWidth = button.dataset.originalMinWidth || '';
        button.style.minHeight = button.dataset.originalMinHeight || '';
        button.disabled = false;
    }
}

function setAdvancedFieldsVisible(visible) {
    advancedFieldsVisible = visible;
    setHidden(elements.proteotypicContainer, !visible);
    setHidden(elements.chargeStateContainer, !visible);
    setHidden(elements.qValueContainer, !visible);
    setHidden(elements.sampleNameContainer, !visible);
    setHidden(elements.samplePickerContainer, !visible);
    if (!visible) {
        setHidden(elements.sampleNameCustomContainer, true);
    }
}

function setStatsControlsVisible(visible) {
    if (elements.showAllStatsButton) {
        setHidden(elements.showAllStatsButton, !visible);
    }
    if (elements.clearStatsButton) {
        setHidden(elements.clearStatsButton, !visible);
    }
}

function setStatsVisibility(visible) {
    setHidden(elements.statsContainer, !visible);
    setStatsControlsVisible(visible);
}

function updateCustomFeaturesUI({ hasCustomFeatures, filename = '', label = '', reportReady = false }) {
    if (hasCustomFeatures) {
        if (filename && elements.customFeaturesFilename) {
            elements.customFeaturesFilename.value = filename;
            setHidden(elements.customFeaturesFilename, false);
        } else if (elements.customFeaturesFilename) {
            setHidden(elements.customFeaturesFilename, true);
        }
        setHidden(elements.customFeaturesFile, true);
        setFileInputLocked(elements.customFeaturesFile, true);
        if (elements.customFeaturesPlaceholder) {
            setHidden(elements.customFeaturesPlaceholder, true);
        }
        elements.customFeaturesLabel.value = label;
        elements.customFeaturesLabel.disabled = true;
        setHidden(elements.customFeaturesLabelContainer, false);
        return;
    }

    setHidden(elements.customFeaturesLabelContainer, true);
    if (elements.customFeaturesFilename) {
        elements.customFeaturesFilename.value = '';
        setHidden(elements.customFeaturesFilename, true);
    }
    elements.customFeaturesFile.value = '';
    if (reportReady) {
        setHidden(elements.customFeaturesFile, true);
        if (elements.customFeaturesPlaceholder) {
            setHidden(elements.customFeaturesPlaceholder, false);
        }
        setFileInputLocked(elements.customFeaturesFile, true);
    } else {
        setHidden(elements.customFeaturesFile, false);
        if (elements.customFeaturesPlaceholder) {
            setHidden(elements.customFeaturesPlaceholder, true);
        }
        setFileInputLocked(elements.customFeaturesFile, false);
    }
}

function getSampleCleanupRegex() {
    if (!elements.sampleNameCleanup) {
        return { regex: null, error: '' };
    }
    const mode = elements.sampleNameCleanup.value;
    let pattern = '';
    if (mode === 'split_underscore') {
        pattern = '^([^_]+)';
    } else if (mode === 'custom') {
        pattern = (elements.sampleNameCustomPattern?.value || '').trim();
    }
    if (!pattern) {
        return { regex: null, error: '' };
    }
    try {
        return { regex: new RegExp(pattern), error: '' };
    } catch (error) {
        return { regex: null, error: 'Invalid regex' };
    }
}

function applySampleCleanup(value, regexInfo) {
    if (!regexInfo || !regexInfo.regex) {
        return value;
    }
    const regex = regexInfo.regex;
    regex.lastIndex = 0;
    const match = regex.exec(value);
    if (!match) {
        return value;
    }
    if (match.length > 1 && match[1] !== undefined) {
        return match[1] ?? value;
    }
    return match[0] ?? value;
}

function deriveCleanedRuns(regexInfo) {
    const cleaned = [];
    const seen = new Set();
    rawRuns.forEach(run => {
        const text = String(run);
        const cleanedRun = regexInfo && regexInfo.regex ? applySampleCleanup(text, regexInfo) : text;
        if (!seen.has(cleanedRun)) {
            seen.add(cleanedRun);
            cleaned.push(cleanedRun);
        }
    });
    return cleaned;
}

function updateSampleNamePreview() {
    if (!elements.sampleNamePreview) {
        return;
    }
    if (!rawRuns.length) {
        elements.sampleNamePreview.textContent = '';
        return;
    }
    const regexInfo = getSampleCleanupRegex();
    if (regexInfo.error) {
        elements.sampleNamePreview.textContent = `Preview: ${regexInfo.error}`;
        return;
    }
    const original = String(rawRuns[0]);
    const cleaned = applySampleCleanup(original, regexInfo);
    elements.sampleNamePreview.textContent = `Preview: ${original} -> ${cleaned}`;
}

function updateSamplePickerCount() {
    if (!elements.samplePickerCount) {
        return;
    }
    if (!cleanedRuns.length) {
        elements.samplePickerCount.textContent = '';
        return;
    }
    elements.samplePickerCount.textContent = `${selectedRuns.size} of ${cleanedRuns.length} selected`;
}

function renderSamplePicker(forceReset = false, cleanedOverride = null) {
    if (!elements.samplePickerList || !elements.samplePickerContainer) {
        return;
    }
    let cleanedNext = [];
    if (Array.isArray(cleanedOverride)) {
        cleanedNext = cleanedOverride.map(run => String(run));
    } else {
        const regexInfo = getSampleCleanupRegex();
        cleanedNext = regexInfo.error ? deriveCleanedRuns(null) : deriveCleanedRuns(regexInfo);
    }
    const previousRuns = cleanedRuns;
    const previousSelection = new Set(selectedRuns);
    const hadFullSelection = previousRuns.length > 0 && previousSelection.size === previousRuns.length;

    cleanedRuns = cleanedNext;

    if (forceReset) {
        selectedRuns = selectionClearedManually ? new Set() : new Set(cleanedRuns);
    } else if (!previousSelection.size) {
        if (selectionClearedManually) {
            selectedRuns = new Set();
        } else {
            selectedRuns = new Set(cleanedRuns);
        }
    } else if (hadFullSelection) {
        selectedRuns = new Set(cleanedRuns);
    } else {
        selectedRuns = new Set(cleanedRuns.filter(run => previousSelection.has(run)));
        if (!selectedRuns.size && cleanedRuns.length && !selectionClearedManually) {
            selectedRuns = new Set(cleanedRuns);
        }
    }

    elements.samplePickerList.innerHTML = '';
    if (!cleanedRuns.length) {
        setHidden(elements.samplePickerContainer, true);
        updateSamplePickerCount();
        return;
    }
    if (advancedFieldsVisible) {
        setHidden(elements.samplePickerContainer, false);
    }

    const fragment = document.createDocumentFragment();
    cleanedRuns.forEach(run => {
        const label = document.createElement('label');
        label.className = 'flex items-center gap-2 rounded-md px-2 py-1 hover:bg-white';

        const input = document.createElement('input');
        input.type = 'checkbox';
        input.value = run;
        input.checked = selectedRuns.has(run);
        input.className = 'h-3.5 w-3.5 rounded border-slate-300 text-slate-700 focus-visible:ring-2 focus-visible:ring-slate-400';
        input.addEventListener('change', function() {
            if (input.checked) {
                selectedRuns.add(run);
                selectionClearedManually = false;
            } else {
                selectedRuns.delete(run);
                if (!selectedRuns.size) {
                    selectionClearedManually = true;
                }
            }
            updateSamplePickerCount();
        });

        const span = document.createElement('span');
        span.className = 'truncate';
        span.textContent = run;

        label.append(input, span);
        fragment.append(label);
    });
    elements.samplePickerList.append(fragment);
    updateSamplePickerCount();
}

function requestServerSampleCleanup() {
    if (!sessionId || !rawRuns.length) {
        return;
    }
    if (sampleCleanupAbortController) {
        sampleCleanupAbortController.abort();
    }
    const requestId = ++sampleCleanupRequestId;
    const controller = new AbortController();
    sampleCleanupAbortController = controller;

    const payload = {
        session_id: sessionId,
        cleanup_mode: elements.sampleNameCleanup?.value || 'none',
        custom_pattern: elements.sampleNameCustomPattern?.value || ''
    };

    fetch('/sample_cleanup_preview', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload),
        signal: controller.signal
    })
        .then(response => response.json())
        .then(data => {
            if (requestId !== sampleCleanupRequestId) {
                return;
            }
            if (data.error) {
                if (elements.sampleNamePreview) {
                    elements.sampleNamePreview.textContent = `Preview: ${data.error}`;
                }
            } else if (data.preview && elements.sampleNamePreview) {
                elements.sampleNamePreview.textContent = `Preview: ${data.preview.original} -> ${data.preview.cleaned}`;
            }
            if (Array.isArray(data.cleaned_runs)) {
                renderSamplePicker(false, data.cleaned_runs);
            }
        })
        .catch(error => {
            if (controller.signal.aborted) {
                return;
            }
        });
}

function scheduleSampleCleanupUpdate() {
    if (sampleCleanupTimer) {
        clearTimeout(sampleCleanupTimer);
    }
    sampleCleanupTimer = setTimeout(function() {
        updateSampleNamePreview();
        renderSamplePicker(false);
        requestServerSampleCleanup();
    }, SAMPLE_CLEANUP_DEBOUNCE_MS);
}

function setAvailableRuns(runs, forceReset = true) {
    rawRuns = Array.isArray(runs) ? runs.map(run => String(run)) : [];
    cleanedRuns = [];
    if (forceReset) {
        selectedRuns = new Set();
        selectionClearedManually = false;
    }
    updateSampleNamePreview();
    renderSamplePicker(true);
    requestServerSampleCleanup();
}

function getSelectedRunsPayload() {
    if (!cleanedRuns.length) {
        return null;
    }
    if (!selectedRuns.size) {
        return [];
    }
    if (selectedRuns.size === cleanedRuns.length) {
        return null;
    }
    return Array.from(selectedRuns);
}

function updateReportTypeBadge(reportType, reportLogo) {
    if (!elements.reportTypeBadge) {
        return;
    }
    if (!elements.reportTypeLogo || !elements.reportTypeLogo.isConnected) {
        const logo = document.createElement('img');
        logo.id = 'report_type_logo';
        logo.className = 'h-12 w-12 object-contain';
        logo.alt = '';
        elements.reportTypeBadge.appendChild(logo);
        elements.reportTypeLogo = logo;
    }
    if (!reportType) {
        setHidden(elements.reportTypeBadge, true);
        elements.reportTypeLogo.removeAttribute('src');
        elements.reportTypeLogo.setAttribute('alt', '');
        return;
    }
    const normalized = reportType.toLowerCase();
    let logoFile = reportLogo || REPORT_LOGO_MAP[normalized];
    if (logoFile) {
        const parts = logoFile.split(/[/\\\\]/);
        logoFile = parts[parts.length - 1];
    }
    if (!logoFile) {
        setHidden(elements.reportTypeBadge, true);
        return;
    }
    setHidden(elements.reportTypeBadge, false);
    elements.reportTypeLogo.src = `${staticBase}res/${logoFile}`;
    elements.reportTypeLogo.setAttribute('alt', reportType);
}

function setHtmlAndRunScripts(container, html) {
    const temp = document.createElement('div');
    temp.innerHTML = html;
    container.innerHTML = '';
    while (temp.firstChild) {
        container.appendChild(temp.firstChild);
    }
    const scripts = container.querySelectorAll('script');
    scripts.forEach(script => {
        const newScript = document.createElement('script');
        Array.from(script.attributes).forEach(attr => {
            newScript.setAttribute(attr.name, attr.value);
        });
        if (script.text) {
            newScript.text = script.text;
        }
        document.body.appendChild(newScript);
        script.remove();
    });
}

function applySessionData(data) {
    sessionId = data.session_id;
    updateSessionShare(sessionId);
    if (data.report_filename) {
        elements.reportFilename.value = data.report_filename;
        setHidden(elements.reportFilename, false);
        setHidden(elements.reportFile, true);
    }
    if (data.fasta_filename) {
        elements.fastaFilename.value = data.fasta_filename;
        setHidden(elements.fastaFilename, false);
        setHidden(elements.fastaFile, true);
    }
    if (data.organism) {
        const optionValues = Array.from(elements.organism.options).map(option => option.value);
        if (optionValues.includes(data.organism)) {
            elements.organism.value = data.organism;
            elements.organism.disabled = true;
        }
    }

    hasCustomFeatures = Boolean(data.has_custom_features);
    customFeaturesLabel = data.custom_features_label || '';
    updateCustomFeaturesUI({
        hasCustomFeatures,
        filename: data.custom_features_filename || '',
        label: customFeaturesLabel,
        reportReady: Boolean(data.report_filename && data.fasta_filename)
    });

    const settings = data.settings || {};
    if (elements.proteotypicCheckbox) {
        elements.proteotypicCheckbox.checked = Boolean(settings.proteotypic_only);
    }
    if (elements.chargeStateMode && settings.charge_state_mode) {
        elements.chargeStateMode.value = settings.charge_state_mode;
    }
    if (elements.qValueCutoff && settings.q_value_cutoff !== undefined) {
        elements.qValueCutoff.value = String(settings.q_value_cutoff);
    }
    if (elements.sampleNameCleanup && settings.sample_name_cleanup) {
        elements.sampleNameCleanup.value = settings.sample_name_cleanup;
    }
    if (elements.sampleNameCustomPattern && settings.sample_name_custom_pattern) {
        elements.sampleNameCustomPattern.value = settings.sample_name_custom_pattern;
    }
    if (elements.summaryMode && settings.summary_mode) {
        elements.summaryMode.value = settings.summary_mode;
    }
    setHidden(elements.sampleNameCustomContainer, elements.sampleNameCleanup.value !== 'custom');
    setHidden(elements.uploadButton, true);
    setHidden(elements.searchContainer, false);
    setHidden(elements.mapContainer, false);
    setHidden(elements.newReportButton, false);
    setHidden(elements.mappingBlock, false);
    setHidden(elements.timer, false);
    setAdvancedFieldsVisible(true);
    setAvailableRuns(data.runs || [], true);
    updateReportTypeBadge(data.report_type || '', data.report_logo || '');
    resetTimer();
    updateMapButtonState();
}

function startTimer(duration, display) {
    let timer = duration;
    const interval = setInterval(function() {
        const minutes = String(Math.floor(timer / 60)).padStart(2, '0');
        const seconds = String(timer % 60).padStart(2, '0');
        display.textContent = `Session will expire in ${minutes}:${seconds} minutes`;

        if (--timer < 0) {
            clearInterval(interval);
            if (!sessionId) {
                location.reload();
                return;
            }
            const formData = new FormData();
            formData.append('session_id', sessionId);
            fetch('/flush', { method: 'POST', body: formData })
                .then(response => response.json())
                .then(() => location.reload())
                .catch(() => location.reload());
        }
    }, 1000);
    return interval;
}

function resetTimer() {
    if (!elements.timer) {
        return;
    }
    if (timerInterval) {
        clearInterval(timerInterval);
    }
    timerInterval = startTimer(sessionTime, elements.timer);
}

function stopTimer() {
    if (timerInterval) {
        clearInterval(timerInterval);
        timerInterval = null;
    }
}

function parseGeneEntries(value) {
    if (!value) {
        return [];
    }
    return value.split(/[\s,;]+/).map(token => token.trim()).filter(Boolean);
}

function clearAutocompleteSuggestions() {
    if (elements.geneSuggestions) {
        elements.geneSuggestions.innerHTML = '';
    }
    if (elements.geneSuggestionsStatus) {
        elements.geneSuggestionsStatus.textContent = '';
    }
}

function cancelAutocompleteRequests() {
    if (autocompleteDebounceTimer) {
        clearTimeout(autocompleteDebounceTimer);
        autocompleteDebounceTimer = null;
    }
    if (autocompleteAbortController) {
        autocompleteAbortController.abort();
        autocompleteAbortController = null;
    }
    autocompleteRequestId += 1;
}

function scheduleAutocomplete(query) {
    cancelAutocompleteRequests();
    if (!query || !sessionId) {
        clearAutocompleteSuggestions();
        return;
    }
    const requestId = ++autocompleteRequestId;
    autocompleteDebounceTimer = setTimeout(function() {
        const controller = new AbortController();
        autocompleteAbortController = controller;
        const sessionQuery = sessionId ? `&session_id=${encodeURIComponent(sessionId)}` : '';
        fetch(`/autocomplete?query=${encodeURIComponent(query)}${sessionQuery}`, { signal: controller.signal })
            .then(response => response.json())
            .then(data => {
                if (requestId !== autocompleteRequestId) {
                    return;
                }
                const suggestions = data.suggestions || [];
                const options = suggestions.map(suggestion =>
                    `<option value="${suggestion}"></option>`
                );
                if (elements.geneSuggestions) {
                    elements.geneSuggestions.innerHTML = options.join('');
                }
                if (elements.geneSuggestionsStatus) {
                    elements.geneSuggestionsStatus.textContent = suggestions.length
                        ? `${suggestions.length} match${suggestions.length === 1 ? '' : 'es'}`
                        : 'No matches found';
                }
            })
            .catch(error => {
                if (controller.signal.aborted) {
                    return;
                }
                if (requestId !== autocompleteRequestId) {
                    return;
                }
                if (elements.geneSuggestionsStatus) {
                    elements.geneSuggestionsStatus.textContent = 'No matches found';
                }
            });
    }, AUTOCOMPLETE_DEBOUNCE_MS);
}

function getFilenameStem(filename) {
    if (!filename) {
        return '';
    }
    const lastDot = filename.lastIndexOf('.');
    if (lastDot <= 0) {
        return filename;
    }
    return filename.slice(0, lastDot);
}

function updateHiddenSearchInput() {
    elements.searchInput.value = geneEntries.map(entry => entry.query).join(' ');
    elements.searchLabels.value = JSON.stringify(geneEntries.map(entry => entry.label));
}

function renderGeneList() {
    elements.geneList.innerHTML = '';
    geneEntries.forEach((entry, index) => {
        const row = document.createElement('div');
        row.className = 'group flex items-center gap-2 rounded-full border border-slate-200 bg-slate-100 px-3 py-1 text-xs text-slate-700';

        const input = document.createElement('input');
        input.type = 'text';
        input.value = entry.label;
        input.disabled = true;
        input.title = entry.query;
        input.className = 'w-28 bg-transparent text-xs text-slate-700 focus:outline-none';

        const actions = document.createElement('div');
        actions.className = 'flex items-center gap-2 opacity-0 transition group-hover:opacity-100 group-focus-within:opacity-100';

        const edit = document.createElement('button');
        edit.type = 'button';
        edit.title = 'Edit';
        edit.className = 'text-slate-500 hover:text-slate-900';
        edit.innerHTML = '<i class="fa-solid fa-pen"></i>';

        const remove = document.createElement('button');
        remove.type = 'button';
        remove.title = 'Remove';
        remove.className = 'text-slate-500 hover:text-slate-900';
        remove.innerHTML = '<i class="fa-solid fa-trash"></i>';

        edit.addEventListener('click', function() {
            const isDisabled = input.disabled;
            input.disabled = !isDisabled;
            if (!isDisabled) {
                const value = input.value.trim();
                if (!value) {
                    geneEntries.splice(index, 1);
                    renderGeneList();
                    return;
                }
                geneEntries[index].label = value;
                geneEntries[index].locked = true;
                renderGeneList();
            } else {
                input.focus();
            }
        });

        input.addEventListener('keydown', function(event) {
            if (event.key === 'Enter') {
                event.preventDefault();
                if (!input.disabled) {
                    edit.click();
                }
            }
        });

        remove.addEventListener('click', function() {
            geneEntries.splice(index, 1);
            renderGeneList();
        });

        actions.append(edit, remove);
        row.append(input, actions);
        elements.geneList.append(row);
    });
    updateHiddenSearchInput();
    updateMapButtonState();
}

function dedupeGeneEntries() {
    const seen = new Set();
    geneEntries = geneEntries.filter(entry => {
        const key = entry.query.toUpperCase();
        if (seen.has(key)) {
            return false;
        }
        seen.add(key);
        return true;
    });
    renderGeneList();
}

async function resolveGeneEntries(startIndex) {
    if (!sessionId || startIndex >= geneEntries.length) {
        return;
    }
    const identifiers = geneEntries.slice(startIndex).map(entry => entry.query);
    if (!identifiers.length) {
        return;
    }
    const response = await fetch('/resolve_labels', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ session_id: sessionId, identifiers })
    });
    if (!response.ok) {
        return;
    }
    const data = await response.json();
    const resolved = data.resolved || [];
    resolved.forEach((item, offset) => {
        const index = startIndex + offset;
        if (!geneEntries[index]) {
            return;
        }
        if (item.uniprot_id) {
            geneEntries[index].query = item.uniprot_id;
        }
        if (!geneEntries[index].locked && item.label) {
            geneEntries[index].label = item.label;
        }
    });
    dedupeGeneEntries();
}

function addGeneEntries(entries) {
    if (!entries.length) {
        return Promise.resolve();
    }
    const existing = new Set(geneEntries.map(entry => entry.query.toUpperCase()));
    const startIndex = geneEntries.length;
    entries.forEach(entry => {
        const key = entry.toUpperCase();
        if (!existing.has(key)) {
            geneEntries.push({ query: entry, label: entry, locked: false });
            existing.add(key);
        }
    });
    renderGeneList();
    return resolveGeneEntries(startIndex);
}

function addGeneEntriesWithLabels(entries, labels) {
    if (!entries.length) {
        return Promise.resolve();
    }
    const existing = new Set(geneEntries.map(entry => entry.query.toUpperCase()));
    const startIndex = geneEntries.length;
    entries.forEach((entry, index) => {
        const key = entry.toUpperCase();
        if (existing.has(key)) {
            return;
        }
        const label = (labels && labels[index]) ? String(labels[index]) : entry;
        geneEntries.push({ query: entry, label, locked: Boolean(labels && labels[index]) });
        existing.add(key);
    });
    renderGeneList();
    return resolveGeneEntries(startIndex);
}

function addGeneFromInput() {
    const value = elements.geneEntryInput.value;
    const entries = parseGeneEntries(value);
    cancelAutocompleteRequests();
    clearAutocompleteSuggestions();
    const promise = addGeneEntries(entries);
    elements.geneEntryInput.value = '';
    return promise;
}

function updateMapButtonState() {
    if (!elements.submitButton) {
        return;
    }
    const hasGenes = geneEntries.length > 0;
    elements.submitButton.disabled = !hasGenes;
    if (elements.submitHint) {
        setHidden(elements.submitHint, hasGenes);
    }
}

function resetFormState() {
    const form = document.querySelector('form');
    if (form) {
        form.reset();
    }
    cancelAutocompleteRequests();
    clearAutocompleteSuggestions();
    setHidden(elements.sampleNameCustomContainer, true);
    setAvailableRuns([], true);
    setStatsVisibility(false);

    setHidden(elements.customFeaturesFilename, true);
    elements.customFeaturesFilename.value = '';
    setHidden(elements.customFeaturesFile, false);
    if (elements.customFeaturesPlaceholder) {
        setHidden(elements.customFeaturesPlaceholder, true);
    }
    elements.customFeaturesFile.value = '';
    setFileInputLocked(elements.customFeaturesFile, false);

    elements.customFeaturesLabel.disabled = false;
    elements.customFeaturesLabel.value = '';
    setHidden(elements.customFeaturesLabelContainer, true);

    setHidden(elements.mapContainer, true);
    elements.reportFilename.value = '';
    elements.fastaFilename.value = '';
    setHidden(elements.reportFilename, true);
    setHidden(elements.fastaFilename, true);
    setHidden(elements.reportFile, false);
    setHidden(elements.fastaFile, false);
    hasCustomFeatures = false;
    customFeaturesLabel = '';
    sessionId = null;
    geneEntries = [];
    elements.geneEntryInput.value = '';
    setHidden(elements.searchContainer, true);
    setHidden(elements.newReportButton, true);
    setHidden(elements.uploadButton, false);
    setHidden(elements.mappingBlock, true);
    setHidden(elements.peptidesPlot, true);
    setHidden(elements.featuresPlot, true);
    elements.peptidesPlotContent.innerHTML = '';
    elements.featuresPlot.innerHTML = '';
    elements.peptidesTable.innerHTML = '';
    setHidden(elements.peptidesTable, true);
    if (elements.copyTableButton) {
        setHidden(elements.copyTableButton, true);
    }
    setHidden(elements.timer, true);
    setHidden(elements.reportTypeBadge, true);
    updateReportTypeBadge('', '');
    stopTimer();
    updateSessionShare(null);
    if (elements.chargeStateMode) {
        elements.chargeStateMode.value = 'all';
    }
    if (elements.qValueCutoff) {
        elements.qValueCutoff.value = '0.01';
    }
    if (elements.summaryMode) {
        elements.summaryMode.value = 'per_sample';
    }
    elements.organism.disabled = false;
    if (elements.customOrganismInput) {
        elements.customOrganismInput.value = '';
    }
    setHidden(elements.customOrganismContainer, true);
    setAdvancedFieldsVisible(false);
    renderGeneList();
}
