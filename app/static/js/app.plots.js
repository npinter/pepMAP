function getPersistentHoverLayer(plotDiv) {
    const mainSvgs = plotDiv.querySelectorAll('.main-svg');
    if (mainSvgs.length === 0) {
        return null;
    }

    const targetSvg = mainSvgs[1] || mainSvgs[0];
    const layerAbove = targetSvg.querySelector('.layer-above');
    if (!layerAbove) {
        return null;
    }

    let persistentHoverLayer = layerAbove.querySelector('.persistent-hoverlayer');
    if (!persistentHoverLayer) {
        persistentHoverLayer = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        persistentHoverLayer.classList.add('persistent-hoverlayer');
        const shapelayer = layerAbove.querySelector('.shapelayer');
        if (shapelayer) {
            shapelayer.after(persistentHoverLayer);
        } else {
            layerAbove.appendChild(persistentHoverLayer);
        }
    }

    return persistentHoverLayer;
}

function getSvgElementRect(element) {
    if (!element || typeof element.getBoundingClientRect !== 'function') {
        return null;
    }
    const rect = element.getBoundingClientRect();
    if (!rect || rect.width === 0 || rect.height === 0) {
        return null;
    }
    return rect;
}

function rectsOverlap(first, second) {
    return (
        first.left < second.right
        && first.right > second.left
        && first.top < second.bottom
        && first.bottom > second.top
    );
}

function hoverBoxOverlapsBar(plotDiv, hoverBox) {
    const hoverRect = getSvgElementRect(hoverBox);
    if (!hoverRect) {
        return false;
    }
    const barElements = plotDiv.querySelectorAll('.barlayer .point path, .barlayer path');
    return Array.from(barElements).some(barElement => {
        const barRect = getSvgElementRect(barElement);
        return barRect && rectsOverlap(hoverRect, barRect);
    });
}

function setPinnedHoverBoxSolid(hoverBox) {
    hoverBox.dataset.pinState = 'solid';
    hoverBox.style.opacity = '1';
}

function configurePinnedHoverBox(plotDiv, hoverBox) {
    hoverBox.style.pointerEvents = 'all';
    hoverBox.style.cursor = 'pointer';

    if (hoverBoxOverlapsBar(plotDiv, hoverBox)) {
        hoverBox.dataset.pinState = 'translucent';
        hoverBox.style.opacity = '0.5';
    } else {
        setPinnedHoverBoxSolid(hoverBox);
    }

    hoverBox.addEventListener('click', function(event) {
        event.stopPropagation();
        if (this.dataset.pinState === 'translucent') {
            setPinnedHoverBoxSolid(this);
            return;
        }
        this.remove();
    });
}

function appendPinnedHoverBox(plotDiv, persistentHoverLayer, hoverBox, hoverBoxId) {
    const persistentHoverBox = hoverBox.cloneNode(true);
    persistentHoverBox.setAttribute('id', hoverBoxId);
    persistentHoverLayer.appendChild(persistentHoverBox);
    configurePinnedHoverBox(plotDiv, persistentHoverBox);
}

function pinHoverBox(plotDiv, persistentHoverLayer, eventData) {
    if (window.Plotly && eventData) {
        if (eventData.points) {
            Plotly.Fx.hover(plotDiv, eventData.points);
        } else if (typeof eventData.xpx === "number" && typeof eventData.ypx === "number") {
            Plotly.Fx.hover(plotDiv, eventData);
        }
    }
    setTimeout(function() {
        const hoverLayer = plotDiv.getElementsByClassName('hoverlayer')[0];
        if (hoverLayer) {
            const hoverBox = hoverLayer.getElementsByClassName('hovertext')[0];
            if (hoverBox) {
                const existingHoverBox = Array.from(persistentHoverLayer.children).find(child =>
                    child.textContent === hoverBox.textContent
                );

                if (existingHoverBox) {
                    existingHoverBox.remove();
                } else {
                    appendPinnedHoverBox(plotDiv, persistentHoverLayer, hoverBox, `hover-${Date.now()}`);
                }
            }
        }
        if (window.Plotly) {
            Plotly.Fx.unhover(plotDiv);
        }
    }, 0);
}

function setupPersistentHover(plotDiv) {
    if (!plotDiv || typeof plotDiv.on !== 'function') {
        return false;
    }
    const persistentHoverLayer = getPersistentHoverLayer(plotDiv);
    if (!persistentHoverLayer) {
        return false;
    }

    if (plotDiv.dataset.nativeClickBound !== 'true') {
        plotDiv.addEventListener('click', function(event) {
            const handledClickTs = plotDiv.__handledClickTs || 0;
            const now = Date.now();
            if (now - handledClickTs < 200) {
                return;
            }
            plotDiv.__handledClickTs = now;
            if (event && window.Plotly) {
                const rect = plotDiv.getBoundingClientRect();
                const xpx = event.clientX - rect.left;
                const ypx = event.clientY - rect.top;
                pinHoverBox(plotDiv, persistentHoverLayer, { xpx, ypx });
            }
        });
        plotDiv.dataset.nativeClickBound = 'true';
    }

    if (plotDiv.dataset.plotlyClickBound !== 'true') {
        plotDiv.on('plotly_click', function(eventData) {
            plotDiv.__handledClickTs = Date.now();
            pinHoverBox(plotDiv, persistentHoverLayer, eventData);
        });
        plotDiv.dataset.plotlyClickBound = 'true';
    }

    return true;
}


function bindPersistentHover(container, attemptsLeft = 30) {
    if (!container) {
        return;
    }
    const plotDivs = container.getElementsByClassName("js-plotly-plot");
    let needsRetry = false;
    if (plotDivs.length > 0) {
        Array.from(plotDivs).forEach(div => {
            if (div.dataset.persistentHoverBound === "true") {
                return;
            }
            if (setupPersistentHover(div)) {
                div.dataset.persistentHoverBound = "true";
            } else {
                needsRetry = true;
            }
        });
    } else {
        needsRetry = true;
    }
    if (needsRetry && attemptsLeft > 0) {
        setTimeout(function() {
            bindPersistentHover(container, attemptsLeft - 1);
        }, 100);
    }
}

function clearPersistentHoverLayer(plotDiv) {
    const persistentHoverLayer = getPersistentHoverLayer(plotDiv);
    if (persistentHoverLayer) {
        persistentHoverLayer.innerHTML = '';
    }
}

function cloneHoverBoxesToPersistent(plotDiv, persistentHoverLayer) {
    const hoverLayer = plotDiv.getElementsByClassName('hoverlayer')[0];
    if (!hoverLayer) {
        return;
    }
    const hoverBoxes = hoverLayer.getElementsByClassName('hovertext');
    Array.from(hoverBoxes).forEach(hoverBox => {
        const existingHoverBox = Array.from(persistentHoverLayer.children).find(child =>
            child.textContent === hoverBox.textContent
        );
        if (!existingHoverBox) {
            const hoverBoxId = `hover-${Date.now()}-${Math.floor(Math.random() * 10000)}`;
            appendPinnedHoverBox(plotDiv, persistentHoverLayer, hoverBox, hoverBoxId);
        }
    });
}

function getPeptidePlotDiv() {
    const tabs = document.querySelector('[data-tab-group="pepmap-peptides-tabs"]');
    if (tabs) {
        const activeButton = tabs.querySelector('[data-active="true"]') || tabs.querySelector('[data-tab-target]');
        if (activeButton) {
            const panelId = activeButton.getAttribute('data-tab-target');
            const panel = document.getElementById(panelId);
            if (panel) {
                return panel.querySelector('.js-plotly-plot');
            }
        }
    }
    const plotDivs = elements.peptidesPlotContent.getElementsByClassName('js-plotly-plot');
    return plotDivs[0] || null;
}

function collectPlotPoints(plotDiv) {
    const points = [];
    if (!plotDiv || !plotDiv._fullData) {
        return points;
    }
    plotDiv._fullData.forEach((trace, traceIndex) => {
        const pointCount = trace.x ? trace.x.length : 0;
        for (let i = 0; i < pointCount; i += 1) {
            points.push({ curveNumber: traceIndex, pointNumber: i });
        }
    });
    return points;
}

function showAllStats(plotDiv) {
    if (!plotDiv || !window.Plotly) {
        return;
    }
    const persistentHoverLayer = getPersistentHoverLayer(plotDiv);
    if (!persistentHoverLayer) {
        return;
    }

    clearPersistentHoverLayer(plotDiv);
    const points = collectPlotPoints(plotDiv);
    if (points.length === 0) {
        return;
    }

    let currentIndex = 0;
    const step = function() {
        if (currentIndex >= points.length) {
            Plotly.Fx.unhover(plotDiv);
            return;
        }
        Plotly.Fx.hover(plotDiv, [points[currentIndex]]);
        setTimeout(function() {
            cloneHoverBoxesToPersistent(plotDiv, persistentHoverLayer);
            Plotly.Fx.unhover(plotDiv);
            currentIndex += 1;
            setTimeout(step, 0);
        }, 0);
    };
    step();
}

function applyCrispEdges() {
    const colorbarGroups = document.querySelectorAll('g.colorbar');
    colorbarGroups.forEach(group => {
        const rects = group.querySelectorAll('rect');
        rects.forEach(rect => {
            rect.setAttribute('shape-rendering', 'crispEdges');
        });
    });
    const shapeGroups = document.querySelectorAll('g.shapelayer');
    shapeGroups.forEach(group => {
        const shapes = group.querySelectorAll('path');
        shapes.forEach(shape => {
            shape.setAttribute('shape-rendering', 'crispEdges');
        });
    });
}

function resizeVisiblePlots() {
    if (!window.Plotly) {
        return;
    }
    ['peptides_plot', 'features_plot'].forEach(id => {
        const container = document.getElementById(id);
        if (!container) {
            return;
        }
        const plotDivs = container.getElementsByClassName('js-plotly-plot');
        Array.from(plotDivs).forEach(div => {
            if (div.offsetParent !== null) {
                Plotly.Plots.resize(div);
            }
        });
    });
}

function updateFeatureVisibility() {
    const tabs = document.querySelector('[data-tab-group="pepmap-peptides-tabs"]');
    const panels = elements.featuresPlot.querySelectorAll('.feature-panel');
    if (!tabs || panels.length === 0) {
        panels.forEach(panel => { panel.style.display = 'block'; });
        return;
    }
    const activeButton = tabs.querySelector('[data-active="true"]');
    const activeIndex = activeButton ? activeButton.getAttribute('data-index') : null;
    panels.forEach(panel => {
        const panelIndex = panel.getAttribute('data-index');
        const isEmpty = panel.classList.contains('empty-feature-panel');
        panel.style.display = panelIndex === activeIndex && !isEmpty ? 'block' : 'none';
    });
}

function activateTab(tabGroup, button, shouldResize = true) {
    const buttons = tabGroup.querySelectorAll('[data-tab-target]');
    buttons.forEach(btn => {
        btn.classList.remove('bg-slate-300', 'text-slate-900', 'text-slate-700', 'text-red-600');
        btn.classList.add('bg-white');
        btn.classList.add(btn.dataset.tabTone === 'danger' ? 'text-red-600' : 'text-slate-700');
        btn.removeAttribute('data-active');
    });
    button.classList.remove('bg-white', 'text-slate-700', 'text-slate-900', 'text-red-600');
    button.classList.add('bg-slate-300');
    button.classList.add(button.dataset.tabTone === 'danger' ? 'text-red-600' : 'text-slate-900');
    button.setAttribute('data-active', 'true');

    const panels = tabGroup.querySelectorAll('.tab-panel');
    panels.forEach(panel => panel.classList.add('hidden'));
    const targetId = button.getAttribute('data-tab-target');
    const targetPanel = document.getElementById(targetId);
    if (targetPanel) {
        targetPanel.classList.remove('hidden');
        bindPersistentHover(targetPanel);
    }
    updateFeatureVisibility();
    if (shouldResize) {
        resizeVisiblePlots();
    }
}

function syncPlotTabs() {
    const tabGroup = document.querySelector('[data-tab-group="pepmap-peptides-tabs"]');
    if (!tabGroup) {
        updateFeatureVisibility();
        return;
    }
    const buttons = tabGroup.querySelectorAll('[data-tab-target]');
    buttons.forEach(button => {
        button.addEventListener('click', function() {
            activateTab(tabGroup, button);
        });
    });
    const activeButton = tabGroup.querySelector('[data-active="true"]') || buttons[0];
    if (activeButton) {
        activateTab(tabGroup, activeButton, false);
    }
}

function moveSummaryToTable(attemptsLeft) {
    const summary = document.querySelector('#peptides_plot_content .peptide-summary');
    if (!summary || !elements.peptidesTable) {
        if (attemptsLeft > 0) {
            setTimeout(function() {
                moveSummaryToTable(attemptsLeft - 1);
            }, 100);
        }
        return;
    }
    elements.peptidesTable.innerHTML = '';
    elements.peptidesTable.appendChild(summary);
    setHidden(elements.peptidesTable, false);
    if (elements.copyTableButton) {
        setHidden(elements.copyTableButton, false);
    }
}

function tableToTsv(table) {
    const rows = Array.from(table.querySelectorAll('tr'));
    return rows.map(row => {
        const cells = Array.from(row.querySelectorAll('th,td'));
        return cells.map(cell => cell.textContent.trim()).join('\t');
    }).join('\n');
}

function buildTablesCopyPayload() {
    const tables = Array.from(elements.peptidesTable?.querySelectorAll('table') || []);
    if (!tables.length) {
        return '';
    }
    const parts = [];
    const summary = elements.peptidesTable?.querySelector('.peptide-summary');
    const summaryTitle = summary
        ?.querySelector('.text-sm.font-semibold.text-center')
        ?.textContent
        ?.trim();
    if (summaryTitle) {
        parts.push(`# ${summaryTitle}`);
    }

    tables.forEach(table => {
        const wrapper = table.closest('.overflow-x-auto');
        let label = '';
        let cursor = wrapper ? wrapper.previousElementSibling : null;
        while (cursor) {
            if (
                cursor.classList.contains('text-sm')
                && cursor.classList.contains('font-semibold')
                && !cursor.classList.contains('text-center')
            ) {
                label = cursor.textContent?.trim() || '';
                break;
            }
            cursor = cursor.previousElementSibling;
        }
        if (label) {
            parts.push(`## ${label}`);
        }
        parts.push(tableToTsv(table));
    });
    return parts.join('\n\n');
}

async function requestPlot(url, formData) {
    const response = await fetch(url, { method: 'POST', body: formData });
    const text = await response.text();
    if (response.ok) {
        return { ok: true, data: text };
    }
    let errorMessage = 'Unknown Error';
    try {
        const data = JSON.parse(text);
        if (data.error) {
            errorMessage = data.error;
        }
    } catch (error) {
        if (text) {
            errorMessage = text;
        }
    }
    return { ok: false, error: errorMessage };
}

async function fetchPlots() {
    setHidden(elements.peptidesPlot, false);
    setHidden(elements.featuresPlot, false);
    setStatsControlsVisible(false);
    elements.peptidesPlotContent.innerHTML = renderLoader();
    elements.featuresPlot.innerHTML = renderLoader();
    elements.peptidesTable.innerHTML = '';
    setHidden(elements.peptidesTable, true);
    if (elements.copyTableButton) {
        setHidden(elements.copyTableButton, true);
    }

    const geneInputValue = elements.geneEntryInput.value.trim();
    if (geneInputValue) {
        await addGeneFromInput();
    }
    updateMapButtonState();
    if (!geneEntries.length) {
        return;
    }

    const formData = new FormData();
    if (sessionId) {
        formData.append('session_id', sessionId);
    }
    formData.append('search_input', elements.searchInput.value);
    formData.append('search_labels', elements.searchLabels.value);
    formData.append('proteotypic_checkbox', elements.proteotypicCheckbox?.checked ? 'true' : 'false');
    formData.append('charge_state_mode', elements.chargeStateMode?.value || 'all');
    formData.append('q_value_cutoff', elements.qValueCutoff?.value || '0.01');
    formData.append('sample_name_cleanup', elements.sampleNameCleanup.value);
    formData.append('sample_name_custom_pattern', elements.sampleNameCustomPattern.value);
    const selectedRuns = getSelectedRunsPayload();
    if (selectedRuns) {
        formData.append('selected_runs', JSON.stringify(selectedRuns));
    }
    formData.append('custom_title', '');
    formData.append('summary_mode', elements.summaryMode?.value || 'per_sample');

    setButtonLoading(elements.submitButton, true, 'Mapping...');
    const [peptidesResult, featuresResult] = await Promise.all([
        requestPlot('/plot_peptides', formData),
        requestPlot('/plot_features', formData)
    ]);
    setButtonLoading(elements.submitButton, false);

    if (peptidesResult.ok) {
        setHtmlAndRunScripts(elements.peptidesPlotContent, peptidesResult.data);
        applyCrispEdges();
        bindPersistentHover(elements.peptidesPlotContent);
        setStatsVisibility(true);
        syncPlotTabs();
        moveSummaryToTable(10);
    } else {
        elements.peptidesPlotContent.innerHTML = `<div class="text-red-600">Failed to load peptides plot. (${peptidesResult.error})</div>`;
        setStatsVisibility(false);
    }

    if (featuresResult.ok) {
        setHtmlAndRunScripts(elements.featuresPlot, featuresResult.data);
        if (!elements.featuresPlot.textContent.trim()) {
            setHidden(elements.featuresPlot, true);
            return;
        }
        setHidden(elements.featuresPlot, false);
        bindPersistentHover(elements.featuresPlot);
        syncPlotTabs();
        moveSummaryToTable(10);
    } else {
        elements.featuresPlot.innerHTML = `<div class="text-red-600">Failed to load features plot. (${featuresResult.error})</div>`;
    }
}
