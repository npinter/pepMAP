window.addEventListener('DOMContentLoaded', function() {
    resetFormState();
    if (elements.timer && !elements.timer.classList.contains('hidden')) {
        startTimer(sessionTime, elements.timer);
    }

    const params = new URLSearchParams(window.location.search);
    const providedSessionId = params.get('session_id');
    if (providedSessionId) {
        updateSessionShare(providedSessionId);
        fetch(`/session_info?session_id=${encodeURIComponent(providedSessionId)}`)
            .then(response => {
                if (!response.ok) {
                    return response.json().then(data => {
                        throw new Error(data.error || 'Failed to load session.');
                    });
                }
                return response.json();
            })
            .then(data => {
                applySessionData(data);
                const genesParam = params.get('genes');
                const labelsParam = params.get('labels');

                if (genesParam) {
                    const genes = parseGeneEntries(genesParam);
                    let labels = [];
                    if (labelsParam) {
                        try {
                            const parsed = JSON.parse(labelsParam);
                            if (Array.isArray(parsed)) {
                                labels = parsed.map(item => String(item));
                            } else {
                                labels = parseGeneEntries(labelsParam);
                            }
                        } catch (e) {
                            labels = parseGeneEntries(labelsParam);
                        }
                    }
                    if (genes.length) {
                        if (labels.length) {
                            addGeneEntriesWithLabels(genes, labels);
                        } else {
                            addGeneEntries(genes);
                        }
                    }
                }
                resolveGeneEntries(0);
                if (elements.organism.value === 'CUSTOM') {
                    setHidden(elements.customOrganismContainer, false);
                }
            })
            .catch(error => {
                alert(error.message);
                resetFormState();
                clearSessionUrlParams();
            });
    }

    elements.geneEntryInput.addEventListener('keypress', function(event) {
        if (event.key === 'Enter') {
            event.preventDefault();
            if (!elements.geneEntryInput.value.trim() && geneEntries.length > 0) {
                cancelAutocompleteRequests();
                clearAutocompleteSuggestions();
                fetchPlots();
                resetTimer();
                return;
            }
            addGeneFromInput();
        }
    });

    elements.uploadButton.addEventListener('click', function() {
        const formData = new FormData();
        setButtonLoading(elements.uploadButton, true, 'Uploading...');

        formData.append('report_file', elements.reportFile.files[0]);
        formData.append('fasta_file', elements.fastaFile.files[0]);
        if (elements.organism.value === 'CUSTOM') {
            formData.append('organism', 'CUSTOM');
            formData.append('custom_organism', elements.customOrganismInput.value);
        } else {
            formData.append('organism', elements.organism.value);
        }
        formData.append('custom_features_file', elements.customFeaturesFile.files[0]);
        formData.append('custom_features_label', elements.customFeaturesLabel.value);

        fetch('/upload', { method: 'POST', body: formData })
            .then(response => {
                if (!response.ok) {
                    return response.json().then(data => {
                        throw new Error(data.error || 'Upload failed');
                    });
                }
                return response.json();
            })
            .then(data => {
                sessionId = data.session_id;
                clearSessionUrlParams();
                updateSessionShare(sessionId);
                elements.reportFilename.value = elements.reportFile.files[0].name;
                elements.fastaFilename.value = elements.fastaFile.files[0].name;
                setHidden(elements.reportFilename, false);
                setHidden(elements.fastaFilename, false);
                setHidden(elements.reportFile, true);
                setHidden(elements.fastaFile, true);
                elements.organism.disabled = true;
                setHidden(elements.uploadButton, true);
                setHidden(elements.searchContainer, false);
                setHidden(elements.mapContainer, false);
                setHidden(elements.newReportButton, false);
                setHidden(elements.mappingBlock, false);
                setHidden(elements.timer, false);
                setAdvancedFieldsVisible(true);
                setAvailableRuns(data.runs || [], true);
                updateMapButtonState();
                resetTimer();
                hasCustomFeatures = Boolean(data.custom_features_uploaded);
                customFeaturesLabel = data.custom_features_label || '';
                const uploadedCustomFeaturesFilename = elements.customFeaturesFile?.files?.[0]?.name || '';
                updateCustomFeaturesUI({
                    hasCustomFeatures,
                    filename: uploadedCustomFeaturesFilename,
                    label: customFeaturesLabel,
                    reportReady: true
                });
                updateReportTypeBadge(data.report_type || '', data.report_logo || '');
                resolveGeneEntries(0);
            })
            .catch(error => {
                alert(error.message);
            })
            .finally(() => {
                setButtonLoading(elements.uploadButton, false);
            });
    });

    elements.newReportButton.addEventListener('click', function() {
        const formData = new FormData();
        if (sessionId) {
            formData.append('session_id', sessionId);
        }
        fetch('/flush', { method: 'POST', body: formData })
            .then(() => {
                resetFormState();
                clearSessionUrlParams();
                hasCustomFeatures = false;
                customFeaturesLabel = '';
                sessionId = null;
                updateSessionShare(null);
            });
    });

    elements.customFeaturesFile.addEventListener('change', function() {
        const hasFile = this.files && this.files.length > 0;
        setHidden(elements.customFeaturesLabelContainer, !hasFile);
        if (hasFile && !elements.customFeaturesLabel.value.trim()) {
            elements.customFeaturesLabel.value = getFilenameStem(this.files[0].name);
        }
        if (!hasFile) {
            elements.customFeaturesLabel.value = '';
        }
    });

    elements.submitButton.addEventListener('click', function(event) {
        event.preventDefault();
        fetchPlots();
        resetTimer();
    });

    elements.addGeneButton.addEventListener('click', function(event) {
        event.preventDefault();
        addGeneFromInput();
    });

    elements.geneEntryInput.addEventListener('input', function() {
        const tokens = parseGeneEntries(elements.geneEntryInput.value);
        const query = tokens.length ? tokens[tokens.length - 1] : '';
        if (!query) {
            cancelAutocompleteRequests();
            clearAutocompleteSuggestions();
            return;
        }
        scheduleAutocomplete(query);
    });

    elements.sampleNameCleanup.addEventListener('change', function() {
        const isCustom = elements.sampleNameCleanup.value === 'custom';
        setHidden(elements.sampleNameCustomContainer, !isCustom);
        updateSampleNamePreview();
        renderSamplePicker(true);
    });

    if (elements.sampleNameCustomPattern) {
        elements.sampleNameCustomPattern.addEventListener('input', function() {
            scheduleSampleCleanupUpdate();
        });
    }

    if (elements.sampleSelectAll) {
        elements.sampleSelectAll.addEventListener('click', function() {
            selectedRuns = new Set(cleanedRuns);
            selectionClearedManually = false;
            renderSamplePicker(false);
        });
    }
    if (elements.sampleClearAll) {
        elements.sampleClearAll.addEventListener('click', function() {
            selectedRuns = new Set();
            selectionClearedManually = true;
            renderSamplePicker(false);
        });
    }

    elements.organism.addEventListener('change', function() {
        const isCustom = elements.organism.value === 'CUSTOM';
        setHidden(elements.customOrganismContainer, !isCustom);
        if (!isCustom && elements.customOrganismInput) {
            elements.customOrganismInput.value = '';
        }
    });

    if (elements.qValueCutoff && !elements.qValueCutoff.value) {
        elements.qValueCutoff.value = '0.01';
    }

    elements.showAllStatsButton.addEventListener('click', function(event) {
        event.preventDefault();
        showAllStats(getPeptidePlotDiv());
    });

    elements.clearStatsButton.addEventListener('click', function(event) {
        event.preventDefault();
        clearPersistentHoverLayer(getPeptidePlotDiv());
    });

    if (elements.copyTableButton) {
        elements.copyTableButton.addEventListener('click', async function() {
            const payload = buildTablesCopyPayload();
            if (!payload) {
                return;
            }
            try {
                await navigator.clipboard.writeText(payload);
                elements.copyTableButton.classList.add('text-emerald-600');
                setTimeout(() => elements.copyTableButton.classList.remove('text-emerald-600'), 1200);
            } catch (error) {
                alert('Copy failed. Please select and copy the table manually.');
            }
        });
    }

    const handleSessionShareCopy = async function() {
        if (!elements.sessionShareLink) {
            return;
        }
        const url = elements.sessionShareLink.dataset.url || elements.sessionShareLink.textContent;
        if (!url) {
            return;
        }
        try {
            await navigator.clipboard.writeText(url);
            if (elements.sessionShareHint) {
                elements.sessionShareHint.textContent = 'Copied to clipboard.';
            }
            elements.sessionShareContainer?.classList.add('border-emerald-300', 'bg-emerald-50');
            setTimeout(() => {
                elements.sessionShareContainer?.classList.remove('border-emerald-300', 'bg-emerald-50');
                if (elements.sessionShareHint) {
                    elements.sessionShareHint.textContent = 'Click the link to copy the full URL.';
                }
            }, 1200);
        } catch (error) {
            if (elements.sessionShareHint) {
                elements.sessionShareHint.textContent = 'Copy failed. Press Ctrl+C to copy.';
            }
        }
    };

    if (elements.sessionShareLink) {
        elements.sessionShareLink.addEventListener('click', handleSessionShareCopy);
    }
    if (elements.sessionShareCopy) {
        elements.sessionShareCopy.addEventListener('click', handleSessionShareCopy);
    }

    if (elements.summaryMode) {
        elements.summaryMode.addEventListener('change', function() {
            if (!sessionId) {
                return;
            }
            fetchPlots();
            resetTimer();
        });
    }
});
