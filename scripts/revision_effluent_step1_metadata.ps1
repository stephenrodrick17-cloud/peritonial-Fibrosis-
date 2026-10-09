param(
    [string]$Root = (Join-Path (Get-Location) 'results\revision\effluent_validation')
)

$ErrorActionPreference = 'Stop'

function Read-GeoSoftMetadata {
    param(
        [string]$Accession,
        [string]$Path
    )

    $fileStream = [System.IO.File]::OpenRead((Resolve-Path -LiteralPath $Path))
    $gzipStream = [System.IO.Compression.GzipStream]::new(
        $fileStream,
        [System.IO.Compression.CompressionMode]::Decompress
    )
    $reader = [System.IO.StreamReader]::new($gzipStream)
    $series = [ordered]@{
        accession = $Accession
        title = ''
        pubmed_id = ''
        contributors = [System.Collections.Generic.List[string]]::new()
        processing = [System.Collections.Generic.List[string]]::new()
        relations = [System.Collections.Generic.List[string]]::new()
    }
    $samples = [System.Collections.Generic.List[object]]::new()
    $sample = $null
    $inSampleTable = $false

    try {
        while (($line = $reader.ReadLine()) -ne $null) {
            if ($line -match '^\^SAMPLE\s*=\s*(\S+)') {
                if ($null -ne $sample) {
                    $samples.Add($sample)
                }
                $sample = [pscustomobject]@{
                    accession = $Accession
                    sample_id = $Matches[1]
                    title = ''
                    source_name = ''
                    characteristics = [System.Collections.Generic.List[string]]::new()
                    description = ''
                    platform_id = ''
                    supplementary_files = [System.Collections.Generic.List[string]]::new()
                    patient_id_reported = ''
                    group_from_metadata = ''
                }
                $inSampleTable = $false
                continue
            }

            if ($line -match '^!Sample_table_begin') {
                break
            }

            if ($line -match '^!Series_title\s*=\s*(.*)$') {
                $series.title = $Matches[1]
                continue
            }
            if ($line -match '^!Series_pubmed_id\s*=\s*(.*)$') {
                $series.pubmed_id = $Matches[1]
                continue
            }
            if ($line -match '^!Series_contributor\s*=\s*(.*)$') {
                $series.contributors.Add($Matches[1])
                continue
            }
            if ($line -match '^!Series_overall_design\s*=\s*(.*)$' -or
                $line -match '^!Series_summary\s*=\s*(.*)$') {
                $series.processing.Add($Matches[1])
                continue
            }
            if ($line -match '^!Series_relation\s*=\s*(.*)$') {
                $series.relations.Add($Matches[1])
                continue
            }
            if ($null -eq $sample) {
                continue
            }
            if ($line -match '^!Sample_title\s*=\s*(.*)$') {
                $sample.title = $Matches[1]
            } elseif ($line -match '^!Sample_source_name_ch1\s*=\s*(.*)$') {
                $sample.source_name = $Matches[1]
            } elseif ($line -match '^!Sample_characteristics_ch1\s*=\s*(.*)$') {
                $sample.characteristics.Add($Matches[1])
            } elseif ($line -match '^!Sample_description\s*=\s*(.*)$') {
                $sample.description = $Matches[1]
            } elseif ($line -match '^!Sample_platform_id\s*=\s*(.*)$') {
                $sample.platform_id = $Matches[1]
            } elseif ($line -match '^!Sample_supplementary_file\s*=\s*(.*)$') {
                $sample.supplementary_files.Add($Matches[1])
            }
        }
        if ($null -ne $sample) {
            $samples.Add($sample)
        }
    } finally {
        $reader.Dispose()
        $gzipStream.Dispose()
        $fileStream.Dispose()
    }

    foreach ($row in $samples) {
        $characteristics = $row.characteristics -join '; '
        $allLabels = "$($row.title); $($row.source_name); $characteristics; $($row.description)"
        if ($allLabels -match '(?i)(?:^|[;\s])(?:patient|subject|donor|individual)(?:[\s_-]*id)?\s*:\s*([^;]+)') {
            $row.patient_id_reported = $Matches[1].Trim()
        } else {
            $row.patient_id_reported = 'not reported in sample metadata'
        }
        switch -Regex ($Accession) {
            '^GSE248762$' {
                if ($allLabels -match '(?i)LV[_ -]?NOT[_ -]?UF|long vintage without (?:ultrafiltration failure|UF)') {
                    $row.group_from_metadata = 'LV_NOT_UF'
                } elseif ($allLabels -match '(?i)LV[_ -]?UF|long vintage with (?:ultrafiltration failure|UF)') {
                    $row.group_from_metadata = 'LV_UF'
                } elseif ($allLabels -match '(?i)\bSV\b|short vintage') {
                    $row.group_from_metadata = 'SV'
                }
            }
            '^GSE130888$' {
                if ($allLabels -match '(?i)normal control|normal peritoneal|cells dissociated from normal peritoneum|\bHC\b') {
                    $row.group_from_metadata = 'Normal control'
                } elseif ($allLabels -match '(?i)long[- ]term PD') {
                    $row.group_from_metadata = 'Long-term PD'
                } elseif ($allLabels -match '(?i)short[- ]term PD') {
                    $row.group_from_metadata = 'Short-term PD'
                }
            }
            '^GSE92455$' {
                if ($allLabels -match '(?i)omentum') {
                    $row.group_from_metadata = 'Omentum'
                } elseif ($allLabels -match '(?i)non[- ]?epithelioid') {
                    $row.group_from_metadata = 'Non-epithelioid'
                } elseif ($allLabels -match '(?i)epithelioid') {
                    $row.group_from_metadata = 'Epithelioid'
                }
            }
        }
        $row.characteristics = $characteristics
        $row.supplementary_files = $row.supplementary_files -join '; '
    }

    [pscustomobject]@{
        Series = $series
        Samples = $samples
    }
}

function Split-GsmTabLine {
    param([string]$Line)

    $values = [System.Collections.Generic.List[string]]::new()
    $value = [System.Text.StringBuilder]::new()
    $quoted = $false
    for ($i = 0; $i -lt $Line.Length; $i++) {
        $character = $Line[$i]
        if ($character -eq '"') {
            if ($quoted -and $i + 1 -lt $Line.Length -and $Line[$i + 1] -eq '"') {
                [void]$value.Append('"')
                $i++
            } else {
                $quoted = -not $quoted
            }
        } elseif ($character -eq "`t" -and -not $quoted) {
            $values.Add($value.ToString())
            [void]$value.Clear()
        } else {
            [void]$value.Append($character)
        }
    }
    $values.Add($value.ToString())
    return ,$values.ToArray()
}

function Read-GeoMatrixSampleMetadata {
    param(
        [string]$Accession,
        [string]$Path,
        [string]$Platform
    )

    $fileStream = [System.IO.File]::OpenRead((Resolve-Path -LiteralPath $Path))
    $gzipStream = [System.IO.Compression.GzipStream]::new(
        $fileStream,
        [System.IO.Compression.CompressionMode]::Decompress
    )
    $reader = [System.IO.StreamReader]::new($gzipStream)
    $vectors = @{}
    try {
        while (($line = $reader.ReadLine()) -ne $null) {
            if ($line -match '^!series_matrix_table_begin') {
                break
            }
            if ($line -notmatch '^!(Sample_[^\t]+)\t') {
                continue
            }
            $fields = Split-GsmTabLine -Line $line
            $key = $Matches[1]
            $vectors[$key] = @($fields | Select-Object -Skip 1)
        }
    } finally {
        $reader.Dispose()
        $gzipStream.Dispose()
        $fileStream.Dispose()
    }

    if (-not $vectors.ContainsKey('Sample_geo_accession') -or
        -not $vectors.ContainsKey('Sample_title')) {
        throw "GEO matrix header lacks required sample identity fields: $Path"
    }
    $ids = @($vectors['Sample_geo_accession'])
    $rows = [System.Collections.Generic.List[object]]::new()
    for ($i = 0; $i -lt $ids.Count; $i++) {
        $get = {
            param([string]$Key)
            if ($vectors.ContainsKey($Key) -and $i -lt $vectors[$Key].Count) {
                return [string]$vectors[$Key][$i]
            }
            return ''
        }
        $title = & $get 'Sample_title'
        $source = & $get 'Sample_source_name_ch1'
        $description = & $get 'Sample_description'
        $characteristics = @()
        foreach ($key in ($vectors.Keys | Where-Object { $_ -like 'Sample_characteristics_ch*' } | Sort-Object)) {
            if ($i -lt $vectors[$key].Count) {
                $characteristics += [string]$vectors[$key][$i]
            }
        }
        $allLabels = "$title; $source; $description; $($characteristics -join '; ')"
        $group = ''
        if ($Platform -eq 'GPL6480') {
            if ($title -match '(?i)^Ex Vivo Control Pool') {
                $group = 'Pooled omentum reference'
            } elseif ($title -match '(?i)Non[- ]Epithelioid') {
                $group = 'Non-epithelioid effluent HPMC'
            } elseif ($title -match '(?i)Epithelioid phenotype') {
                $group = 'Epithelioid effluent HPMC'
            }
        } elseif ($Platform -eq 'GPL6848' -and $title -match '(?i)In Vitro HPMCs Treated') {
            if ($source -match '(?i)_Control$') {
                $group = 'In vitro control'
            } elseif ($source -match '(?i)_Treated$') {
                $group = 'In vitro treated'
            }
        }
        $patientId = 'not reported in sample metadata'
        if ($allLabels -match '(?i)HLP\d+') {
            $patientId = $Matches[0].ToUpperInvariant()
        } elseif ($allLabels -match '(?i)(?:^|[;\s])(?:patient|subject|donor|individual)(?:[\s_-]*id)?\s*:\s*([^;]+)') {
            $patientId = $Matches[1].Trim()
        }
        $rows.Add([pscustomobject]@{
            accession = $Accession
            platform = $Platform
            sample_id = $ids[$i]
            title = $title
            source_name = $source
            patient_id_reported = $patientId
            group_from_metadata = $group
            characteristics = $characteristics -join '; '
            description = $description
        })
    }
    return $rows
}

$accessions = @('GSE248762', 'GSE130888', 'GSE92455', 'GSE125498')
$seriesRows = [System.Collections.Generic.List[object]]::new()
foreach ($accession in $accessions) {
    $dir = Join-Path $Root $accession
    $softPath = Join-Path $dir "${accession}_family.soft.gz"
    if (-not (Test-Path -LiteralPath $softPath)) {
        throw "Required GEO family SOFT file is missing: $softPath"
    }
    $parsed = Read-GeoSoftMetadata -Accession $accession -Path $softPath
    if ($parsed.Samples.Count -gt 0) {
        $parsed.Samples | Export-Csv -LiteralPath (Join-Path $dir 'sample_metadata.csv') -NoTypeInformation -Encoding UTF8
    }
    $seriesRows.Add([pscustomobject]@{
        accession = $accession
        title = $parsed.Series.title
        pubmed_id = $parsed.Series.pubmed_id
        contributors = $parsed.Series.contributors -join '; '
        processing_or_design = $parsed.Series.processing -join ' | '
        relations = $parsed.Series.relations -join '; '
        family_soft_bytes = (Get-Item -LiteralPath $softPath).Length
        sample_records = $parsed.Samples.Count
    })
}
$seriesRows | Export-Csv -LiteralPath (Join-Path $Root 'series_metadata_audit.csv') -NoTypeInformation -Encoding UTF8

$matrixSamples = [System.Collections.Generic.List[object]]::new()
foreach ($platform in @('GPL6480', 'GPL6848')) {
    $matrixPath = Join-Path (Join-Path $Root 'GSE92455') "GSE92455-${platform}_series_matrix.txt.gz"
    if (-not (Test-Path -LiteralPath $matrixPath)) {
        throw "Required processed matrix is missing: $matrixPath"
    }
    foreach ($row in (Read-GeoMatrixSampleMetadata -Accession 'GSE92455' -Path $matrixPath -Platform $platform)) {
        $matrixSamples.Add($row)
    }
}
$matrixSamples |
    Sort-Object platform, sample_id -Unique |
    Export-Csv -LiteralPath (Join-Path (Join-Path $Root 'GSE92455') 'sample_metadata.csv') -NoTypeInformation -Encoding UTF8

$groupSummaries = [System.Collections.Generic.List[object]]::new()
foreach ($accession in @('GSE248762', 'GSE130888', 'GSE92455')) {
    $metadataPath = Join-Path (Join-Path $Root $accession) 'sample_metadata.csv'
    $rows = @(Import-Csv -LiteralPath $metadataPath)
    foreach ($group in ($rows | Group-Object group_from_metadata)) {
        $reportedPatients = @(
            $group.Group |
                Where-Object { $_.patient_id_reported -ne 'not reported in sample metadata' } |
                Select-Object -ExpandProperty patient_id_reported -Unique
        )
        $groupSummaries.Add([pscustomobject]@{
            accession = $accession
            platform = (($group.Group | Select-Object -First 1).platform)
            group = $group.Name
            sample_records_or_arrays = $group.Count
            distinct_patient_codes_present_in_GEO = $reportedPatients.Count
            patient_ids_complete = ($reportedPatients.Count -eq $group.Count)
        })
    }
}
$groupSummaries |
    Export-Csv -LiteralPath (Join-Path $Root 'sample_group_counts.csv') -NoTypeInformation -Encoding UTF8

$fileAudit = [System.Collections.Generic.List[object]]::new()
foreach ($accession in @('GSE248762', 'GSE130888', 'GSE92455', 'GSE125498')) {
    $dir = Join-Path $Root $accession
    foreach ($file in (Get-ChildItem -LiteralPath $dir -File)) {
        $state = switch -Regex ($file.Name) {
            '_family\.soft\.gz$' { 'GEO family metadata; no expression table'; break }
            '^filelist\.txt$' { 'GEO supplementary-file manifest'; break }
            '^GPL6480\.annot\.gz$' {
                'GEO platform annotation; probe-to-gene identifiers, not expression values'
                break
            }
            '_series_matrix\.txt\.gz$' {
                if ($accession -eq 'GSE92455') {
                    'Processed expression matrix; GEO states background subtraction, LOESS/LOWESS and quantile normalization'
                } else {
                    'Sample metadata only; no series expression table'
                }
                break
            }
            default { 'Metadata audit output or unclassified file' }
        }
        if ($file.Name -match '_family\.soft\.gz$' -or
            $file.Name -eq 'filelist.txt' -or
            $file.Name -eq 'GPL6480.annot.gz' -or
            $file.Name -match '_series_matrix\.txt\.gz$') {
            $fileAudit.Add([pscustomobject]@{
                accession = $accession
                file_name = $file.Name
                expected_bytes_from_GEO_manifest = ''
                local_bytes = $file.Length
                download_status = 'downloaded'
                processing_state = $state
            })
        }
    }
}
foreach ($accession in @('GSE248762', 'GSE130888', 'GSE92455')) {
    $dir = Join-Path $Root $accession
    $manifest = Join-Path $dir 'filelist.txt'
    if (-not (Test-Path -LiteralPath $manifest)) {
        continue
    }
    foreach ($line in (Get-Content -LiteralPath $manifest | Select-Object -Skip 1)) {
        $fields = $line -split "`t"
        if ($fields.Count -lt 4 -or $fields[0] -ne 'Archive') {
            continue
        }
        $name = $fields[1]
        $expectedBytes = [long]($fields[3] -replace '[^0-9]', '')
        $localArchive = Join-Path $dir $name
        $localBytes = ''
        $status = 'listed by GEO; archive not downloaded'
        if (Test-Path -LiteralPath $localArchive) {
            $localBytes = (Get-Item -LiteralPath $localArchive).Length
            if ([long]$localBytes -eq $expectedBytes) {
                $status = 'downloaded complete'
            } else {
                $status = 'incomplete download; not usable'
            }
        }
        $fileAudit.Add([pscustomobject]@{
            accession = $accession
            file_name = $name
            expected_bytes_from_GEO_manifest = $expectedBytes
            local_bytes = $localBytes
            download_status = $status
            processing_state = if ($accession -eq 'GSE92455') {
                'RAW archive listed by GEO; separate processed matrices downloaded'
            } else {
                'RAW archive contains 10x feature-barcode count matrices'
            }
        })
    }
}
$fileAudit |
    Export-Csv -LiteralPath (Join-Path $Root 'step1_file_download_audit.csv') -NoTypeInformation -Encoding UTF8
