SAMPLES = [
    "wat_dose_response",
]

rule all:
    input:
        expand("data/input/{sample}_filtered_matches.idXML", sample=SAMPLES),
        expand("data/input/{sample}.mzML", sample=SAMPLES),

rule convert:
    input:
        "data/input/{sample}.raw"
    output:
        "data/pipeline/{sample}/{sample}.mzML"
    shell:
        """
        docker run -it --rm --volume $PWD:/workspace ghcr.io/radusuciu/yan_etal_natcomm_2024 \
        mono /opt/ThermoRawFileParser/ThermoRawFileParser.exe \
            --input /workspace/{input} \
            --output_file /workspace/{output} \
            --format 2 \
            --metadata 1 \
            --noiseData
        """

rule search:
    input:
        "data/pipeline/{sample}/{sample}.mzML"
    output:
        "data/pipeline/{sample}/{sample}.pep.xml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        """
        docker run -it --rm --volume $PWD:/workspace ghcr.io/radusuciu/yan_etal_natcomm_2024 \
        comet \
            -P/workspace/data/params/comet_tmt.params \
            -D/workspace/data/fasta/uniprot_mouse_2020_10_10_long-ddhd1_contam_revcat.fasta \
            "/workspace/{input}"
        """

rule add_mod_descriptions:
    input:
        "data/pipeline/{sample}/{sample}.pep.xml"
    output:
        "data/pipeline/{sample}/{sample}.fixed_mods.pep.xml"
    run:
        import pathlib
        from utils import add_mod_descriptions_to_comet_output
        add_mod_descriptions_to_comet_output(pathlib.Path(input[0]))

rule convert_ids:
    input:
        "data/pipeline/{sample}/{sample}.fixed_mods.pep.xml"
    output:
        "data/pipeline/{sample}/{sample}.idXML"
    shell:
        """
        docker run -it --rm --volume $PWD:/workspace --volume $PWD/data/params/unimod.xml:/opt/OpenMS/share/OpenMS/CHEMISTRY/unimod.xml ghcr.io/radusuciu/yan_etal_natcomm_2024 \
        IDFileConverter \
            -in "/workspace/{input}" \
            -out "/workspace/{output}"
        """

rule merge_ids:
    input:
        "data/pipeline/{sample}/{sample}.idXML"
    output:
        "data/pipeline/{sample}/search_results.idXML"
    shell:
        """
        docker run -it --rm --volume $PWD:/workspace --volume $PWD/data/params/unimod.xml:/opt/OpenMS/share/OpenMS/CHEMISTRY/unimod.xml ghcr.io/radusuciu/yan_etal_natcomm_2024 \
        IDMerger \
            -in "/workspace/{input}" \
            -out "/workspace/{output}" \
            -annotate_file_origin "true"
        """

rule index_peptides:
    input:
        "data/pipeline/{sample}/search_results.idXML"
    output:
        "data/pipeline/{sample}/peptides_indexed.idXML"
    shell:
        """
        docker run -it --rm --volume $PWD:/workspace --volume $PWD/data/params/unimod.xml:/opt/OpenMS/share/OpenMS/CHEMISTRY/unimod.xml ghcr.io/radusuciu/yan_etal_natcomm_2024 \
        PeptideIndexer \
            -in "/workspace/{input}" \
            -out "/workspace/{output}" \
            -fasta "/workspace/data/fasta/uniprot_mouse_2020_10_10_long-ddhd1_contam_revcat.fasta" \
            -decoy_string "Reverse_" \
            -decoy_string_position "prefix" \
            -enzyme:specificity "semi" \
            -missing_decoy_action "warn" \
            -write_protein_description
        """

rule extract_psm_features:
    input:
        "data/pipeline/{sample}/peptides_indexed.idXML"
    output:
        "data/pipeline/{sample}/psm_features.idXML"
    shell:
        """
        docker run -it --rm --volume $PWD:/workspace --volume $PWD/data/params/unimod.xml:/opt/OpenMS/share/OpenMS/CHEMISTRY/unimod.xml ghcr.io/radusuciu/yan_etal_natcomm_2024 \
        PSMFeatureExtractor \
            -in "/workspace/{input}" \
            -out "/workspace/{output}"
        """

rule percolator:
    input:
        "data/pipeline/{sample}/psm_features.idXML"
    output:
        "data/pipeline/{sample}/percolated.idXML"
    shell:
        """
        docker run -it --rm --volume $PWD:/workspace --volume $PWD/data/params/unimod.xml:/opt/OpenMS/share/OpenMS/CHEMISTRY/unimod.xml ghcr.io/radusuciu/yan_etal_natcomm_2024 \
        PercolatorAdapter \
            -in "/workspace/{input}" \
            -out "/workspace/{output}" \
            -percolator_executable /usr/local/bin/percolator
        """

rule filter_ids:
    input:
        "data/pipeline/{sample}/percolated.idXML"
    output:
        "data/pipeline/{sample}/filtered_matches.idXML"
    shell:
        """
        docker run -it --rm --volume $PWD:/workspace --volume $PWD/data/params/unimod.xml:/opt/OpenMS/share/OpenMS/CHEMISTRY/unimod.xml ghcr.io/radusuciu/yan_etal_natcomm_2024 \
        IDFilter \
            -in "/workspace/{input}" \
            -out "/workspace/{output}" \
            -score:pep 0.01 \
            -remove_decoys \
            -delete_unreferenced_peptide_hits
        """

rule symlink_mzml:
    input:
        "data/pipeline/{sample}/{sample}.mzML"
    output:
        "data/input/{sample}.mzML"
    shell:
        """
        ln -sf $(realpath {input}) {output}
        """

rule copy_results:
    input:
        "data/pipeline/{sample}/filtered_matches.idXML"
    output:
        "data/input/{sample}_filtered_matches.idXML"
    shell:
        """
        cp {input} {output}
        """
