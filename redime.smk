SAMPLES = [
    "brain_conc1.0_rep1",
    "brain_conc1.0_rep2",
    "brain_conc1.0_rep3",
    "brain_conc10.0_rep1",
    "brain_conc10.0_rep2",
    "brain_conc10.0_rep3",
    "kidney_conc1.0_rep1",
    "kidney_conc1.0_rep2",
    "kidney_conc1.0_rep3",
    "kidney_conc10.0_rep1",
    "kidney_conc10.0_rep2",
    "kidney_conc10.0_rep3",
]

rule all:
    input:
        expand('data/input/{sample}.raw_output_rt_10_sn_2.5.to_excel', sample=SAMPLES)

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
            -P/workspace/data/params/comet_redime.params \
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

rule export_protein_texts:
    input:
        "data/pipeline/{sample}/filtered_matches.idXML"
    output:
        "data/pipeline/{sample}/proteins.csv"
    shell:
        """
        docker run -it --rm --volume $PWD:/workspace --volume $PWD/data/params/unimod.xml:/opt/OpenMS/share/OpenMS/CHEMISTRY/unimod.xml ghcr.io/radusuciu/yan_etal_natcomm_2024 \
        TextExporter \
            -in "/workspace/{input}" \
            -out "/workspace/{output}" \
            -id:proteins_only \
            -quoting "double" \
            -feature:add_metavalues 100 \
            -id:add_metavalues 100 \
            -id:add_hit_metavalues 100
        """

rule export_psms_texts:
    input:
        "data/pipeline/{sample}/filtered_matches.idXML"
    output:
        "data/pipeline/{sample}/psms.csv"
    shell:
        """
        docker run -it --rm --volume $PWD:/workspace --volume $PWD/data/params/unimod.xml:/opt/OpenMS/share/OpenMS/CHEMISTRY/unimod.xml ghcr.io/radusuciu/yan_etal_natcomm_2024 \
        TextExporter \
            -in "/workspace/{input}" \
            -out "/workspace/{output}" \
            -id:peptides_only \
            -quoting "double" \
            -feature:add_metavalues 100 \
            -id:add_metavalues 100 \
            -id:add_hit_metavalues 100
        """

rule prepare_cimage_inputs:
    input:
        "data/pipeline/{sample}/psms.csv",
        "data/pipeline/{sample}/proteins.csv",
        "data/pipeline/{sample}/{sample}.mzML",
    output:
        "data/pipeline/{sample}/cimage/dta/ipi_name.table",
        "data/pipeline/{sample}/cimage/dta/all_scan.table",
        "data/pipeline/{sample}/cimage/dta/cross_scan.table",
        "data/pipeline/{sample}/cimage/{sample}_01.mzML",
    run:
        import pathlib
        from utils import create_input_files, setup_cimage_folder

        setup_cimage_folder(pathlib.Path(f'data/pipeline/{wildcards.sample}'))
        print(wildcards.sample)
        create_input_files(
            experiment_name=wildcards.sample,
            cimage_folder=pathlib.Path(f'data/pipeline/{wildcards.sample}/cimage'),
        )

rule prepare_cimage_params:
    output:
        "data/pipeline/{sample}/cimage/cimage.params",
        "data/pipeline/{sample}/cimage/cimage_heavy.table",
        "data/pipeline/{sample}/cimage/cimage_light.table",
    shell:
        """
        for f in {output} ; do
            source_path=data/params/$(basename $f)
            cp $source_path $f
        done

        sed -i "s/cimage_light.table/\/workspace\/data\/pipeline\/{wildcards.sample}\/cimage\/cimage_light.table/g" data/pipeline/{wildcards.sample}/cimage/cimage.params
        sed -i "s/cimage_heavy.table/\/workspace\/data\/pipeline\/{wildcards.sample}\/cimage\/cimage_heavy.table/g" data/pipeline/{wildcards.sample}/cimage/cimage.params
        """

rule cimage:
    input:
        "data/pipeline/{sample}/cimage/cimage.params",
        "data/pipeline/{sample}/cimage/cimage_heavy.table",
        "data/pipeline/{sample}/cimage/cimage_light.table",
        "data/pipeline/{sample}/cimage/dta/ipi_name.table",
        "data/pipeline/{sample}/cimage/dta/all_scan.table",
        "data/pipeline/{sample}/cimage/dta/cross_scan.table",
        "data/pipeline/{sample}/cimage/{sample}_01.mzML",
    output:
        "data/pipeline/{sample}/cimage/dta/output/output_rt_10_sn_2.5.to_excel.txt"
    shell:
        """
        docker run -it --rm --volume $PWD:/workspace --workdir /workspace/data/pipeline/{wildcards.sample}/cimage/dta ghcr.io/radusuciu/yan_etal_natcomm_2024 \
            /opt/cimage/cimageFromOpenMS.R ../cimage.params {wildcards.sample}
        """

rule copy_results:
    input:
        "data/pipeline/{sample}/cimage/dta/output/output_rt_10_sn_2.5.to_excel.txt"
    output:
        "data/input/{sample}.raw_output_rt_10_sn_2.5.to_excel"
    shell:
        """
        cp {input} {output}
        """
