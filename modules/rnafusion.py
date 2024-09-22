"""Module for running nf-core/rnafusion."""

from json import loads
from logging import LoggerAdapter
from pathlib import Path

from cellophane import Checkpoints, Config, Executor, Samples, output, runner
from mpire.async_result import AsyncResult
from ruamel.yaml import YAML

from modules.common import nf_config
from modules.nextflow import nextflow

# Taken from https://github.com/nf-core/rnafusion/blob/3.0.1/conf/modules.config
# The whole string needs to be included here only to override the limitSjdbInsertNsj
# parameter.
rnafusion_nf_config = """\
process {{
    withName: 'STAR_FOR_STARFUSION' {{
        ext.args = '--twopassMode Basic \
        --outReadsUnmapped None \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --outSAMunmapped Within \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 8 \
        --chimOutJunctionFormat 1 \
        --alignSJDBoverhangMin 10 \
        --alignMatesGapMax 100000 \
        --alignIntronMax 100000 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimMultimapScoreRange 3 \
        --chimScoreJunctionNonGTAG -4 \
        --chimMultimapNmax 20 \
        --chimNonchimScoreDropMin 10 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --alignInsertionFlush Right \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 \
        --chimOutType Junctions \
        --quantMode GeneCounts \
        --limitSjdbInsertNsj {config.limitSjdbInsertNsj}'
    }}
    withName: 'STAR_FOR_ARRIBA' {{
        ext.args = '--readFilesCommand zcat \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outBAMcompression 0 \
        --outFilterMultimapNmax 50 \
        --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 \
        --chimOutType WithinBAM HardClip \
        --chimJunctionOverhangMin 10 \
        --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3 \
        --chimMultimapNmax 50 \
        --limitSjdbInsertNsj {config.limitSjdbInsertNsj}'
    }}
}}
"""


def _patch_fusionreport(
    samples: Samples,
    workdir: Path,
    logger: LoggerAdapter,
    log_tag: str,
):
    """
    Patch fusionreport html to keep fusion data separate from the report itself.

    This is done to make it easier to access the report from the directory listing.
    """
    logger.info(f"Patching fusionreport html ({log_tag})")
    for id_, group in samples.split(by="id"):
        original = workdir / f"fusionreport/{id_}/{id_}_fusionreport_index.html"
        patched = workdir / f"{id_}.fusionreport.html"

        try:
            index = original.read_text(encoding="utf-8").replace(
                "${fusion.replace('--','_')}.html",
                "fusionreport/${fusion.replace('--','_')}.html",
            )

            patched.write_text(index, encoding="utf-8")
        except Exception as exception:  # pylint: disable=broad-except
            logger.error(f"Failed to patch fusionreport for {id_}: {exception}")
            for sample in group:
                sample.fail(f"Failed to patch fusionreport: {exception}")


def _validate_inputs(config: Config, root: Path, logger: LoggerAdapter):
    if not config.rnafusion.genomes_base.exists():
        logger.error("Missing required reference files for nf-core/rnafusion.")
        logger.error(
            "Instructions for downloading reference files can be found at "
            "https://nf-co.re/rnafusion/3.0.1/docs/usage#download-and-build-references"
        )
        logger.error(
            f"Use {root}/dependencies/nf-core/rnafusion/main.nf when downloading "
            "the reference files"
        )
        raise SystemExit(1)


def _pipeline_args(config: Config, workdir: Path, nf_samples: Path, /):
    return [
        f"--outdir {workdir}",
        f"--input {nf_samples}",
        f"--genomes_base {config.rnafusion.genomes_base}",
        f"--read_length {config.read_length}",
        f"--tools_cutoff {config.rnafusion.tools_cutoff}",
        f"--fusioncatcher_limitSjdbInsertNsj {config.limitSjdbInsertNsj}",
        f"--fusioninspector_limitSjdbInsertNsj {config.limitSjdbInsertNsj}",
        "--all",
    ]


def _standalone_arriba_visualisation(
    samples: Samples,
    config: Config,
    executor: Executor,
    workdir: Path,
    root: Path,
    logger: LoggerAdapter,
    **kwargs,
) -> None:
    logger.info("Generating standalone arriba visualisations")

    try:
        params_glob = (workdir / "pipeline_info").glob("params_*.json")
        versions_glob = (workdir / "pipeline_info").glob("software_versions.yml")
        params_json = next(params_glob)
        versions_yml = next(versions_glob)
        params = loads(params_json.read_text())
        versions = YAML(typ="safe").load(versions_yml)
    except Exception:
        params = {}
        versions = {}

    cytobands, protein_domains, annotation, version = (
        config.rnafusion.arriba_standalone.get("cytobands")
        or params.get("arriba_ref_cytobands"),
        config.rnafusion.arriba_standalone.get("protein_domains")
        or params.get("arriba_ref_protein_domains"),
        config.rnafusion.arriba_standalone.get("annotation") or params.get("gtf"),
        config.rnafusion.arriba_standalone.get("version")
        or versions.get("ARRIBA_VISUALISATION", {}).get("arriba"),
    )

    if not (cytobands and protein_domains and annotation and version):
        logger.warning("Unable to detect params for standalone arriba visualisation")
        return

    results: list[AsyncResult] = []

    for _id in samples.unique_ids:
        if not (workdir / "arriba" / f"{_id}.arriba.fusions.tsv").exists():
            logger.warning(f"No arriba fusions found for sample {_id}")
            continue

        result, _ = executor.submit(
            str(root / "scripts" / "arriba_draw_fusions.sh"),
            env={
                "_ARRIBA_STANDALONE_FUSIONS": f"{workdir}/arriba/{_id}.arriba.fusions.tsv",
                "_ARRIBA_STANDALONE_BAM": f"{workdir}/star_for_arriba/{_id}.Aligned.out.bam",
                "_ARRIBA_STANDALONE_OUTPUT": f"{workdir}/arriba_visualisation/{_id}_standalone_arriba_visualisation.pdf",
                "_ARRIBA_STANDALONE_CYTOBANDS": str(cytobands),
                "_ARRIBA_STANDALONE_PROTEIN_DOMAINS": str(protein_domains),
                "_ARRIBA_STANDALONE_ANNOTATION": str(annotation),
                "_ARRIBA_STANDALONE_THREADS": config.rnafusion.arriba_standalone.threads,
            },
            workdir=workdir,
            name="standalone_arriba_visualisation",
            conda_spec={"dependencies": [f"arriba ={version}", "samtools =1.16"]},
            cpus=config.rnafusion.arriba_standalone.threads,
            **kwargs,
        )
        results.append(result)

    for result in results:
        result.get()


@output(
    "arriba_visualisation/{sample.id}_combined_fusions_arriba_visualisation.pdf",
    dst_dir="{sample.id}_{sample.last_run}",
    checkpoint="DUMMY",
)
@output(
    "arriba_visualisation/{sample.id}_standalone_arriba_visualisation.pdf",
    dst_dir="{sample.id}_{sample.last_run}/arriba",
)
@output(
    "arriba/{sample.id}.*",
    dst_dir="{sample.id}_{sample.last_run}/arriba",
)
@output(
    "fusioncatcher/{sample.id}.*",
    dst_dir="{sample.id}_{sample.last_run}/fusioncatcher",
)
@output(
    "starfusion/{sample.id}.*",
    dst_dir="{sample.id}_{sample.last_run}/starfusion",
)
@output(
    "fusionreport/{sample.id}",
    dst_name="{sample.id}_{sample.last_run}/fusionreport",
)
@output(
    "{sample.id}.fusionreport.html",
    dst_dir="{sample.id}_{sample.last_run}",
)
@output(
    "fusioninspector/{sample.id}.*",
    dst_dir="{sample.id}_{sample.last_run}/fusioninspector",
)
@output(
    "star_for_starfusion/{sample.id}.Aligned.sortedByCoord.out.bam",
    dst_dir="{sample.id}_{sample.last_run}",
)
@output(
    "star_for_starfusion/{sample.id}.Aligned.sortedByCoord.out.bam.bai",
    dst_dir="{sample.id}_{sample.last_run}",
)
@output(
    "pipeline_info",
    dst_name="{sample.id}_{sample.last_run}/pipeline_info/rnafusion",
)
@runner(split_by="id")
def rnafusion(
    samples: Samples,
    config: Config,
    workdir: Path,
    executor: Executor,
    logger: LoggerAdapter,
    root: Path,
    checkpoints: Checkpoints,
    **_,
) -> Samples:
    """Run nf-core/rnafusion."""
    log_tag = samples[0].id if (n := len(samples.unique_ids)) == 1 else f"{n} samples"
    if config.rnafusion.skip:
        samples.output = set()
        return samples

    if checkpoints.main.check(rnafusion_nf_config):
        logger.info(f"Using previous nf-core/rnafusion output ({log_tag})")
        return samples

    _validate_inputs(
        config=config,
        root=root,
        logger=logger,
    )

    nf_config(
        template=rnafusion_nf_config,
        location=workdir / "nextflow.config",
        include=config.nextflow.get("config"),
        config=config,
    )

    logger.info(f"Running nf-core/rnafusion ({log_tag})")

    sample_sheet = samples.nfcore_samplesheet(
        location=workdir,
        strandedness=config.strandedness,
    )

    nextflow(
        root / "dependencies" / "nf-core" / "rnafusion" / "main.nf",
        *_pipeline_args(config, workdir, sample_sheet),
        nxf_config=workdir / "nextflow.config",
        config=config,
        name="rnafusion",
        workdir=workdir,
        resume=True,
        executor=executor,
    )

    logger.debug(f"nf-core/rnafusion finished ({log_tag})")

    _patch_fusionreport(
        samples=samples,
        workdir=workdir,
        logger=logger,
        log_tag=log_tag,
    )

    _standalone_arriba_visualisation(
        samples=samples,
        workdir=workdir,
        config=config,
        executor=executor,
        root=root,
        logger=logger,
    )

    checkpoints.main.store(rnafusion_nf_config)

    return samples
