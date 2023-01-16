"""Module for running nf-core/rnafusion."""

import time
from logging import LoggerAdapter
from pathlib import Path

from cellophane import cfg, slims, modules, sge

def get_output(outdir: Path, config: cfg.Config):
    """Return a dictionary of output files for the rnafusion module."""
    return {
        "arriba": [*(outdir / config.rnaseq.aligner).glob("*.tsv")],
        "arrriba_visualisation": [*(outdir / "arriba_visualisation").glob("*.pdf")],
        "fusioncatcher": [*(outdir / "fusioncatcher").glob("*.txt")],
        "fusioninspector": [*(outdir / "fusioninspector").glob("*")],
        "fusionreport": [*(outdir / "fusionreport").glob("**/*.html")],
        "pizzly": [*(outdir / "pizzly").glob("*")],
        "squid": [*(outdir / "squid").glob("*.txt")],
        "starfusion": [*(outdir / "starfusion").glob("*.tsv")],
    }

@modules.runner()
def rnafusion(
    label: str,
    samples: slims.Samples,
    config: cfg.Config,
    logger: LoggerAdapter,
    scripts_path: Path,
) -> None:
    """Run nf-core/rnafusion."""
    timestamp = time.strftime("%y%m%d-%H%M%S")
    outdir = config.outdir / timestamp / label
    

    if not config.rnafusion.skip:
        logger.info("Running nf-core/rnafusion")
        logger.debug(f"Output will be written to {outdir}")

        return

        outdir.mkdir(parents=True, exist_ok=True)
        sample_sheet = samples.nfcore_samplesheet(
            location=outdir,
            strandedness=config.rnafusion.strandedness,
        )
        
        if "workdir" in config.nextflow:
            config.nextflow.workdir.mkdir(parents=True, exist_ok=True)
        sge.submit(
            str(scripts_path / "nextflow.sh"),
            f"-log {outdir / 'logs' / 'rnafusion.log'}",
            (
                f"-config {config.nextflow.config}"
                if "config" in config.nextflow
                else ""
            ),
            f"run {config.rnafusion.nf_main}",
            "-ansi-log false",
            f"-work-dir {config.nextflow.workdir or outdir / 'work'}",
            ("-resume" if config.nextflow.resume else ""),
            f"-with-report {outdir / 'logs' / 'rnafusion-execution.html'}",
            f"-profile {config.nextflow.profile}",
            f"--outdir {outdir}",
            f"--input {sample_sheet}",
            f"--genomes_base {config.rnafusion.genomes_base}",
            f"--arriba_ref {config.rnafusion.arriba_ref}",
            f"--arriba_ref_blacklist {config.rnafusion.arriba_blacklist}",
            f"--arriba_ref_protein_domain {config.rnafusion.arriba_protein_domain}",
            f"--fusionreport_tool_cutoff {config.rnafusion.fusionreport_tool_cutoff}",
            f"--read_length {config.rnafusion.read_length}",
            "--all",
            "--fusioninspector_filter",
            env={"_MODULES_INIT": config.modules_init},
            queue=config.nextflow.sge_queue,
            pe=config.nextflow.sge_pe,
            slots=config.nextflow.sge_slots,
            check=True,
            name="rnafusion",
            stderr=config.logdir / f"rnafusion.{timestamp}.err",
            stdout=config.logdir / f"rnafusion.{timestamp}.out",
            cwd=outdir,
        )

        output = get_output(outdir, config)

        logger.debug(output)