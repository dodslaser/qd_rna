"""Module for fetching files from HCP."""

from pathlib import Path
from typing import Optional

from cellophane import cfg, sge


def nextflow(
    main: Path,
    *args,
    config: cfg.Config,
    env: dict[str, str] = {},
    log: Optional[Path] = None,
    report: Optional[Path] = None,
    workdir: Optional[Path] = None,
    nf_config: Optional[Path] = None,
    ansi_log: bool = False,
    resume: bool = False,
    **kwargs,
):
    if "workdir" in config.nextflow:
        config.nextflow.workdir.mkdir(parents=True, exist_ok=True)

    sge.submit(
        str(Path(__file__).parent / "scripts" / "nextflow.sh"),
        f"-log {log}",
        (
            f"-config {config.nextflow.config}"
            if "config" in config.nextflow
            else nf_config or ""
        ),
        f"run {main}",
        "-ansi-log false" if not ansi_log or config.nextflow.ansi_log else "",
        f"-work-dir {workdir}" if workdir else "",
        "-resume" if resume else "",
        f"-with-report {report}" if report else "",
        f"-profile {config.nextflow.profile}",
        *args,
        env={
            "_NEXTFLOW_MODULE": config.nextflow.nf_module,
            "_JAVA_MODULE": config.nextflow.java_module,
            **env,
        },
        queue=config.nextflow.sge_queue,
        pe=config.nextflow.sge_pe,
        slots=config.nextflow.sge_slots,
        **kwargs,
    )
