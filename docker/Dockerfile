# CPU-only image with the real scientific stack (RDKit, fpocket, Vina, OpenBabel).
# GPU tools (BioEmu, GNINA, REINVENT4) are intentionally excluded -- see
# docs/DEPENDENCIES.md for the isolated-GPU-worker strategy.
FROM continuumio/miniconda3:latest

WORKDIR /app

COPY environment.yml .
RUN conda env create -f environment.yml && conda clean -afy

SHELL ["conda", "run", "-n", "conform", "/bin/bash", "-c"]
ENV PATH=/opt/conda/envs/conform/bin:$PATH
ENV PYTHONPATH=/app

COPY . .

RUN python -c "import rdkit, vina, MDAnalysis; print('scientific stack ok')" \
 && fpocket -h >/dev/null 2>&1 || true

CMD ["python", "scripts/run_experiment.py", "validate"]
