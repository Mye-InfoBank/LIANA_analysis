FROM continuumio/miniconda3:24.3.0-0

WORKDIR /app

COPY environment.yml /app/environment.yml
RUN conda env create -f /app/environment.yml && conda clean -afy

COPY shiny /app/shiny
COPY results_IBD /app/results_IBD

ENV PORT=8080
ENV DATA_DIR=/app/results_IBD/colon

WORKDIR /app/shiny

EXPOSE 8080

CMD ["bash", "-lc", "conda run --no-capture-output -n liana_dashboard python app.py --host 0.0.0.0 --port ${PORT} ${DATA_DIR}"]